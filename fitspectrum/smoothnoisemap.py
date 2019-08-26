#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Make noise maps from the variance maps, smooth them, and then work out the new variance map
#
# Mike Peel    03 Sep 2017    Start
# Mike Peel    07 Sep 2017    Bug fixes / tidying running order
# Mike Peel    17 Sep 2017    ud_grade to use constant nsides
# Mike Peel    04 Oct 2017    Optimising and debugging. Added multiple nside support.
# Mike Peel    06 Oct 2017    Return to older variance calculation, add rescale param and change output name format
# Mike Peel    05 Jun 2019    v0.4 Add taper, cope with cut sky maps

import numpy as np
import healpy as hp
from astrocode.fitspectrum.smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits
import os
from scipy import optimize

def gaussfit(x, param):
    return hp.gauss_beam(np.radians(param/60.0),300)


def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
    return newmap


def smoothnoisemap(indir, outdir, runname, inputmap, mapnumber=2, fwhm=0.0, numrealisations=10, sigma_0 = 0.0, nside=[512], windowfunction = [], rescale=1.0,usehealpixfits=False,taper=False,lmin_taper=350,lmax_taper=600,taper_gauss=False):
    ver = "0.4"

    if (os.path.isfile(indir+"/"+runname+"_actualvariance.fits")):
        print("You already have a file with the output name " + indir+"/"+runname+"_actualvariance.fits" + "! Not going to overwrite it. Move it, or set a new output filename, and try again!")
        # exit()
        return

    # Read in the input map
    inputfits = fits.open(indir+"/"+inputmap)
    cols = inputfits[1].columns
    col_names = cols.names
    nmaps = len(cols)
    maps = []
    if usehealpixfits:
        maps = hp.read_map(indir+inputmap,field=None)
    else:
        for i in range(0,nmaps):
            maps.append(inputfits[1].data.field(i))
    nside_in = hp.get_nside(maps)

    # Check to see whether we have nested data, and switch to ring if that is the case.
    if (inputfits[1].header['ORDERING'] == 'NESTED'):
        maps = hp.reorder(maps,n2r=True)

    maps[mapnumber][maps[mapnumber]<-1e10] = hp.UNSEEN
    maps[mapnumber][~np.isfinite(maps[mapnumber])] = hp.UNSEEN
    map_before = maps[mapnumber].copy()

    if sigma_0 != 0.0:
        # If we have a value for sigma_0, then we have an Nobs map and need to convert it.
        maps[mapnumber][maps[mapnumber]<=0.0] = hp.UNSEEN
        maps[mapnumber] = conv_nobs_variance_map(maps[mapnumber], sigma_0)

    # We want to sqrt it to get a noise rms map
    noisemap = np.sqrt(maps[mapnumber])
    noisemap = noisemap * rescale

    # Write the variance map to disk so we can compare to it later.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=np.square(noisemap)))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.BinTableHDU.from_columns(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['NSIDE']=nside_in
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(outdir+"/"+runname+"_actualvariance.fits")

    # Also save the input nobs map, as a cross-check.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=conv_nobs_variance_map(np.square(noisemap), sigma_0)))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.BinTableHDU.from_columns(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['NSIDE']=nside_in
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(outdir+"/"+runname+"_actualnobs.fits")

    # Calculate the window function
    conv_windowfunction = hp.gauss_beam(np.radians(fwhm/60.0),4*nside_in)
    if (windowfunction != []):
        window_len = len(conv_windowfunction)
        beam_len = len(windowfunction)
        if (beam_len > window_len):
            windowfunction  = windowfunction[0:len(conv_windowfunction)]
        else:
            windowfunction = np.pad(windowfunction, (0, window_len - beam_len), 'constant')
        conv_windowfunction[windowfunction!=0] /= windowfunction[windowfunction!=0]
        conv_windowfunction[windowfunction==0] = 0.0
    conv_windowfunction /= conv_windowfunction[0]
    # If needed, apply a taper
    if taper:
        conv_windowfunction[lmin_taper:lmax_taper] = conv_windowfunction[lmin_taper:lmax_taper] * np.cos((np.pi/2.0)*((np.arange(lmin_taper,lmax_taper)-lmin_taper)/(lmax_taper-lmin_taper)))
        conv_windowfunction[lmax_taper:] = 0.0
    if taper_gauss:
        trip = 0
        val = 0
        beam1 = hp.gauss_beam(np.radians(fwhm/60.0),len(conv_windowfunction))
        param_est, cov_x = optimize.curve_fit(gaussfit, range(0,299), windowfunction[0:301], 60.0)
        beam2 = hp.gauss_beam(np.radians(param_est[0]/60.0),len(conv_windowfunction))
        # plt.plot(windowfunction[0:301])
        # plt.plot(beam2)
        # plt.savefig(outdir+'temp.pdf')
        print(param_est[0])
        # exit()
        for l in range(1,len(conv_windowfunction)):
            if trip == 1:
                # conv_windowfunction[l] = val * np.exp(-0.5*(np.radians(fwhm_arcmin/60.0)**2-np.radians(param_est[0]/60.0)**2)*l*(l+1))
                conv_windowfunction[l] = val * (beam1[l]/beam2[l])
            elif (conv_windowfunction[l]-conv_windowfunction[l-1]) > 0.0:
                print(l)
                trip = 1
                # val = conv_windowfunction[l]/np.exp(-0.5*(np.radians(fwhm_arcmin/60.0)**2-np.radians(param_est[0]/60.0)**2)*l*(l+1))
                val = conv_windowfunction[l-1]/(beam1[l-1]/beam2[l-1])
                conv_windowfunction[l] = val * (beam1[l]/beam2[l])

    # Now generate the noise realisations
    numpixels = len(noisemap)
    returnmap = np.zeros(numpixels)
    # print(np.min(noisemap))
    noisemap[map_before == hp.UNSEEN] = 0.0
    hp.write_map(outdir+"/"+runname+"_noisemap.fits",noisemap,overwrite=True)
    # alms = hp.map2alm(noisemap)#,lmax=4*nside_in)
    # alms = hp.almxfl(alms, conv_windowfunction)
    # newnoisemap = hp.alm2map(alms, nside_in)#,lmax=4*nside_in)
    # hp.write_map(outdir+"/"+runname+"_newnoisemap.fits",newnoisemap,overwrite=True)
    # exit()
    for i in range(0,numrealisations):
        if i%10==0:
            print(i)
        # Generate the noise realisation
        newmap = noiserealisation(noisemap, numpixels)
        # smooth it
        alms = hp.map2alm(newmap)#,lmax=4*nside_in)
        alms = hp.almxfl(alms, conv_windowfunction)
        newmap = hp.alm2map(alms, nside_in)#,lmax=4*nside_in)
        returnmap = returnmap + np.square(newmap)
    returnmap = returnmap/(numrealisations-1)
    returnmap[map_before == hp.UNSEEN] = hp.UNSEEN
    
    # All done - now just need to write it to disk.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=returnmap))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.BinTableHDU.from_columns(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['NSIDE']=nside_in
    bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's smoothnoisemap version "+ver +"."
    bin_hdu.writeto(outdir+"/"+runname+"_variance.fits")

    # Also do an Nobs map for a consistency check.
    nobs_map = conv_nobs_variance_map(returnmap, sigma_0)
    cols = []
    cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.BinTableHDU.from_columns(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['NSIDE']=nside_in
    bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
    bin_hdu.writeto(outdir+"/"+runname+"_nobs.fits")

    # Do ud_graded versions
    num_nside = len(nside)
    for i in range(0,num_nside):
        # ud_grade it using power=0 (assuming correlated pixels)
        returnmap_ud = hp.ud_grade(returnmap, nside[i], power=0)

        # Output the variance map
        cols = []
        cols.append(fits.Column(name='II_cov', format='E', array=returnmap_ud))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.BinTableHDU.from_columns(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        bin_hdu.header['NSIDE']=nside[i]
        bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
        bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_variance.fits")

        # Also do an Nobs map for a consistency check.
        nobs_map = conv_nobs_variance_map(returnmap_ud, sigma_0)
        cols = []
        cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.BinTableHDU.from_columns(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        bin_hdu.header['NSIDE']=nside[i]
        bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
        bin_hdu.writeto(outdir+"/"+str(nside[i])+"_"+runname+"_nobs.fits")



    return