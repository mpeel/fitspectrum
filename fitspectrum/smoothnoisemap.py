#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Make noise maps from the variance maps, smooth them, and then work out the new variance map
#
# Mike Peel    03 Sep 2017    Start
# Mike Peel    07 Sep 2017    Bug fixes / tidying running order
# Mike Peel    17 Sep 2017    ud_grade to use constant nsides
# Mike Peel    04 Oct 2017    Optimising and debugging. Added multiple nside support.
# Mike Peel    06 Oct 2017    Return to older variance calculation, add rescale param and change output name format

import numpy as np
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
    return newmap


def smoothnoisemap(indir, runname, inputmap, mapnumber=2, fwhm=0.0, numrealisations=10, sigma_0 = 0.0, nside=[512], windowfunction = [], rescale=1.0):
    ver = "0.3"

    # Read in the input map
    inputfits = fits.open(indir+"/"+inputmap)
    cols = inputfits[1].columns
    col_names = cols.names
    nmaps = len(cols)
    maps = []
    for i in range(0,nmaps):
        maps.append(inputfits[1].data.field(i))
    nside_in = hp.get_nside(maps)

    # Check to see whether we have nested data, and switch to ring if that is the case.
    if (inputfits[1].header['ORDERING'] == 'NESTED'):
        maps = hp.reorder(maps,n2r=True)

    if sigma_0 != 0.0:
        # If we have a value for sigma_0, then we have an Nobs map and need to convert it.
        maps[mapnumber] = conv_nobs_variance_map(maps[mapnumber], sigma_0)

    # We want to sqrt it to get a noise rms map
    noisemap = np.sqrt(maps[mapnumber])
    noisemap = noisemap * rescale

    # Write the variance map to disk so we can compare to it later.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=np.square(noisemap)))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(indir+"/"+runname+"_actualvariance.fits")

    # Also save the input nobs map, as a cross-check.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=conv_nobs_variance_map(np.square(noisemap), sigma_0)))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(indir+"/"+runname+"_actualnobs.fits")

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

    # Now generate the noise realisations
    numpixels = len(noisemap)
    returnmap = np.zeros(numpixels)
    for i in range(0,numrealisations):
        if i%10==0:
            print i
        # Generate the noise realisation
        newmap = noiserealisation(noisemap, numpixels)
        # smooth it
        alms = hp.map2alm(newmap,lmax=4*nside_in)
        alms = hp.almxfl(alms, conv_windowfunction)
        newmap = hp.alm2map(alms, nside_in,lmax=4*nside_in)
        returnmap = returnmap + np.square(newmap)
    returnmap = returnmap/(numrealisations-1)
    
    # All done - now just need to write it to disk.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=returnmap))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's smoothnoisemap version "+ver +"."
    bin_hdu.writeto(indir+"/"+runname+"_variance.fits")

    # Also do an Nobs map for a consistency check.
    nobs_map = conv_nobs_variance_map(returnmap, sigma_0)
    cols = []
    cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
    bin_hdu.writeto(indir+"/"+runname+"_nobs.fits")

    # Do ud_graded versions
    num_nside = len(nside)
    for i in range(0,num_nside):
        # ud_grade it using power=0 (assuming correlated pixels)
        returnmap_ud = hp.ud_grade(returnmap, nside[i], power=0)

        # Output the variance map
        cols = []
        cols.append(fits.Column(name='II_cov', format='E', array=returnmap_ud))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.new_table(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
        bin_hdu.writeto(indir+"/"+str(nside[i])+"_"+runname+"_variance.fits")

        # Also do an Nobs map for a consistency check.
        nobs_map = conv_nobs_variance_map(returnmap_ud, sigma_0)
        cols = []
        cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.new_table(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's smoothnoisemap version "+ver +"."
        bin_hdu.writeto(indir+"/"+str(nside[i])+"_"+runname+"_nobs.fits")



    return