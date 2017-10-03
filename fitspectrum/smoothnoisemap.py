#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Make noise maps from the variance maps, smooth them, and then work out the new variance map
#
# Mike Peel    03 Sep 2017    Start
# Mike Peel    07 Sep 2017    Bug fixes / tidying running order
# Mike Peel    17 Sep 2017    ud_grade to use constant nsides

import numpy as np
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
#    for i in range(0,numpixels):
#        newmap[i] = np.random.normal(scale=inputmap[i])
    return newmap


def smoothnoisemap(indir, runname, inputmap, mapnumber=2, fwhm=0.0, numrealisations=10, sigma_0 = 0.0, nside=[512]):
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

    # Write the variance map to disk so we can compare to it later.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=np.square(noisemap)))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    # bin_hdu.header['NSIDE']=nside_out
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
    # bin_hdu.header['NSIDE']=nside_out
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(indir+"/"+runname+"_actualnobs.fits")

    numpixels = len(noisemap)
    # nside = hp.get_nside(maps)

    returnmap = np.zeros(numpixels)
    conv_windowfunction = hp.gauss_beam(np.radians(fwhm/60.0),3*nside_in-1)
    print conv_windowfunction[0]
    conv_windowfunction /= conv_windowfunction[0]
    print conv_windowfunction
    print fwhm
    for i in range(0,numrealisations):
        print i
        # Generate the noise realisation
        newmap = noiserealisation(noisemap, numpixels)
        # smooth it
        # newmap = hp.smoothing(newmap, fwhm=np.radians(fwhm/60.0))
        alms = hp.map2alm(newmap)
        # alms = hp.almxfl(alms, conv_windowfunction)
        newmap2 = hp.alm2map(alms, nside_in,verbose=False)
        returnmap = returnmap + np.square(newmap2)
        print np.median(returnmap/i)


    returnmap = returnmap/(numrealisations-1)
    print np.median(returnmap)
    print np.median(np.square(noisemap))
    print np.median(np.square(noisemap))/np.median(returnmap)

    # All done - now just need to write it to disk.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=returnmap))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    # bin_hdu.header['NSIDE']=nside_out
    bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's code for testing purposes only."
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
    # bin_hdu.header['NSIDE']=nside_out
    bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's code for testing purposes only."
    bin_hdu.writeto(indir+"/"+runname+"_nobs.fits")

    # Do ud_graded versions
    num_nside = len(nside)
    for i in range(0,num_nside):
        # ud_grade it
        returnmap_ud = hp.ud_grade(returnmap, nside[i], power=0)

        cols = []
        cols.append(fits.Column(name='II_cov', format='E', array=returnmap_ud))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.new_table(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        # bin_hdu.header['NSIDE']=nside_out
        bin_hdu.header['COMMENT']="Smoothed variance map calculated by Mike Peel's code for testing purposes only."
        bin_hdu.writeto(indir+"/"+runname+"_"+str(nside[i])+"_variance.fits")

        # Also do an Nobs map for a consistency check.
        nobs_map = conv_nobs_variance_map(returnmap_ud, sigma_0)
        cols = []
        cols.append(fits.Column(name='II_nobs', format='E', array=nobs_map))
        cols = fits.ColDefs(cols)
        bin_hdu = fits.new_table(cols)
        bin_hdu.header['ORDERING']='RING'
        bin_hdu.header['POLCONV']='COSMO'
        bin_hdu.header['PIXTYPE']='HEALPIX'
        # bin_hdu.header['NSIDE']=nside_out
        bin_hdu.header['COMMENT']="Smoothed Nobs map calculated by Mike Peel's code for testing purposes only."
        bin_hdu.writeto(indir+"/"+runname+"_"+str(nside[i])+"_nobs.fits")



    return