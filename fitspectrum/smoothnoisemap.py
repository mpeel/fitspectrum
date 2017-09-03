#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Make noise maps from the variance maps, smooth them, and then work out the new variance map
#
# Mike Peel    03 Sep 2017    Start


import numpy as np
import healpy as hp
from smoothmap import smoothmap, conv_nobs_variance_map
import astropy.io.fits as fits

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    test = 0
    for i in range(0,numpixels):
        newmap[i] = np.random.normal(scale=inputmap[i])
        test += 1
    # print test, numpixels, len(newmap)
    return newmap


def smoothnoisemap(indir, runname, inputmap, mapnumber=2, fwhm=0.0, numrealisations=10, sigma_0 = 0.0):
    # Read in the input map
    inputfits = fits.open(indir+"/"+inputmap)
    cols = inputfits[1].columns
    col_names = cols.names
    nmaps = len(cols)
    maps = []
    for i in range(0,nmaps):
        maps.append(inputfits[1].data.field(i))
    # Check to see whether we have nested data, and switch to ring if that is the case.
    if (inputfits[1].header['ORDERING'] == 'NESTED'):
        maps = hp.reorder(maps,n2r=True)

    # If we have a value for sigma_0, then we have an Nobs map and need to convert it.
    if sigma_0 != 0.0:
        maps[mapnumber] = np.sqrt(conv_nobs_variance_map(maps[mapnumber], sigma_0))

    # Write the variance map to disk so we can compare to it later.
    cols = []
    cols.append(fits.Column(name='II_cov', format='E', array=maps[mapnumber]))
    cols = fits.ColDefs(cols)
    bin_hdu = fits.new_table(cols)
    bin_hdu.header['ORDERING']='RING'
    bin_hdu.header['POLCONV']='COSMO'
    bin_hdu.header['PIXTYPE']='HEALPIX'
    # bin_hdu.header['NSIDE']=nside_out
    bin_hdu.header['COMMENT']="Input variance map - for test purposes only."
    bin_hdu.writeto(indir+"/"+runname+"_actualvariance.fits")

    # We want to work with the std not the variance
    maps[mapnumber] = np.sqrt(maps[mapnumber])

    numpixels = len(maps[mapnumber])
    nside = hp.get_nside(maps)
    # print 'Number of pixels should be: '
    # print hp.nside2npix(nside)

    returnmap = np.zeros(numpixels)
    for i in range(0,numrealisations):
        print i
        # Generate the noise realisation
        newmap = noiserealisation(maps[mapnumber], numpixels)
        # smooth it
        conv_windowfunction = hp.gauss_beam(np.radians(fwhm/60.0),3*nside)
        # newmap = hp.smoothing(newmap, fwhm=np.radians(fwhm))
        alms = hp.map2alm(newmap)
        alms = hp.almxfl(alms, conv_windowfunction)
        newmap = hp.alm2map(alms, nside,verbose=False)
        returnmap = returnmap + np.square(np.abs(newmap))

    returnmap = returnmap/numrealisations

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

    return