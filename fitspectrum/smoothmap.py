#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Read in a map, smooth it, and write it out
#
# History:
# Mike Peel   10-Jan-2016   Initial version.
# Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
# Mike Peel   6-Sep-2016    Switch to using pyfits; use alms and window functions to smooth; add unit conversion functionality
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from spectra import *
import astropy.io.fits as fits

def smoothmap(input, output, fwhm_arcmin, nside_out=0,nobsmap=-1,maxnummaps=-1, frequency=0, units_in='JyPix',units_out='JyPix', windowfunction = []):

	# Read in the fits map, and put it into the format Healpix expects
	inputfits = fits.open(input)
	cols = inputfits[1].columns
	col_names = cols.names
	nmaps = len(cols)
	maps = []
	for i in range(0,nmaps):
		maps.append(inputfits[1].data.field(i))
	# Check to see whether we have nested data, and switch to ring if that is the case.
	if (inputfits[1].header['ORDERING'] == 'NESTED'):
		print 'Reordering!'
		maps = hp.reorder(maps,n2r=True)

	if maxnummaps != -1:
		nmaps = maxnummaps
		maps = maps[0:nmaps]

	if (nobsmap >=0 ):
		nobs_sum = np.sum(maps[nobsmap])

	# Calculate the unit conversion factor
	const = get_spectrum_constants()
	pix_area = hp.nside2pixarea(nside_out)
	unit_factor = convertunits(const, units_in, units_out, frequency, pix_area)

	nside = hp.get_nside(maps)
	conv_windowfunction = hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside)
	if (windowfunction != []):
		conv_windowfunction /= windowfunction

	smoothed_map = maps
	for i in range(0,nmaps):
		# Calculate the alm's, multiply them by the window function, and convert back to the map
		alms = hp.map2alm(maps[i])
		alms = hp.almxfl(alms, conv_windowfunction)
		smoothed_map[i][:] = hp.alm2map(alms, nside)

	if (nside_out == 0):
		nside_out = nside
	else:
		for i in range (0,nmaps):
			smoothed_map[i] = hp.ud_grade(smoothed_map[i], nside_out)
			smoothed_map[i] *= unit_factor

			if (i == nobsmap): # If we have an nobs map, renormalise to account for the reduction in pixel numbers.
				smoothed_map[i] *= (nside/nside_out)^2
				nobs_sum2 = np.sum(smoothed_map[nobsmap])
				print 'Smoothing nobs map: total before was ' + str(nobs_sum) + ', now is ' + str(nobs_sum2)

	cols = []
	for i in range(0,nmaps):
		cols.append(fits.Column(name=col_names[i], format='E', array=smoothed_map[i]))
	cols = fits.ColDefs(cols)
	bin_hdu = fits.new_table(cols)
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['NSIDE']=nside
	bin_hdu.writeto(output)


	
# Change the resolution of an nobs map. Simplest just to convert it into a variance map, and smooth that.
# def udgrade_nobs(map_in, nside_out):
# 	return 1.0/udgrade_variance(1.0/map_in, nside_out)

# # Change the resolution of an variance map, assuming that beams are Gaussian
# def udgrade_variance(map_in, nside_out):
# 	nside_in = hp.get_nside(map)
# 	pixarea_in = hp.nside2pixarea(nside_in)
# 	pixarea_out = hp.nside2pixarea(nside_out)
# 	Avb = 1
# 	bv = 1
# 	return Avb * bv * map_in
# 	#σ C2 = A v b b v ∗ σ I2 ,

