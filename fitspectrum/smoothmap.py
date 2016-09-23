#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Read in a map, smooth it, and write it out
#
# History:
# v0.1 Mike Peel   10-Jan-2016   Initial version.
# v0.2 Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
# v0.3 Mike Peel   6-Sep-2016    Switch to using pyfits; use alms and window functions to smooth; add unit conversion functionality
# v0.4 Mike Peel   23-Sep-2016   Carry fits headers through to the output file. Make fwhm_arcmin optional so that the function can be used solely to ud_grade maps.
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from spectra import *
import astropy.io.fits as fits

def smoothmap(input, output, fwhm_arcmin=-1, nside_out=0,nobsmap=-1,maxnummaps=-1, frequency=100.0, units_in='',units_out='', windowfunction = []):
	ver = "0.4"

	# Read in the fits map, and put it into the format Healpix expects
	inputfits = fits.open(input)
	cols = inputfits[1].columns
	# print cols
	# exit()
	col_names = cols.names
	nmaps = len(cols)
	maps = []
	for i in range(0,nmaps):
		maps.append(inputfits[1].data.field(i))
	# Check to see whether we have nested data, and switch to ring if that is the case.
	if (inputfits[1].header['ORDERING'] == 'NESTED'):
		maps = hp.reorder(maps,n2r=True)

	if maxnummaps != -1:
		nmaps = maxnummaps
		maps = maps[0:nmaps]

	if (nobsmap >=0 ):
		nobs_sum = np.sum(maps[nobsmap])

	# Calculate the unit conversion factor
	const = get_spectrum_constants()

	pix_area = hp.nside2pixarea(nside_out)

	nside = hp.get_nside(maps)
	if (fwhm_arcmin != -1):
		conv_windowfunction = hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside)
		# conv_windowfunction_variance = hp.gauss_beam(np.radians(fwhm_arcmin/60.0)/np.sqrt(2),3*nside)
		# sigma = 2.0*np.sqrt(2.0*np.log(2.0))/np.radians(fwhm_arcmin/60.0)
		# variance_Avb = hp.nside2pixarea(nside)/(4.0*const['pi']*sigma**2)

		if (windowfunction != []):
			conv_windowfunction /= windowfunction

	smoothed_map = maps
	for i in range(0,nmaps):
		# Check that we actually want to do smoothing, as opposed to udgrading
		if fwhm_arcmin != -1:
			# Calculate the alm's, multiply them by the window function, and convert back to the map
			alms = hp.map2alm(maps[i])
			alms = hp.almxfl(alms, conv_windowfunction)
			smoothed_map[i][:] = hp.alm2map(alms, nside,verbose=False)

	if (nside_out == 0):
		nside_out = nside
	else:
		for i in range (0,nmaps):
			smoothed_map[i] = hp.ud_grade(smoothed_map[i], nside_out)
			if (units_out != ''):
				if (units_in == ''):
					unit = inputfits[1].header['TUNIT'+str(i+1)]
					unit = unit.strip()
				else:
					# Assume the user is right to have specified different input units from what is in the file.
					unit = units_in
				smoothed_map[i] *= 	convertunits(const, unit, units_out, frequency, pix_area)

			if (i == nobsmap): # If we have an nobs map, renormalise to account for the reduction in pixel numbers.
				smoothed_map[i] *= (nside/nside_out)^2
				nobs_sum2 = np.sum(smoothed_map[nobsmap])
				print 'Smoothing nobs map: total before was ' + str(nobs_sum) + ', now is ' + str(nobs_sum2)

	cols = []
	for i in range(0,nmaps):
		cols.append(fits.Column(name=col_names[i], format='E', array=smoothed_map[i]))
	cols = fits.ColDefs(cols)
	bin_hdu = fits.new_table(cols)
	bin_hdu.header = inputfits[1].header.copy(strip=False)
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['NSIDE']=nside_out
	bin_hdu.header['COMMENT']="Smoothed using Mike Peel's smoothmap.py version "+ver
	for i in range (0,nmaps):
		if (units_out != ''):
			bin_hdu.header['TUNIT'+str(i+1)] = units_out

	bin_hdu.writeto(output)


	
# Change the resolution of an nobs map. Simplest just to convert it into a variance map, and smooth that.
# def udgrade_nobs(map_in, nside_out):
	# return 1.0/np.pow(udgrade_variance(1.0/np.sqrt(map_in), nside_out), 2.0)

# # Change the resolution of an variance map, assuming that beams are Gaussian
# def udgrade_variance(map_in, nside_out):
# 	nside_in = hp.get_nside(map_in)
# 	pixarea_in = 
# 	pixarea_out = hp.nside2pixarea(nside_out)
# 	Avb = 1
# 	bv = 1
# 	return Avb * bv * map_in
# 	#σ C2 = A v b b v ∗ σ I2 ,

