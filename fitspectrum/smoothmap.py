#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Read in a map, smooth it, and write it out
#
# History:
# v0.1 Mike Peel   10-Jan-2016   Initial version.
# v0.2 Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
# v0.3 Mike Peel   6-Sep-2016    Switch to using pyfits; use alms and window functions to smooth; add unit conversion functionality
# v0.4 Mike Peel   23-Sep-2016   Carry fits headers through to the output file. Make fwhm_arcmin optional so that the function can be used solely to ud_grade maps.
# v0.5 Mike Peel   23-Sep-2016   Adding calc_variance_windowfunction (port of Paddy Leahy's IDL code) to properly smooth variance maps.
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import scipy as sp
import matplotlib.pyplot as plt
from spectra import *
import astropy.io.fits as fits
from scipy import special

def smoothmap(input, output, fwhm_arcmin=-1, nside_out=0,nobsmap=-1,maxnummaps=-1, frequency=100.0, units_in='',units_out='', windowfunction = []):
	ver = "0.5"

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
		if (windowfunction != []):
			conv_windowfunction /= windowfunction

		print conv_windowfunction
		# Check whether we'll need to smooth variances too.
		test = False
		for i in range(0,nmaps):
			if 'cov' in inputfits[1].header['TTYPE'+str(i+1)]:
				test = True
		if test:
			print 'Covariance maps detected. Calculating variance window function (this may take a short while)'
			conv_windowfunction_variance = calc_variance_windowfunction(conv_windowfunction)
			print 'Done! Onwards...'

	smoothed_map = maps
	for i in range(0,nmaps):
		# Check that we actually want to do smoothing, as opposed to udgrading
		if fwhm_arcmin != -1:
			# Calculate the alm's, multiply them by the window function, and convert back to the map
			alms = hp.map2alm(maps[i])
			if 'cov' in inputfits[1].header['TTYPE'+str(i+1)]:
				print 'Column '+str(i)+'is a covariance matrix ('+inputfits[1].header['TUNIT'+str(i+1)]+') - smoothing appropriately.'
				alms = hp.almxfl(alms, conv_windowfunction_variance)
			else:
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

# Code to replicate IDL's INT_TABULATED function
# From http://stackoverflow.com/questions/14345001/idls-int-tabulate-scipy-equivalent
def int_tabulated(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = sp.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret

def calc_variance_windowfunction(conv_windowfunction):
	# Calculate the window function for variance maps.
	# Based on IDL code by Paddy Leahy, 'write_cvbl_planck3.pro'

	const = get_spectrum_constants()
	nbl = len(conv_windowfunction)
	ll = np.linspace(0,1,nbl)
	# Choose scale to match size of beam. First find half-power point
	lhalf = [ n for n,i in enumerate(conv_windowfunction) if i<0.5 ][0]

	# Calculate beam out to about 40 * half-power radius, roughly (note that
	# this gives 100 points out to the half-power point).
	numelements = 4000
	rad = np.linspace(0,1,numelements)*10.0*const['pi']/(float(lhalf)*float(numelements))

	x = np.cos(rad)
	sinrad = np.sin(rad)

	lgndr = np.zeros((numelements,nbl))
	for i in range(0,nbl):
		lgndr[:,i] = special.lpmv(0, i, x)
	
	# Generate radial profile of convolving beam:
	conva = np.zeros(numelements)
	for j in range(0,numelements):
		conva[j] = np.sum((ll+0.5)*conv_windowfunction*lgndr[j,:])

	conva = conva / (2.0*const['pi'])

	# print 'Peak of convolving beam is ' + str(conva[0]) + + " (check: " + str(np.max(conva)) + ")"

	# Square convolving beam and convert back to window function
	mult = sinrad*conva**2
	cvbl = np.zeros(nbl)
	for l in range(0,nbl):
		cvbl[l] = int_tabulated(rad,mult*lgndr[:,l])

	# Put in 2pi normalization factor:
	cvbl = 2.0*const['pi']*cvbl

	return cvbl

# EOF
