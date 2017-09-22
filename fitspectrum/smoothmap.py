#!/usr/bin/env python
# -*- coding: utf-8  -*-
# Read in a map, smooth it, and write it out
#
# History:
# v0.1  Mike Peel   10-Jan-2016   Initial version.
# v0.2  Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
# v0.3  Mike Peel   6-Sep-2016    Switch to using pyfits; use alms and window functions to smooth; add unit conversion functionality
# v0.4  Mike Peel   23-Sep-2016   Carry fits headers through to the output file. Make fwhm_arcmin optional so that the function can be used solely to ud_grade maps.
# v0.5  Mike Peel   23-Sep-2016   Adding calc_variance_windowfunction (port of Paddy Leahy's IDL code) to properly smooth variance maps.
# v0.6  Mike Peel   26-Sep-2016   Support Nobs maps, and conversion to variance maps, plus debugging/tidying.
# v0.7  Mike Peel   27-Jan-2017   Bug fix - ud_grading the variance maps should include a factor of (nside_orig/nside_new)^2
# V0.8  Adam Barr   Various       Various tweaks
# v0.9  Mike Peel   17 Sep 2017   Tweaks, add rescale parameter
# v0.9a Mike Peel   21 Sep 2017   More bug fixes
# v0.9b Mike Peel   22 Sep 2017   Add nosmooth and outputmaps parameters
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import scipy as sp
import math as m
import matplotlib.pyplot as plt
from spectra import *
import astropy.io.fits as fits
from scipy import special
import os.path

def smoothmap(input, output, fwhm_arcmin=-1, nside_out=0,maxnummaps=-1, frequency=100.0, units_in='',units_out='', windowfunction = [],nobs_out=False,variance_out=True, sigma_0 = -1, sigma_0_unit='', rescale=1.0, nosmooth=[], outputmaps=[]):
	ver = "0.9a"

	if (os.path.isfile(output)):
		print "You already have a file with the output name " + output + "! Not going to overwrite it. Move it, or set a new output filename, and try again!"
		exit()

	# Check to see if we have a sigma_0 value to use when converting from Nobs maps and back.
	no_sigma_0 = False
	if (sigma_0 == -1):
		no_sigma_0 = True

	# Read in the fits map, and put it into the format Healpix expects
	inputfits = fits.open(input)
	cols = inputfits[1].columns
	col_names = cols.names
	nmaps = len(cols)
	nmaps_orig = nmaps
	maps = []
	for i in range(0,nmaps):
		maps.append(inputfits[1].data.field(i))
	# Check to see whether we have nested data, and switch to ring if that is the case.
	if (inputfits[1].header['ORDERING'] == 'NESTED'):
		maps = hp.reorder(maps,n2r=True)

	# Crop to just have the maps we want to output
	if maxnummaps != -1:
		nmaps = maxnummaps
	if outputmaps == []:
		outputmaps = range(0,nmaps)
	maps=maps[outputmaps]
	noutputmaps = len(outputmaps)
	newheader = inputfits[1].header.copy(strip=False)
	for i in range(0,nmaps_orig):
		if i < noutputmaps-1:
			newheader['TTYPE'+str(i+1)] = newheader['TTYPE'+str(outputmaps[i]+1)]
			newheader['TFORM'+str(i+1)] = newheader['TFORM'+str(outputmaps[i]+1)]
			newheader['TUNIT'+str(i+1)] = newheader['TUNIT'+str(outputmaps[i]+1)]
		else:
			del newheader['TTYPE'+str(i+1)]
			del newheader['TFORM'+str(i+1)]
			del newheader['TUNIT'+str(i+1)]
	newheader['TFIELDS'] = noutputmaps
	nmaps = noutputmaps
	print newheader

	# Calculate the unit conversion factor
	const = get_spectrum_constants()

	pix_area = hp.nside2pixarea(nside_out)

	nside = hp.get_nside(maps)
	if (fwhm_arcmin != -1):
		conv_windowfunction = hp.gauss_beam(np.radians(fwhm_arcmin/60.0),3*nside)
		if (windowfunction != []):
			window_len = len(conv_windowfunction)
			beam_len = len(windowfunction)

			if (beam_len > window_len):
				windowfunction  = windowfunction[0:len(conv_windowfunction)]
			else:
				conv_windowfunction = conv_windowfunction[0:len(windowfunction)]

			conv_windowfunction /= windowfunction
		# Normalise window function
		conv_windowfunction /= conv_windowfunction[0]

		# Check whether we'll need to smooth variances too.
		test = False
		for i in range(0,nmaps):
			if ('cov' in newheader['TTYPE'+str(i+1)]) or ('N_OBS' in newheader['TTYPE'+str(i+1)]):
				test = True
		if test:
			print 'Covariance maps detected. Calculating variance window function (this may take a short while)'
			conv_windowfunction_variance = calc_variance_windowfunction(conv_windowfunction)
			conv_windowfunction_variance /= conv_windowfunction_variance[0]
			print conv_windowfunction_variance[0]
			print 'Done! Onwards...'

			print conv_windowfunction
			print conv_windowfunction_variance
			plt.xscale('log')
			plt.yscale('log')
			plt.plot(conv_windowfunction,label='Window function')
			plt.plot(conv_windowfunction_variance,label='Variance window function')
			# plt.legend()
			plt.savefig('test_plotwindowfunction.png')
			# # exit()

	# Do the smoothing
	print "Smoothing the maps"
	smoothed_map = maps
	for i in range(0,nmaps):
		# Check that we actually want to do smoothing, as opposed to udgrading. Also check to see if this is in the list of maps to not smooth
		if fwhm_arcmin != -1 and (i not in nosmooth):
			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				print 'Column '+str(i)+' is an N_OBS map ('+newheader['TUNIT'+str(i+1)]+') - converting to variance map.'
				print np.sum(maps[i])
				maps[i] = conv_nobs_variance_map(maps[i], sigma_0)
				print np.sum(maps[i])
				if (nobs_out == False and no_sigma_0 == False):
					# We don't want to convert back later.
					print 'test'
					newheader['TTYPE'+str(i+1)] = 'cov'
					newheader['TUNIT'+str(i+1)] = '('+sigma_0_unit+')^2'

			# Calculate the alm's, multiply them by the window function, and convert back to the map
			alms = hp.map2alm(maps[i])
			if ('cov' in newheader['TTYPE'+str(i+1)]) or ('N_OBS' in newheader['TTYPE'+str(i+1)]):
				print 'Column '+str(i)+' is a covariance matrix ('+newheader['TUNIT'+str(i+1)]+') - smoothing appropriately.'
				alms = hp.almxfl(alms, conv_windowfunction_variance)
			else:
				alms = hp.almxfl(alms, conv_windowfunction)
			smoothed_map[i][:] = hp.alm2map(alms, nside,verbose=False)
			print np.sum(smoothed_map[i])

			if ('N_OBS' in newheader['TTYPE'+str(i+1)]) and (nobs_out or no_sigma_0):
				print 'You\'ve either asked for an N_OBS map to be returned, or not set sigma_0, so you will get an N_OBS map returned in your data!'
				print np.sum(smoothed_map[i])
				smoothed_map[i] = conv_nobs_variance_map(smoothed_map[i], sigma_0)
				print np.sum(smoothed_map[i])
				newheader['TTYPE'+str(i+1)] = 'N_OBS'

	# Do the ud_grading
	print "ud_grading the maps (if needed)"
	nobs_sum = 0
	if (nside_out == 0):
		nside_out = nside
	else:
		for i in range (0,nmaps):
			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				nobs_sum = np.sum(smoothed_map[i])

			# Check to see which type of map we have, and adjust the factor of (nside/nside_out)^power appropriately
			power = 0
			if ('cov' in newheader['TTYPE'+str(i+1)]):
				power = 2
			elif ('N_OBS' in newheader['TTYPE'+str(i+1)]) or ('Hits' in newheader['TTYPE'+str(i+1)]):
				power = -2
			print power
			smoothed_map[i] = hp.ud_grade(smoothed_map[i], nside_out, power=power)

			if 'N_OBS' in newheader['TTYPE'+str(i+1)]:
				null = 0
			else:
				# If we don't have an N_OBS map, then we might want to convert the units.
				if (units_out != ''):
					if (units_in == ''):
						unit = newheader['TUNIT'+str(i+1)]
						unit = unit.strip()
					else:
						# Assume the user is right to have specified different input units from what is in the file.
						unit = units_in

					bin_hdu.header['TUNIT'+str(i+1)] = units_out
					power = 1.0
					if ('^2' in unit):
						power = 2.0
						unit = unit.replace(")^2",'').replace('(','')
						bin_hdu.header['TUNIT'+str(i+1)] = '('+units_out+")^2"
					print unit + " " + str(power)
					conversion = convertunits(const, unit, units_out, frequency, pix_area)
					print conversion
					smoothed_map[i] = smoothed_map[i] * conversion**power

	# All done - now just need to write it to disk.
	print "Writing maps to disk: " + output
	cols = []
	for i in range(0,nmaps):
		if ('cov' in newheader['TTYPE'+str(i+1)]):
			smoothed_map[i] = smoothed_map[i] * rescale**2
		else:
			smoothed_map[i] = smoothed_map[i] * rescale

		cols.append(fits.Column(name=col_names[i], format='E', array=smoothed_map[i]))
		
	cols = fits.ColDefs(cols)
	bin_hdu = fits.new_table(cols)
	# print newheader
	bin_hdu.header = newheader
	# print bin_hdu.header
	bin_hdu.header['ORDERING']='RING'
	bin_hdu.header['POLCONV']='COSMO'
	bin_hdu.header['PIXTYPE']='HEALPIX'
	bin_hdu.header['NSIDE']=nside_out
	bin_hdu.header['COMMENT']="Smoothed using Mike Peel's smoothmap.py version "+ver +" modified by Adam Barr"
	# print bin_hdu.header
	
	bin_hdu.writeto(output)

	return

def conv_nobs_variance_map(inputmap, sigma_0):
	newmap = sigma_0**2 / inputmap
	return newmap

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
	ll = np.arange(0,nbl,1)
	
	# Choose scale to match size of beam. First find half-power point
	lhalf = [ n for n,i in enumerate(conv_windowfunction) if i<0.5 ][0]
	
	# Calculate beam out to about 40 * half-power radius, roughly (note that
	# this gives 100 points out to the half-power point).
	numelements = 4000
	rad = np.arange(0,numelements,1)*10.0*const['pi']/(float(lhalf)*float(numelements))
	
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
	print 'Peak of convolving beam is ' + str(conva[0]) + " (check: " + str(np.max(conva)) + ")"

	# Square convolving beam and convert back to window function
	mult = sinrad*conva**2
	cvbl = np.zeros(nbl)
	print nbl
	for l in range(0,nbl):
		cvbl[l] = int_tabulated(rad,mult*lgndr[:,l])

	# Put in 2pi normalization factor:
	cvbl = 2.0*const['pi']*cvbl

	return cvbl

# EOF
