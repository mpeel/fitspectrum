#!/usr/bin/env python
# Read in a map, smooth it, and write it out
#
# History:
# Mike Peel   10-Jan-2016   Initial version.
# Mike Peel   8-Mar-2016    Start to generalise to cover other numbers of maps
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def smoothmap(input, output, fwhm_arcmin, pol=False,nside_out=0,nobsmap=-1):
	# Slightly convoluted in order to also process the header...
	maps = hp.read_map(input, field=None, h=True)
	hdr = maps[-1]
	nmaps = len(maps)
	map = maps[0:nmaps-1]
	del maps
	nmaps -= 1

	if (nobsmap >=0 ):
		nobs_sum = np.sum(map[nobsmap])

	nside = hp.get_nside(map)
	smoothed_map = map
	for i in range (0,nmaps):
		smoothed_map[i][:] = hp.sphtfunc.smoothing(map[i], fwhm=np.radians(fwhm_arcmin/60.0))

	if (nside_out == 0):
		nside_out = nside

	if nside_out != nside:
		for i in range (0,nmaps):
			smoothed_map[i] = hp.ud_grade(smoothed_map[i], nside_out)
			if (i == nobsmap): # If we have an nobs map, renormalise to account for the reduction in pixel numbers.
				smoothed_map[i] *= (nside/nside_out)^2
				nobs_sum2 = np.sum(smoothed_map[nobsmap])
				print 'Smoothing nobs map: total before was ' + str(nobs_sum) + ', now is ' + str(nobs_sum2)

	# Need to modify this to write out all N maps.
	# if (pol):
		# hp.write_map(output, smoothed_map[0:3])
	# else:
		# hp.write_map(output, smoothed_map[i])
	hp.write_map(output, smoothed_map)
	
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

