#!/usr/bin/env python
# Functions for templatefit.py

import numpy as np
from scipy import special
from astropy.coordinates import SkyCoord

def templatefit(covar, templates, data):
	covar_inv = np.linalg.inv(covar)
	return templatefit_givencovarinv(covar, templates, data, covar_inv)

def templatefit_givencovarinv(covar, templates, data, covar_inv):
	a_arr = np.dot(templates, np.dot(covar_inv, templates.T)) # Bottom half of equation
	a_inv = np.linalg.inv(a_arr) # Since we're dividing
	b_arr = np.dot(data, np.dot(covar_inv, templates.T)) # Top half of the equation
	result = np.dot(b_arr, a_inv) # Combine the top and bottom halves of the equation
	a_err = np.sqrt(np.diag(a_inv))
	diff = data - result.dot(templates)
	chisq = (diff.dot(covar_inv)).dot(diff.T)

	return result, a_err, chisq


def calc_cmbcovar(const, npix, pixel_positions, lrange, cmbspectrum, beam_windowfunctionsq, pixel_windowfunctionsq):

	cmb_covar = np.zeros((npix, npix))
	# Calculate all of the separations between pixels in one big array
	pos_coord = SkyCoord(frame="galactic", l=pixel_positions[1][:]*(180.0/const['pi']), b=90.0-(pixel_positions[0][:]*180.0/const['pi']),unit="deg")
	separations = np.array([np.cos(pos_coord[i_temp].separation(pos_coord).rad) for i_temp in range(npix)])
	separations = np.reshape(separations, (npix, npix))

	for i_cov in range(0,npix):
		for j_cov in range (0,npix):
			Pl = special.lpmv(0, lrange, separations[i_cov,j_cov])
			cmb_covar[i_cov,j_cov] = np.sum((2.0*lrange+1.0)*cmbspectrum*(beam_windowfunctionsq)*(pixel_windowfunctionsq)*Pl)

	cmb_covar /= (4.0*const['pi'])
	del separations
	del Pl

	return cmb_covar

def chisq_minimize_indices(indices, num_maps, npix_region, num_templates, templates_masked, data_frequencies, covar, data_masked, covar_inv):
	# print indices
	# Create the scaled templates for these indices
	templates_masked_scaled = np.zeros((num_templates, num_maps * npix_region))
	for j_cov in range(0,num_maps):
		for i_cov in range(0,npix_region):
			for k_cov in range(0,num_templates):
				templates_masked_scaled[k_cov][j_cov*npix_region+i_cov] = templates_masked[k_cov][i_cov] * (data_frequencies[j_cov]/data_frequencies[0])**indices[k_cov]

	# Do the template fitting
	a, a_err, chisq = templatefit_givencovarinv(covar, templates_masked_scaled, data_masked, covar_inv)

	# Return the chisq
	return chisq