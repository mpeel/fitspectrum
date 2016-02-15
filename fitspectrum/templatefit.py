#!/usr/bin/env python
# Template fitting analysis for Planck data
# Aim: find magnetic dust
#
# History:
# Mike Peel   06-Nov-2015   Initial version.
# Mike Peel   04-Jan-2016   Formatting
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from spectra import *


###
# Configuration options
###
template_filenames = ["data/haslam408_dsds_Remazeilles2014.fits", "data/512_60.00smoothed_HFI_SkyMap_857_2048_R2.00_full.fits"]
cmbspectrum_filename = "data/wmap_tt_spectrum_7yr_v4p1_nohead.txt";
nside = 4
regions = [[118, 135], [20, 37]]
outputdir = "output"


### 
# Set up derived parameters and memory allocations in advance
###
npix = hp.nside2npix(nside)
num_templates = len(template_filenames)
num_regions = len(regions)
const = get_spectrum_constants()
# Arrays
templates = np.zeros((num_templates, npix))
coefficients = np.zeros(num_templates)


###
# Prepare the templates
###
for i,filenames in enumerate(template_filenames):

	# Read in the templates
	template = hp.read_map(template_filenames[i])

	# Change nside if need be
	template_nside = hp.get_nside(template)
	if template_nside != nside:
		template = hp.ud_grade(template, nside)

	# Copy to the main array
	templates[i] = template

###
# Prepare the masks
###



###
# Read in the data
###

# Define a test map for now
testmap = 10.0*templates[0] + 0.1*templates[1]



###
# Define the CMB covariance matrix
###
# Just an identity matrix for now
covar = np.identity(npix)

# CMB covariance matrix
cmbspectrum = np.loadtxt(cmbspectrum_filename)

# Invert the covariance matrix
covar_inv = np.linalg.inv(covar)


###
# Calculate the coefficients
###
sigma = (templates.dot(covar_inv)).dot(templates.T)
sigma_inv = np.linalg.inv(sigma)
a = ((covar_inv.dot(testmap)).dot(templates.T)).dot(sigma_inv)
diff = testmap - a.dot(templates)
chisq = (diff.dot(covar_inv)).dot(diff.T)

abar = sum(a.dot(sigma_inv)) / sum(sigma_inv)
a_err = np.sqrt(np.diag(sigma_inv))

###
# Output the results
###
print a
print chisq
print abar
print a_err

# That's all, folks!
# EOF