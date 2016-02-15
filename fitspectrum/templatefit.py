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
from scipy import special
from spectra import *
from astroutils import *
from astropy.coordinates import SkyCoord

###
# Configuration options
###
# Input data
data_filenames = ['PR2/512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits', 'PR2/512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits']
data_frequencies = [30.0, 44.0]
data_units = ['K_CMB', 'K_CMB']
data_units_use = 'mK_CMB'
# Templates. If you want an offset too, include a blank map here.
template_filenames = ["data/haslam408_dsds_Remazeilles2014.fits", "data/512_60.00smoothed_HFI_SkyMap_857_2048_R2.00_full.fits"]
template_units = ['K_CMB', 'K_CMB'] # Not set properly yet!
template_units_use = ['K_CMB', 'K_CMB'] # Not set properly yet!
# CMB
cmbsub = 0 # Set to 1 to subtract CMB map
cmbmap_filename = "PR2/512_60.00smoothed_COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits"
cmbmap_units = 'mK_CMB'
cmbspectrum_filename = "data/wmap_tt_spectrum_7yr_v4p1_nohead.txt";
# Masks
usemasks = 0
mask_filenames = ["PR2/512_60.00smoothed_COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits"]
# Map configuration
nside = 4
resolution = 60.0 # Arcmin. Assumes that the maps are already at this resolution.
regions = [[[-90, 90], [-10, 10]], [[118, 135], [20, 37]]]
# Output directory
outdir = "templatefit_output/"


### 
# Set up derived parameters and memory allocations in advance
###
npix = hp.nside2npix(nside)
lmax = 3*nside-1
num_templates = len(template_filenames)
num_maps = len(data_filenames)
num_masks = len(mask_filenames)
num_regions = len(regions)
const = get_spectrum_constants()
# Arrays
templates = np.zeros((num_templates, npix))
maps = np.zeros((num_maps, npix))
overall_mask = np.ones(npix)
region_maps = np.zeros((num_regions, npix))
coefficients = np.zeros((num_regions, num_templates))
# Output directories
ensure_dir(outdir)


###
# Prepare the masks
###
if (usemasks == 1):
	for i,filenames in enumerate(mask_filenames):

		# Read in the templates
		mask = hp.read_map(mask_filenames[i])

		# Change nside if need be
		mask_nside = hp.get_nside(mask)
		if mask_nside != nside:
			mask = hp.ud_grade(mask, nside)

		# Copy to the main mask array
		overall_mask *= mask


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
# Prepare optional CMB subtraction
###
if (cmbsub == 1):
	cmb_map = hp.read_map(cmbmap_filename)
	# Change nside if need be
	cmb_map_nside = hp.get_nside(cmb_map)
	if cmb_map_nside != nside:
		cmb_map = hp.ud_grade(cmb_map, nside)
	# Convert the units
	cmb_map *= convertunits(const, cmbmap_units, data_units_use, 1.0)


###
# Prepare the data
###
for i,filenames in enumerate(data_filenames):

	# Read in the templates
	data_map = hp.read_map(data_filenames[i])

	# Change nside if need be
	data_nside = hp.get_nside(data_map)
	if data_nside != nside:
		data_map = hp.ud_grade(data_map, nside)
	# Convert the units
	data_map *= convertunits(const, data_units[i], data_units_use, data_frequencies[i])

	# Subtract the CMB if we're doing that.
	if (cmbsub == 1):
		data_map -= cmb_map

	# Copy to the main array
	maps[i] = data_map


###
# Prepare the region masks
###
for i in range(0,num_regions):
	region_map = healpixmask(nside, regions[i][0][0], regions[i][0][1], regions[i][1][0], regions[i][1][1])
	hp.write_map(outdir + "mask_"+str(i)+".fits", region_map)
	region_maps[i] = region_map

###
# Define the CMB covariance matrix
###
# Just an identity matrix for now
covar = np.identity(npix)

# CMB covariance matrix - want this where we're not subtracting the CMB.
# NOT TESTED AS WORKING YET!
if (cmbsub == 0):
	cmbspectrum = np.loadtxt(cmbspectrum_filename)
	pixel_windowfunction = hp.sphtfunc.pixwin(nside)
	beam_windowfunction = hp.sphtfunc.gauss_beam(np.radians(resolution/60.0))
	cmb_covar = np.zeros((npix, npix))
	for i in range(0,npix):
		pos = hp.pixelfunc.pix2ang(nside, i)
		pos_i = SkyCoord("galactic", l=pos[1]*(180.0/pi), b=90.0-(pos[0]*180.0/pi),unit="deg")

		for j in range (0,npix):
			pos = hp.pixelfunc.pix2ang(nside, j)
			pos_j = SkyCoord("galactic", l=pos[1]*(180.0/pi), b=90.0-(pos[0]*180.0/pi),unit="deg")
			cmb_covar[i,j] = 0.0
			separation = np.cos(pos_i.separation(pos_j).rad)
			for l in range (0,lmax):
				# Pl = np.polynomial.legendre.legval(separation, l)
				Pl = special.lpmv(0, l, separation)
				cmb_covar[i,j] += (2.0*l+1)*cmbspectrum[l,0]*(beam_windowfunction[l]**2)*(pixel_windowfunction[l]**2)*Pl
			# print i,j,cmb_covar[i,j]

	cmb_covar /= (4.0*const['pi'])

	covar *= cmb_covar

# Invert the covariance matrix
covar_inv = np.linalg.inv(covar)


# Define a test map for now
testmap = 10.0*templates[0] + 0.1*templates[1] + cmb_map
hp.write_map(outdir + "testmap.fits", testmap)

for i in range(0,num_regions):

	###
	# Calculate the coefficients
	###
	sigma = (templates.dot(covar_inv)).dot(templates.T) # Bottom half of equation
	sigma_inv = np.linalg.inv(sigma) # Since we're dividing
	a = ((covar_inv.dot(testmap)).dot(templates.T)).dot(sigma_inv) # Top half of equation

	# Calculate the difference map between the test map and the templates with coefficients
	diff = testmap - a.dot(templates)
	hp.write_map(outdir + "diff.fits", diff)

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