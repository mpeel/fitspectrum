#!/usr/bin/env python
# Template fitting analysis for Planck data
# Aim: find magnetic dust
#
# History:
# Mike Peel   06-Nov-2015   Initial version.
# Mike Peel   04-Jan-2016   Formatting
# Mike Peel   08-Mar-2016   Major expansion, use regions throughout.
#
# Requirements:
# Numpy, healpy, matplotlib, scipy, astropy

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
# Type of template fit to do
tpl_style = 1 # 1 = I, 2 = Q&U, 3 = P
debug = 1 # 0 = silent, 1 = prints debug statements
simulation = 1 # 1 = do a simulation rather than using real data, 0 = use the real data.
# Input data
data_filenames = ['PR2/512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits', 'PR2/512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits']
data_frequencies = [30.0, 44.0]
data_units = ['K_CMB', 'K_CMB']
data_units_use = ['mK_CMB', 'mK_CMB']
# Templates. If you want an offset too, include a blank map here.
template_filenames = ["data/haslam408_dsds_Remazeilles2014.fits", "data/512_60.00smoothed_HFI_SkyMap_857_2048_R2.00_full.fits"]
data_frequencies = [0.408, 857.0] # Only needed if converting units below
template_units = ['K_CMB', 'K_CMB']
template_units_use = ['K_CMB', 'K_CMB']
# CMB
cmbsub = 0 # Set to 1 to subtract CMB map
cmbmap_filename = "PR2/512_60.00smoothed_COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits"
cmbmap_units = 'mK_CMB' # The CMB map will automatically be converted to data units as needed for the subtraction.
cmbspectrum_filename = "data/wmap_tt_spectrum_7yr_v4p1_nohead.txt";
# Masks
usemask = 0
mask_filenames = ["PR2/512_60.00smoothed_COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits"] # If you specify multiple masks here, then they will all be multiplied together.
regions = [ [[-90, 90],[-10, 10]], [[118, 135],[20, 37]] ]
# Map configuration
nside = 8
resolution = 60.0 # Arcmin. Assumes that the maps are already at this resolution.
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
coefficients = np.zeros((num_regions, num_templates))
lrange = np.arange(0,lmax)
i_pix = np.arange(0,npix)
pixel_positions = hp.pixelfunc.pix2ang(nside, i_pix)

# Output directories
ensure_dir(outdir)


###
# Prepare the masks
###
if (usemask == 1):
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
if (debug):
	print 'Reading in templates'
for i,filenames in enumerate(template_filenames):

	# Read in the templates
	template = hp.read_map(template_filenames[i])

	# Change nside if need be
	template_nside = hp.get_nside(template)
	if template_nside != nside:
		template = hp.ud_grade(template, nside)

	# Convert the units
	if (template_units[i] != template_units_use[i]):
		template *= convertunits(const, template_units[i], template_units_use[i], template_frequencies[i])

	# Copy to the main array
	templates[i] = template


###
# Prepare optional CMB subtraction
###
if (cmbsub == 1):
	if (debug):
		print 'Reading in CMB map for subtraction'
	cmb_map = hp.read_map(cmbmap_filename)
	# Change nside if need be
	cmb_map_nside = hp.get_nside(cmb_map)
	if cmb_map_nside != nside:
		cmb_map = hp.ud_grade(cmb_map, nside)

# Read in the full CMB spectrum
cmbspectrum_full = np.loadtxt(cmbspectrum_filename)


###
# Prepare the data
###
if (debug):
	print 'Reading in the data'
for i,filenames in enumerate(data_filenames):

	# Read in the templates
	data_map = hp.read_map(data_filenames[i])

	# Change nside if need be
	data_nside = hp.get_nside(data_map)
	if data_nside != nside:
		data_map = hp.ud_grade(data_map, nside)

	# Convert the units
	data_map *= convertunits(const, data_units[i], data_units_use[i], data_frequencies[i])

	# Subtract the CMB if we're doing that, including converting units for the CMB map as needed.
	if (cmbsub == 1):
		data_map -= cmb_map * convertunits(const, cmbmap_units, data_units_use[i], data_frequencies[i])

	# Copy to the main array
	maps[i] = data_map

	# If we're doing a simulation, use templates instead of the actual data.
	if (simulation):
		maps[i] = 10.0 * templates[0] + 0.1 * templates[1]


###
# Start the main loop over regions
###
for i in range(0,num_regions):
	# Prepare the region masks
	region_map = healpixmask(nside, regions[i][0][0], regions[i][0][1], regions[i][1][0], regions[i][1][1])
	hp.write_map(outdir + "mask_"+str(i)+".fits", region_map*overall_mask)

	# Need to mask the templates, datasets and pixel positions here
	mask = region_map * overall_mask
	npix_region = int(np.sum(mask))
	templates_masked = np.zeros((num_templates, npix_region))
	for j in range(0,num_templates):
		templates_masked[j] = templates[j][mask == 1]
	data_masked = np.zeros((num_maps, npix_region))
	for j in range(0,num_maps):
		data_masked[j] = maps[j][mask == 1]
	positions_masked = np.zeros((2,npix_region))
	positions_masked[0] = pixel_positions[0][mask == 1]
	positions_masked[1] = pixel_positions[1][mask == 1]
	print np.shape(data_masked)
	###
	# Start loop over data maps
	###
	for j in range(0,num_maps):

		###
		# Define the CMB covariance matrix
		###
		# Just an identity matrix for now
		covar = np.identity(npix_region)

		# CMB covariance matrix - want this where we're not subtracting the CMB.
		if (cmbsub == 0):
			cmbspectrum = cmbspectrum_full[:lmax, 1]
			pixel_windowfunction = hp.sphtfunc.pixwin(nside)
			pixel_windowfunction = pixel_windowfunction[:lmax]
			pixel_windowfunctionsq = pixel_windowfunction**2
			beam_windowfunction = hp.sphtfunc.gauss_beam(np.radians(resolution/60.0))
			beam_windowfunction = beam_windowfunction[:lmax]
			beam_windowfunctionsq = beam_windowfunction**2
			cmb_covar = np.zeros((npix_region, npix_region))
			# Calculate all of the separations between pixels in one big array
			pos_coord = SkyCoord(frame="galactic", l=positions_masked[1][:]*(180.0/pi), b=90.0-(positions_masked[0][:]*180.0/pi),unit="deg")
			separations = np.array([np.cos(pos_coord[i].separation(pos_coord).rad) for i in range(npix_region)])
			separations = np.reshape(separations, (npix_region, npix_region))
			print np.shape(separations)

			for i_cov in range(0,npix_region):
				for j_cov in range (0,npix_region):
					cmb_covar[i_cov,j_cov] = 0.0
					Pl = special.lpmv(0, lrange, separations[i_cov,j_cov])
					cmb_covar[i_cov,j_cov] = np.sum((2.0*lrange+1)*cmbspectrum*(beam_windowfunctionsq)*(pixel_windowfunctionsq)*Pl)
					# print i,j,cmb_covar[i,j]

			cmb_covar /= (4.0*const['pi'])

			covar += cmb_covar
			del separations # We don't need this any more.
			del cmb_covar

		# Invert the covariance matrix
		covar_inv = np.linalg.inv(covar)

		###
		# Calculate the coefficients
		###
		sigma = (templates_masked.dot(covar_inv)).dot(templates_masked.T) # Bottom half of equation
		sigma_inv = np.linalg.inv(sigma) # Since we're dividing
		a = ((covar_inv.dot(data_masked[j])).dot(templates_masked.T)).dot(sigma_inv) # Top half of equation

		# Calculate the difference map between the test map and the templates with coefficients
		diff = data_masked[j] - a.dot(templates_masked)
		# hp.write_map(outdir + "diff.fits", diff)

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