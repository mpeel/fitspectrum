#!/usr/bin/env python
# Template fitting analysis for Planck data
# Aim: find magnetic dust
#
# History:
# Mike Peel   06-Nov-2015   Initial version.
# Mike Peel   04-Jan-2016   Formatting
# Mike Peel   08-Mar-2016   Major expansion, use regions throughout.
# Mike Peel   29-Mar-2016   Add num_runs and simulation values. Debugging.
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
num_runs = 1000 # Set to 1 for normal use, or more than that if doing a simulation.
simulation_add_cmb = 0 # Use a simulated CMB map too?
cmbsub = 0 # Set to 1 to subtract CMB map, or 0 to use the covariance matrix
cmb_use_covar = 0 # Set to 1 to use the covariance matrix, or 0 to not use it.
usemask = 1
outdir = "templatefit_davies_wmap1_sim_output/" # Output directory. Will be created if it doesn't already exist.
# Input data
# data_filenames = ['wmap9/512_60.00smoothed_wmap_band_iqumap_r9_9yr_K_v5.fits', 'wmap9/512_60.00smoothed_wmap_band_iqumap_r9_9yr_Ka_v5.fits', 'wmap9/512_60.00smoothed_wmap_band_iqumap_r9_9yr_Q_v5.fits', 'wmap9/512_60.00smoothed_wmap_band_iqumap_r9_9yr_V_v5.fits', 'wmap9/512_60.00smoothed_wmap_band_iqumap_r9_9yr_W_v5.fits'] # WMAP 9-year data
data_filenames = ['wmap1/512_60.00smoothed_map_k_imap_yr1_v1.fits']#, 'wmap1/512_60.00smoothed_map_ka_imap_yr1_v1.fits', 'wmap1/512_60.00smoothed_map_q_imap_yr1_v1.fits', 'wmap1/512_60.00smoothed_map_v_imap_yr1_v1.fits', 'wmap1/512_60.00smoothed_map_w_imap_yr1_v1.fits'] # WMAP 1-year data
#data_filenames = ['PR2/512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits', 'PR2/512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits']
# data_variance_column = [3, 3, 3, 3, 3] # Starting from 0
data_variance_column = [1, 1, 1, 1, 1] # Starting from 0
data_sigma_0 = [1.418, 1.429, 2.105, 2.854, 5.253] # Set to -1 if the variance maps are actually variance maps, or the appropriate value if they are Nobs maps. If useful values, needs to be in the same units as the input maps.
data_frequencies = [22.8, 33.0, 40.7, 60.8, 93.5]
data_units = ['mK_CMB', 'mK_CMB', 'mK_CMB', 'mK_CMB', 'mK_CMB']
data_units_use = ['uK_RJ', 'uK_RJ', 'uK_RJ', 'uK_RJ', 'uK_RJ']
# Templates. If you want an offset too, include a blank map here.
template_names = ['sync', 'halpha', 'dust', 'offset']
template_filenames = ["2p3ghz/64_60.00smoothed_lambda_haslam408_dsds_halpha_ddd3.fits", "2p3ghz/64_60.00smoothed_ha_correct_33_h256.fits", "2p3ghz/64_60.00smoothed_lambda_fds_dust_94GHz.fits", "2p3ghz/64_const.fits"]
#template_filenames = ["2p3ghz/64_60.00smoothed_map_2300mhz_1deg_halpha_ddd3.fits", "2p3ghz/64_60.00smoothed_ha_correct_33_h256.fits", "2p3ghz/64_60.00smoothed_lambda_fds_dust_94GHz.fits"]
#template_filenames = ["data/haslam408_dsds_Remazeilles2014.fits", "data/512_60.00smoothed_HFI_SkyMap_857_2048_R2.00_full.fits"]
template_frequencies = [0.408, 1.0, 94.0, 1.0] # Only needed if converting units below
template_units = ['K_CMB', 'R', 'K_CMB', 'K_CMB']
template_units_use = ['K_CMB', 'R', 'K_CMB', 'K_CMB']
simulation_values = [10.0, 1.0, 0.0, 0.0]
simulation_noise = [1.0] # In whichever units are used for the templates
# CMB
cmbmap_filename = "PR2/512_60.00smoothed_COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits"
cmbmap_units = 'mK_CMB' # The CMB map will automatically be converted to data units as needed for the subtraction.
# cmbspectrum_filename = "data/wmap_tt_spectrum_7yr_v4p1_nohead.txt"; # WMAP 7-year
#cmbspectrum_units = 'uK_CMB'
#cmbspectrum_Cl = False # If the CMB spectrum is already in Cl, set this to true, otherwise false and we'll assume l(l+1)C_l/2pi
cmbspectrum_filename = "wmap1/wmap_1yr_kband_powerspectrum.fits.dat"; # WMAP 1-year
cmbspectrum_units = 'uK_CMB'
cmbspectrum_Cl = True # If the CMB spectrum is already in Cl, set this to true, otherwise false and we'll assume l(l+1)C_l/2pi
cmbspectrum_minl = 2
# Masks
mask_filenames = ["wmap1/64_60.00smoothed_map_kp2_mask_yr1_v1_2.fits"] # If you specify multiple masks here, then they will all be multiplied together.
#mask_filenames = ["2p3ghz/64_60.00smoothed_wmap_ext_temperature_analysis_mask_r9_7yr_v4_2.fits"] # If you specify multiple masks here, then they will all be multiplied together.
regions = [ [[245, 260],[21, 31]]]#, [[140, 155],[15, 20]], [[200, 230],[-48, -41]], [[250, 260],[-35, -25]], [[90, 97],[-30, -13]], [[118, 135],[20, 37]], [[300, 315],[35, 45]], [[227, 237],[12, 18]], [[145, 165],[-38, -30]], [[300, 320],[-40, -30]], [[33, 45],[50, 70]], [[270, 310],[55, 70]], [[350, 365],[-50, -35]], [[70, 90],[20, 30]], [[76, 84],[-50, -30]]] # Davies et al. (2006) 15 regions
# Map configuration
nside = 64 # Maps will be ud_graded to this as needed.
resolution = 60.0 # Arcmin. Assumes that the maps are already at this resolution.


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
variance = np.zeros((num_maps, npix))
overall_mask = np.ones(npix)
coefficients = np.zeros((num_regions, num_templates))
lrange = np.arange(0,lmax)
i_pix = np.arange(0,npix)
pixel_positions = hp.pixelfunc.pix2ang(nside, i_pix)
a = np.zeros((num_runs, num_regions, num_maps, num_templates))
abar = np.zeros((num_runs, num_regions, num_maps, num_templates))
a_err = np.zeros((num_runs, num_regions, num_maps, num_templates))
chisq = np.zeros((num_runs, num_regions, num_maps))
num_region_pixels = np.zeros((num_regions))

# Output directories
ensure_dir(outdir)


###
# Prepare the masks
###
if (debug):
	print ''
	print '** Preparing masks'
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
	print ''
	print '** Reading in templates'
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
		print ''
		print '** Reading in CMB map for subtraction'
	cmb_map = hp.read_map(cmbmap_filename)
	# Change nside if need be
	cmb_map_nside = hp.get_nside(cmb_map)
	if cmb_map_nside != nside:
		cmb_map = hp.ud_grade(cmb_map, nside)


# Read in the full CMB spectrum. Also calculate pixel and beam window functions here.
if (cmb_use_covar):
	print (convertunits(const, cmbspectrum_units, data_units_use[0], data_frequencies[0]))**2
	cmbspectrum_full = np.loadtxt(cmbspectrum_filename)
	cmbspectrum = cmbspectrum_full[:lmax, 1] * (convertunits(const, cmbspectrum_units, data_units_use[0], data_frequencies[0]))**2
	cmbspectrum = np.roll(cmbspectrum, cmbspectrum_minl)
	if (cmbspectrum_Cl == False):
		cmbspectrum *= (2.0 * const['pi'] / (lrange*(lrange+1.0)))
	for i in range(0,cmbspectrum_minl):
		cmbspectrum[i] = 0.0
	pixel_windowfunction = hp.sphtfunc.pixwin(nside)
	pixel_windowfunction = pixel_windowfunction[:lmax]
	pixel_windowfunctionsq = pixel_windowfunction**2
	beam_windowfunction = hp.sphtfunc.gauss_beam(np.radians(resolution/60.0))
	beam_windowfunction = beam_windowfunction[:lmax]
	beam_windowfunctionsq = beam_windowfunction**2

###
# Loop over runs
###
for r in range(0,num_runs):
	###
	# Prepare the data
	###
	if (debug):
		print ''
		print '* Run ' + str(r) + ' of ' + str(num_runs)
		print '** Reading in the data'

	for i,filenames in enumerate(data_filenames):
		print i
		# Read in the templates
		data_map = hp.read_map(data_filenames[i],field=None)

		# Change nside if need be
		data_nside = hp.get_nside(data_map)
		if data_nside != nside:
			data_map = hp.ud_grade(data_map, nside)
			rescale = (data_nside/nside)**2

		# Copy the variance map over
		variance[i] = data_map[data_variance_column[i]]
		if (data_sigma_0[i] != -1):
			print "Average variance:" + str(np.mean(variance[i]) * rescale)
			print "Average value:" + str((data_sigma_0[i] * convertunits(const, data_units[i], data_units_use[i], data_frequencies[i])))
			variance[i] = (data_sigma_0[i] * convertunits(const, data_units[i], data_units_use[i], data_frequencies[i]))**2 / (variance[i] * rescale)
			print "Average variance: " + str(np.mean(variance[i]))

		# Convert the units of the map data
		data_map[0] *= convertunits(const, data_units[i], data_units_use[i], data_frequencies[i])

		# Subtract the CMB if we're doing that, including converting units for the CMB map as needed.
		if (cmbsub == 1):
			data_map[0] -= cmb_map * convertunits(const, cmbmap_units, data_units_use[i], data_frequencies[i])

		# Copy to the main array
		maps[i] = data_map[0][:]

		# If we're doing a simulation, use templates instead of the actual data.
		if (simulation):
			maps[i] = simulation_values[0] * templates[0]
			for j in range(1,num_templates):
				maps[i] += simulation_values[j] * templates[j]
			maps[i] += np.random.randn(npix) * simulation_noise[i]
			if (simulation_add_cmb == 1):
				sim_cmb_map = hp.sphtfunc.synfast(cmbspectrum, nside,lmax, fwhm=np.radians(resolution/60.0))
				hp.write_map(outdir + "cmb_simulation_"+str(r)+".fits", sim_cmb_map)
				maps[i] += sim_cmb_map
				# maps[i] += cmb_map * convertunits(const, cmbmap_units, data_units_use[i], data_frequencies[i])
			variance[i] = simulation_noise[i]**2

	###
	# Start the main loop over regions
	###
	if (debug):
		print ''
		print '** Looping over regions'
	for i in range(0,num_regions):
		if (debug):
			print ''
			print '*** Region ' + str(i)

		# Prepare the region masks
		region_map = healpixmask(nside, regions[i][0][0], regions[i][0][1], regions[i][1][0], regions[i][1][1])
		hp.write_map(outdir + "mask_"+str(i)+".fits", region_map*overall_mask)

		# Need to mask the templates, datasets and pixel positions here
		mask = region_map * overall_mask
		npix_region = int(np.sum(mask))
		num_region_pixels[i] = npix_region
		templates_masked = np.zeros((num_templates, npix_region))
		for j in range(0,num_templates):
			templates_masked[j] = templates[j][mask == 1]
		data_masked = np.zeros((num_maps, npix_region))
		for j in range(0,num_maps):
			data_masked[j] = maps[j][mask == 1]
		variance_masked = np.zeros((num_maps, npix_region))
		for j in range(0,num_maps):
			variance_masked[j] = variance[j][mask == 1]

		positions_masked = np.zeros((2,npix_region))
		positions_masked[0] = pixel_positions[0][mask == 1]
		positions_masked[1] = pixel_positions[1][mask == 1]
		if (debug):
			print "Number of pixels: " + str(np.shape(data_masked))

		outputfile = open(outdir+'results_'+str(r)+'_'+str(i)+'.txt', "w")
		outputfile.write("Run " + str(r) + ", region " + str(i))
		outputfile.write('\n')
		np.savetxt(outputfile, ["freq"] + template_names + ["chisq"], fmt="%s", newline=" ")
		outputfile.write('\n')
		outputfile_unc = open(outdir+'results_'+str(r)+'_'+str(i)+'_unc.txt', "w")
		outputfile_unc.write("Run " + str(r) + ", region " + str(i) + " - Uncertainty")
		outputfile_unc.write('\n')
		np.savetxt(outputfile_unc, ["freq"] + template_names + ["chisq"], fmt="%s", newline=" ")
		outputfile_unc.write('\n')
		###
		# Start loop over data maps
		###
		if (debug):
			print ''
			print '*** Looping over maps'
		for j in range(0,num_maps):
			if (debug):
				print ''
				print '**** Map ' + str(j)

			###
			# Define the CMB covariance matrix
			###
			# Just an identity matrix for now
			covar = np.identity(npix_region)

			# Covariance matrix containing noise information
			for i_cov in range(0,npix_region):
				covar[i_cov,i_cov]=np.sqrt(variance_masked[j][i_cov])

			# CMB covariance matrix - want this where we're not subtracting the CMB.
			if (cmb_use_covar == 1):
				if (debug):
					print ''
					print '**** Calculating CMB covariance matrix'
				cmb_covar = np.zeros((npix_region, npix_region))
				# Calculate all of the separations between pixels in one big array
				pos_coord = SkyCoord(frame="galactic", l=positions_masked[1][:]*(180.0/const['pi']), b=90.0-(positions_masked[0][:]*180.0/const['pi']),unit="deg")
				separations = np.array([np.cos(pos_coord[i_temp].separation(pos_coord).rad) for i_temp in range(npix_region)])
				separations = np.reshape(separations, (npix_region, npix_region))
				# print np.shape(separations)

				for i_cov in range(0,npix_region):
					for j_cov in range (0,npix_region):
						Pl = special.lpmv(0, lrange, separations[i_cov,j_cov])
						# print lrange[10]
						# print separations[i_cov,j_cov]
						# print Pl[10]
						cmb_covar[i_cov,j_cov] = np.sum((2.0*lrange+1.0)*cmbspectrum*(beam_windowfunctionsq)*(pixel_windowfunctionsq)*Pl)
						# print i,j,cmb_covar[i,j]

				cmb_covar /= (4.0*const['pi'])
				covar += cmb_covar
				del separations # We don't need this any more.
				del cmb_covar

			# Invert the covariance matrix
			if (debug):
				print ''
				print '**** Inverting covariance matrix'
			covar_inv = np.linalg.inv(covar)

			###
			# Calculate the coefficients
			###
			if (debug):
				print ''
				print '**** Calculating coefficients'

			sigma = (templates_masked.dot(covar_inv)).dot(templates_masked.T) # Bottom half of equation
			sigma_inv = np.linalg.inv(sigma) # Since we're dividing
			a[r][i][j] = ((covar_inv.dot(data_masked[j])).dot(templates_masked.T)).dot(sigma_inv) # Top half of equation

			# Calculate the difference map between the test map and the templates with coefficients
			diff = data_masked[j] - a[r][i][j].dot(templates_masked)
			if (simulation == 0):
				diffmap = (maps[j] - a[r][i][j].dot(templates)) * mask
				hp.write_map(outdir + "diff_"+str(i)+"_"+str(j)+".fits", diffmap)
				origmap = (maps[j]) * mask
				hp.write_map(outdir + "orig_"+str(i)+"_"+str(j)+".fits", origmap)

			chisq[r][i][j] = (diff.dot(covar_inv)).dot(diff.T)

			abar[r][i][j] = sum(a[r][i][j].dot(sigma_inv)) / sum(sigma_inv)

			a_err[r][i][j] = np.sqrt(np.diag(sigma_inv))

			###
			# Output the results
			###
			print str(["freq"] + template_names + ["chisq"])
			print a[r][i][j]
			# print abar[r][i][j]
			print a_err[r][i][j]
			print chisq[r][i][j]
			# for t in range(0,num_templates):
			outputfile.write("%.2f val %s %.4g\n" % (data_frequencies[j], str(a[r][i][j])[1:-1], chisq[r][i][j]))
			outputfile_unc.write("%.2f err %s %.4g\n" % (data_frequencies[j], str(a_err[r][i][j])[1:-1], chisq[r][i][j]))
			# np.savetxt(outputfile, (data_frequencies[j], a_err[j], chisq[j]), fmt="%.2f", newline=" ")

		outputfile.close()
		outputfile_unc.close()


if (num_runs > 1):
	for i in range(0,num_regions):
		for j in range(0,num_maps):
			# Let's do some statistics
			print "i: "+str(i) + ", j: " + str(j)
			print 'Sim: ' + str(simulation_values)
			print 'Noise: ' + str(simulation_noise[0])
			print 'Noise/sqrt(Npix): ' + str(simulation_noise[0]/np.sqrt(num_region_pixels[i]))
			print template_names
			print 'Average: ' + str(np.mean(a, axis=0))
			# print 'Abar avg: ' + str(np.mean(abar, axis=0))
			print 'std: ' + str(np.std(a, axis=0))
			print 'avg err: ' + str(np.mean(a_err, axis=0))

			for k in range(0,num_templates):
				plt.hist(a[:,i,j,k])
				plt.title("Histogram")
				plt.xlabel("Value")
				plt.ylabel("Frequency")
				plt.savefig(outdir+"hist_region_"+str(i)+"_map_"+str(j)+"_template_"+str(k)+".png")
				plt.close()

			# print '#	Input	std	Avg	std:'
			# for k in range(0, num_templates):
			# 	print str(k) + "	" + str(simulation_values[k]) + "	" + str([0])
			# 	print a.shape
			# 	print  "	" + str(np.mean(a[:][i][j][k]))
			# 	print  "	" + str(np.std(a[:][i][j][k]))


# That's all, folks!
# EOF