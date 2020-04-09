# Do a quick comparison of the WMAP K-band and Planck 28.4GHz maps, and generate a mask for the latter based on the difference in the polarisation amplitude.
# MP, 19 March 2020
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Set configuration
basedir = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
wmap_filename = basedir+'512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits'
planck_filename = basedir+'512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits'
index = -3.0
threshold = 0.03
apodise_threshold=0.6

# Get the maps
wmap_map = hp.read_map(wmap_filename,field=None)
wmap = np.sqrt(wmap_map[1]**2+wmap_map[2]**2)

planck_map = hp.read_map(planck_filename,field=None)
planck = np.sqrt(planck_map[1]**2+planck_map[2]**2)

# Rescale Planck to WMAP frequency
planck_rescale = planck.copy()*(22.8/28.4)**index
# Difference them
difference = wmap-planck_rescale

# Plot them
hp.mollview(wmap,max=0.3)
plt.savefig('compare_wmap_planck_polmap_wmap.png')
plt.clf()
hp.mollview(planck_rescale,max=0.3)
plt.savefig('compare_wmap_planck_polmap_planck.png')
plt.clf()
hp.mollview(difference,min=-threshold,max=threshold)
plt.savefig('compare_wmap_planck_polmap_diff.png')
plt.clf()

# Smooth the difference a bit, to reduce noise for the mask
difference_smth = hp.sphtfunc.smoothing(difference, fwhm=1.0*np.pi/180.0)

# Create a mask
mask = np.ones(len(wmap))
mask[np.abs(difference_smth) > threshold] = 0
hp.mollview(mask)
plt.savefig('compare_wmap_planck_polmap_mask.png')
plt.clf()

# Do a bit of smoothing
mask = hp.sphtfunc.smoothing(mask, fwhm=1.0*np.pi/180.0)
hp.mollview(mask)
plt.savefig('compare_wmap_planck_polmap_mask_smth.png')
plt.clf()
mask[mask < apodise_threshold] = 0.0
mask[mask > 0.5*apodise_threshold] = 1.0
hp.mollview(mask)
plt.savefig('compare_wmap_planck_polmap_mask_smth_ap.png')
plt.clf()

# Plot the masked data
difference[mask == 0] = hp.UNSEEN
hp.mollview(difference,min=-threshold,max=threshold)
plt.savefig('compare_wmap_planck_polmap_diff_mask.png')
plt.clf()

# Finally, write out the mask
hp.write_map('compare_wmap_planck_polmap_mask.fits',mask,overwrite=True)


# EOF