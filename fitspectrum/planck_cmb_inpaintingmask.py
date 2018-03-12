#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Plot the maps
# 
# Version history:
#
# 01-Mar-2018  M. Peel       Started

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astroutils import *
from smoothmap import smoothmap

directory = 'cmb_dx12v3/'
plot = True # Change to True to output plots of the initial maps, False to not do so
smooth = 120.0 # arcmin, just for the pol maps, set to False to disable


## Read in the maps, and have a look at them
commander = hp.read_map(directory+'dx12_v3_commander_cmb_005a_2048.fits',field=None)
commander_noise = hp.read_map(directory+'dx12_v3_commander_cmb_hmhd_005a_2048.fits',field=None)
nilc = hp.read_map(directory+'dx12_v3_nilc_cmb_005a_2048.fits',field=None)
nilc_noise = hp.read_map(directory+'dx12_v3_nilc_cmb_hmhd_005a_2048.fits',field=None)
sevem = hp.read_map(directory+'dx12_v3_sevem_cmb_005a_2048.fits',field=None)
sevem_noise = hp.read_map(directory+'dx12_v3_sevem_cmb_hmhd_005a_2048.fits',field=None)
smica = hp.read_map(directory+'dx12_v3_smica_cmb_005a_2048.fits',field=None)
smica_noise = hp.read_map(directory+'dx12_v3_smica_cmb_hmhd_005a_2048.fits',field=None)
smica_mask = hp.read_map(directory+'dx12_v3_smica_mask_int_processing_raw.fits')
smica_mask_raw = hp.read_map(directory+'dx12_v3_smica_mask_int_raw.fits')
smica_mask_pol = hp.read_map(directory+'dx12_v3_smica_mask_pol_processing_raw.fits')
common_mask = hp.read_map(directory+'dx12_v3_common_mask_int_005a_2048.fits')

nummaps = len(commander)
numpix = len(commander[0])
print numpix
nside = hp.npix2nside(numpix)

# Do we want to smooth the maps?
if smooth:
	# Do the smoothing
	smoothmap(directory,directory,'dx12_v3_commander_cmb_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_commander_cmb_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_nilc_cmb_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_nilc_cmb_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_sevem_cmb_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_sevem_cmb_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_smica_cmb_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_smica_cmb_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_commander_cmb_hmhd_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_commander_cmb_hmhd_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_nilc_cmb_hmhd_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_nilc_cmb_hmhd_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_sevem_cmb_hmhd_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_sevem_cmb_hmhd_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	smoothmap(directory,directory,'dx12_v3_smica_cmb_hmhd_005a_2048.fits',str(smooth)+'smoothed_dx12_v3_smica_cmb_hmhd_005a_2048.fits', np.sqrt(smooth**2-(5.0)**2),nside_out=nside,usehealpixfits=True)
	# Read in the smoothed maps
	commander_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_commander_cmb_005a_2048.fits',field=None)
	nilc_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_nilc_cmb_005a_2048.fits',field=None)
	sevem_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_sevem_cmb_005a_2048.fits',field=None)
	smica_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_smica_cmb_005a_2048.fits',field=None)
	commander_noise_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_commander_cmb_hmhd_005a_2048.fits',field=None)
	nilc_noise_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_nilc_cmb_hmhd_005a_2048.fits',field=None)
	sevem_noise_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_sevem_cmb_hmhd_005a_2048.fits',field=None)
	smica_noise_smooth = hp.read_map(directory+str(smooth)+'smoothed_dx12_v3_smica_cmb_hmhd_005a_2048.fits',field=None)

if plot:
	for i in range(0,3):
		# Overall plots
		hp.mollview(commander[i],min=-0.0003,max=0.0003, xsize=2000)
		plt.savefig(directory+'dx12_v3_commander_cmb_005a_2048_'+str(i)+'.pdf')
		hp.mollview(commander_noise[i],min=-0.00002,max=0.00002, xsize=2000)
		plt.savefig(directory+'dx12_v3_commander_cmb_hmhd_005a_2048_'+str(i)+'.pdf')
		hp.mollview(nilc[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'dx12_v3_nilc_cmb_005a_2048_'+str(i)+'.pdf')
		hp.mollview(nilc_noise[i],min=-0.0002,max=0.0002, xsize=2000)
		plt.savefig(directory+'dx12_v3_nilc_cmb_hmhd_005a_2048_'+str(i)+'.pdf')
		hp.mollview(sevem[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'dx12_v3_sevem_cmb_005a_2048_'+str(i)+'.pdf')
		hp.mollview(sevem_noise[i],min=-0.0002,max=0.0002, xsize=2000)
		plt.savefig(directory+'dx12_v3_sevem_cmb_hmhd_005a_2048_'+str(i)+'.pdf')
		hp.mollview(smica[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'dx12_v3_smica_cmb_005a_2048_'+str(i)+'.pdf')
		hp.mollview(smica_noise[i],min=-0.0002,max=0.0002, xsize=2000)
		plt.savefig(directory+'dx12_v3_smica_cmb_hmhd_005a_2048_'+str(i)+'.pdf')
		plt.close('all')
		if i == 0:
			hp.mollview(common_mask, xsize=4000)
			plt.savefig(directory+'mask_common_'+str(i)+'.pdf')
			hp.mollview(smica_mask, xsize=4000)
			plt.savefig(directory+'mask_smica_processing_'+str(i)+'.pdf')
			hp.mollview(smica_mask_raw, xsize=4000)
			plt.savefig(directory+'mask_smica_raw_'+str(i)+'.pdf')
		else:
			hp.mollview(smica_mask_pol, xsize=4000)
			plt.savefig(directory+'mask_smica_processing_'+str(i)+'.pdf')
		# Difference plots
		hp.mollview(commander[i]-nilc[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_nilc_'+str(i)+'.pdf')
		hp.mollview(commander[i]-sevem[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_sevem_'+str(i)+'.pdf')
		hp.mollview(commander[i]-smica[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_smica_'+str(i)+'.pdf')
		plt.close('all')
		# M31 plots
		hp.gnomview(commander[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_commander_'+str(i)+'.pdf')
		hp.gnomview(nilc[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_nilc_'+str(i)+'.pdf')
		hp.gnomview(sevem[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_sevem_'+str(i)+'.pdf')
		hp.gnomview(smica[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_smica_'+str(i)+'.pdf')
		plt.close('all')

## Create the mask
mask = np.ones((nummaps+1, numpix))
# Straight cuts on the maximum in the difference
cut = [0.0001, 0.0000015, 0.0000015]
for i in range(0,3):
	if i == 0 or smooth == False:
		print str(i) + "raw"
		mask[i][commander[i]-nilc[i] > cut[i]] = 0
		mask[i][commander[i]-nilc[i] < -cut[i]] = 0
		mask[i][commander[i]-sevem[i] > cut[i]] = 0
		mask[i][commander[i]-sevem[i] < -cut[i]] = 0
		mask[i][commander[i]-smica[i] > cut[i]] = 0
		mask[i][commander[i]-smica[i] < -cut[i]] = 0
		# and on the noise levels
		mask[i][commander_noise[i] > cut[i]] = 0
		mask[i][commander_noise[i] < -cut[i]] = 0
		mask[i][nilc_noise[i] > cut[i]] = 0
		mask[i][nilc_noise[i] < -cut[i]] = 0
		mask[i][sevem_noise[i] > cut[i]] = 0
		mask[i][sevem_noise[i] < -cut[i]] = 0
		mask[i][smica_noise[i] > cut[i]] = 0
		mask[i][smica_noise[i] < -cut[i]] = 0
	else:
		print str(i) + "smooth"
		mask[i][commander_smooth[i]-nilc_smooth[i] > cut[i]] = 0
		mask[i][commander_smooth[i]-nilc_smooth[i] < -cut[i]] = 0
		mask[i][commander_smooth[i]-sevem_smooth[i] > cut[i]] = 0
		mask[i][commander_smooth[i]-sevem_smooth[i] < -cut[i]] = 0
		mask[i][commander_smooth[i]-smica_smooth[i] > cut[i]] = 0
		mask[i][commander_smooth[i]-smica_smooth[i] < -cut[i]] = 0
		# and on the noise levels
		mask[i][commander_noise_smooth[i] > cut[i]] = 0
		mask[i][commander_noise_smooth[i] < -cut[i]] = 0
		mask[i][nilc_noise_smooth[i] > cut[i]] = 0
		mask[i][nilc_noise_smooth[i] < -cut[i]] = 0
		mask[i][sevem_noise_smooth[i] > cut[i]] = 0
		mask[i][sevem_noise_smooth[i] < -cut[i]] = 0
		mask[i][smica_noise_smooth[i] > cut[i]] = 0
		mask[i][smica_noise_smooth[i] < -cut[i]] = 0


	if i == 0:
		mask[i][smica_mask == 0] = 0
	else:
		mask[i][smica_mask_pol == 0] = 0

# Also combine the masks
mask[3][mask[1] == 0] = 0
mask[3][mask[2] == 0] = 0

for i in range(0,4):
	hp.mollview(mask[i], xsize=4000)
	plt.savefig(directory+'mask_inpainting_'+str(i)+'.pdf')
	print str(100.0 - np.sum(mask[i])/len(mask[i]) * 100) + "% masked in the inpainting mask " + str(i)
	if i == 0:
		print str(100.0 - np.sum(common_mask)/len(common_mask) * 100) + "% masked in the common mask"
		print str(100.0 - np.sum(smica_mask)/len(smica_mask) * 100) + "% masked in the smica mask " + str(i)
		print str(100.0 - np.sum(smica_mask_raw)/len(smica_mask_raw) * 100) + "% masked in the smica mask"
	else:
		print str(100.0 - np.sum(smica_mask_pol)/len(smica_mask_pol) * 100) + "% masked in the smica mask " + str(i)

## Output the mask
hp.write_map(directory+"mask_inpainting.fits", mask,overwrite=True)
plt.close('all')


## Plot the masked maps
for i in range(0,3):
	# Main maps
	hp.mollview(commander[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(nilc[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(sevem[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(smica[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048_'+str(i)+'.pdf')
	plt.close('all')
	hp.mollview(commander_smooth[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_smooth_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(nilc_smooth[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_smooth_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(sevem_smooth[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_smooth_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(smica_smooth[i]*mask[i], xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_smooth_cmb_005a_2048_'+str(i)+'.pdf')
	plt.close('all')
	if i > 0:
		hp.mollview(commander_smooth[i]*mask[3], xsize=4000)
		plt.savefig(directory+'masked_dx12_v3_commander_smooth_cmb_005a_2048_'+str(i)+'_3.pdf')
		hp.mollview(nilc_smooth[i]*mask[3], xsize=4000)
		plt.savefig(directory+'masked_dx12_v3_nilc_smooth_cmb_005a_2048_'+str(i)+'_3.pdf')
		hp.mollview(sevem_smooth[i]*mask[3], xsize=4000)
		plt.savefig(directory+'masked_dx12_v3_sevem_smooth_cmb_005a_2048_'+str(i)+'_3.pdf')
		hp.mollview(smica_smooth[i]*mask[3], xsize=4000)
		plt.savefig(directory+'masked_dx12_v3_smica_smooth_cmb_005a_2048_'+str(i)+'_3.pdf')
		plt.close('all')

	# Noise maps
	hp.mollview(commander_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(nilc_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(sevem_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(smica_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(commander_noise_smooth[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_smooth_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(nilc_noise_smooth[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_smooth_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(sevem_noise_smooth[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_smooth_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(smica_noise_smooth[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_smooth_cmb_005a_2048_noise_'+str(i)+'.pdf')
	plt.close('all')
	# M31
	hp.gnomview(commander[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+'m31_commander_mask_'+str(i)+'.pdf')
	hp.gnomview(nilc[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+'m31_nilc_mask_'+str(i)+'.pdf')
	hp.gnomview(sevem[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+'m31_sevem_mask_'+str(i)+'.pdf')
	hp.gnomview(smica[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+'m31_smica_mask_'+str(i)+'.pdf')
	plt.close('all')

# EOF