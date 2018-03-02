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

directory = 'cmb_dx12v3/'
plot = True # Change to True to output plots of the initial maps, False to not do so

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
common_mask = hp.read_map(directory+'dx12_v3_common_mask_int_005a_2048.fits')
smica_mask_raw = hp.read_map(directory+'dx12_v3_smica_mask_int_raw.fits')

nummaps = len(commander)
numpix = len(commander[0])

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
		hp.mollview(smica_mask[i], xsize=4000)
		plt.savefig(directory+'mask_smica_processing_'+str(i)+'.pdf')
		hp.mollview(common_mask[i], xsize=4000)
		plt.savefig(directory+'mask_common_'+str(i)+'.pdf')
		hp.mollview(smica_mask_raw[i], xsize=4000)
		plt.savefig(directory+'mask_smica_raw_'+str(i)+'.pdf')
		# Difference plots
		hp.mollview(commander[i]-nilc[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_nilc_'+str(i)+'.pdf')
		hp.mollview(commander[i]-sevem[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_sevem_'+str(i)+'.pdf')
		hp.mollview(commander[i]-smica[i],min=-0.002,max=0.002, xsize=2000)
		plt.savefig(directory+'commander_sub_smica_'+str(i)+'.pdf')
		# M31 plots
		hp.gnomview(commander[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_commander_'+str(i)+'.pdf')
		hp.gnomview(nilc[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_nilc_'+str(i)+'.pdf')
		hp.gnomview(sevem[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_sevem_'+str(i)+'.pdf')
		hp.gnomview(smica[i],rot=[10.68, 41.27],coord=['G','C'])
		plt.savefig(directory+'m31_smica_'+str(i)+'.pdf')

## Create the mask
mask = np.ones((nummaps, numpix))
# Straight cuts on the maximum in the difference
cut = 0.0001
for i in range(0,3):
	mask[i][commander[i]-nilc[i] > cut] = 0
	mask[i][commander[i]-nilc[i] < -cut] = 0
	mask[i][commander[i]-sevem[i] > cut] = 0
	mask[i][commander[i]-sevem[i] < -cut] = 0
	mask[i][commander[i]-smica[i] > cut] = 0
	mask[i][commander[i]-smica[i] < -cut] = 0
	# and on the noise levels
	mask[i][commander_noise[i] > cut] = 0
	mask[i][commander_noise[i] < -cut] = 0
	mask[i][nilc_noise[i] > cut] = 0
	mask[i][nilc_noise[i] < -cut] = 0
	mask[i][sevem_noise[i] > cut] = 0
	mask[i][sevem_noise[i] < -cut] = 0
	mask[i][smica_noise[i] > cut] = 0
	mask[i][smica_noise[i] < -cut] = 0

	mask[i][smica_mask[i] == 0] = 0

	hp.mollview(mask[i], xsize=4000)
	plt.savefig(directory+'mask_inpainting_'+str(i)+'.pdf')
	print str(100.0 - np.sum(mask[i])/len(mask[i]) * 100) + "% masked in the inpainting mask " + str(i)
	print str(100.0 - np.sum(smica_mask[i])/len(smica_mask[i]) * 100) + "% masked in the smica mask " + str(i)
	print str(100.0 - np.sum(common_mask[i])/len(common_mask[i]) * 100) + "% masked in the common mask"
	print str(100.0 - np.sum(smica_mask_raw[i])/len(smica_mask_raw[i]) * 100) + "% masked in the smica mask"

## Output the mask
hp.write_map(directory+"mask_inpainting.fits", mask,overwrite=True)


## Plot the masked maps
for i in range(0,3):
	# Main maps
	hp.mollview(commander[i]*mask[i],min=-0.0003,max=0.0003, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(nilc[i]*mask[i],min=-0.0003,max=0.0003, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048_'+str(i)+'.pdf')
	hp.mollview(sevem[i]*mask[i],min=-0.0003,max=0.0003, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048_'+str(i)+'pdf')
	hp.mollview(smica[i]*mask[i],min=-0.0003,max=0.0003, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048_'+str(i)+'.pdf')
	# Noise maps
	hp.mollview(commander_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(nilc_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(sevem_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048_noise_'+str(i)+'.pdf')
	hp.mollview(smica_noise[i]*mask[i],min=-0.00002,max=0.00002, xsize=4000)
	plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048_noise_'+str(i)+'.pdf')
	# M31
	hp.gnomview(commander[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+"m31_commander_mask_'+str(i)+'.pdf")
	hp.gnomview(nilc[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+"m31_nilc_mask_'+str(i)+'.pdf")
	hp.gnomview(sevem[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+"m31_sevem_mask_'+str(i)+'.pdf")
	hp.gnomview(smica[i]*mask[i],rot=[10.68, 41.27],coord=['G','C'])
	plt.savefig(directory+"m31_smica_mask_'+str(i)+'.pdf")

# EOF