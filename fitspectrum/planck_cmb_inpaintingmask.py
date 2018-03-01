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

## Read in the maps, and have a look at them
commander = hp.read_map(directory+'dx12_v3_commander_cmb_005a_2048.fits')
hp.mollview(commander,min=-0.0003,max=0.0003, xsize=2000)
plt.savefig(directory+'dx12_v3_commander_cmb_005a_2048.pdf')
commander_noise = hp.read_map(directory+'dx12_v3_commander_cmb_hmhd_005a_2048.fits')
hp.mollview(commander_noise,min=-0.00002,max=0.00002, xsize=2000)
plt.savefig(directory+'dx12_v3_commander_cmb_hmhd_005a_2048.pdf')
exit()
nilc = hp.read_map(directory+'dx12_v3_nilc_cmb_005a_2048.fits')
# hp.mollview(nilc,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'dx12_v3_nilc_cmb_005a_2048.pdf')
nilc_noise = hp.read_map(directory+'dx12_v3_nilc_cmb_hmhd_005a_2048.fits')
# hp.mollview(nilc_noise,min=-0.0002,max=0.0002, xsize=2000)
# plt.savefig(directory+'dx12_v3_nilc_cmb_hmhd_005a_2048.pdf')

sevem = hp.read_map(directory+'dx12_v3_sevem_cmb_005a_2048.fits')
# hp.mollview(sevem,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'dx12_v3_sevem_cmb_005a_2048.pdf')
sevem_noise = hp.read_map(directory+'dx12_v3_sevem_cmb_hmhd_005a_2048.fits')
# hp.mollview(sevem_noise,min=-0.0002,max=0.0002, xsize=2000)
# plt.savefig(directory+'dx12_v3_sevem_cmb_hmhd_005a_2048.pdf')

smica = hp.read_map(directory+'dx12_v3_smica_cmb_005a_2048.fits')
# hp.mollview(smica,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'dx12_v3_smica_cmb_005a_2048.pdf')
smica_noise = hp.read_map(directory+'dx12_v3_smica_cmb_hmhd_005a_2048.fits')
# hp.mollview(smica_noise,min=-0.0002,max=0.0002, xsize=2000)
# plt.savefig(directory+'dx12_v3_smica_cmb_hmhd_005a_2048.pdf')

# hp.mollview(commander-nilc,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'commander_sub_nilc.pdf')

# hp.mollview(commander-sevem,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'commander_sub_sevem.pdf')

# hp.mollview(commander-smica,min=-0.002,max=0.002, xsize=2000)
# plt.savefig(directory+'commander_sub_smica.pdf')


hp.gnomview(commander,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_commander.pdf")
hp.gnomview(nilc,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_nilc.pdf")
hp.gnomview(sevem,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_sevem.pdf")
hp.gnomview(smica,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_smica.pdf")

## Create the mask
mask = np.ones(len(commander))
# Straight cuts on the maximum in the difference
cut = 0.0001
mask[commander-nilc > cut] = 0
mask[commander-nilc < -cut] = 0
mask[commander-sevem > cut] = 0
mask[commander-sevem < -cut] = 0
mask[commander-smica > cut] = 0
mask[commander-smica < -cut] = 0
# and on the noise levels
mask[commander_noise > cut] = 0
mask[commander_noise < -cut] = 0
mask[nilc_noise > cut] = 0
mask[nilc_noise < -cut] = 0
mask[sevem_noise > cut] = 0
mask[sevem_noise < -cut] = 0
mask[smica_noise > cut] = 0
mask[smica_noise < -cut] = 0


## Output the mask
hp.write_map(directory+"mask_inpainting.fits", mask,overwrite=True)
hp.mollview(mask, xsize=4000)
plt.savefig(directory+'mask_inpainting.pdf')
print str(100.0 - np.sum(mask)/len(mask) * 100) + "% masked in the inpainting mask"

# nside = hp.npix2nside(len(mask))
# numaffected = 0
# for i in range(0,len(mask)):
# 	neighbours = hp.get_all_neighbours(nside, i)
# 	if mask[i] == 1:
# 		if np.sum(mask[neighbours]) < 3:
# 			mask[i] = 0
# 			numaffected += 1
# 			print i
# print numaffected
# hp.write_map(directory+"mask_inpainting_infill.fits", mask,overwrite=True)
# hp.mollview(mask, xsize=4000)
# plt.savefig(directory+'mask_inpainting_infill.pdf')
# print str(100.0 - np.sum(mask)/len(mask) * 100) + "% masked in the inpainting (infilled) mask"



## Output masked maps
hp.mollview(commander*mask,min=-0.0003,max=0.0003, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048.pdf')
hp.mollview(nilc*mask,min=-0.0003,max=0.0003, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048.pdf')
hp.mollview(sevem*mask,min=-0.0003,max=0.0003, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048.pdf')
hp.mollview(smica*mask,min=-0.0003,max=0.0003, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048.pdf')

hp.mollview(commander_noise*mask,min=-0.00002,max=0.00002, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_commander_cmb_005a_2048_noise.pdf')
hp.mollview(nilc_noise*mask,min=-0.00002,max=0.00002, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_nilc_cmb_005a_2048_noise.pdf')
hp.mollview(sevem_noise*mask,min=-0.00002,max=0.00002, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_sevem_cmb_005a_2048_noise.pdf')
hp.mollview(smica_noise*mask,min=-0.00002,max=0.00002, xsize=4000)
plt.savefig(directory+'masked_dx12_v3_smica_cmb_005a_2048_noise.pdf')

## Also plot the existing masks
common_mask = hp.read_map(directory+'dx12_v3_common_mask_int_005a_2048.fits')
hp.mollview(common_mask, xsize=4000)
plt.savefig(directory+'mask_common.pdf')
print str(100.0 - np.sum(common_mask)/len(common_mask) * 100) + "% masked in the common mask"

smica_mask = hp.read_map(directory+'dx12_v3_smica_mask_int_raw.fits')
hp.mollview(smica_mask, xsize=4000)
plt.savefig(directory+'mask_smica.pdf')
print str(100.0 - np.sum(smica_mask)/len(smica_mask) * 100) + "% masked in the smica mask"

smica_mask = hp.read_map(directory+'dx12_v3_smica_mask_int_processing_raw.fits')
hp.mollview(smica_mask, xsize=4000)
plt.savefig(directory+'mask_smica_processing.pdf')
print str(100.0 - np.sum(smica_mask)/len(smica_mask) * 100) + "% masked in the smica mask"


hp.gnomview(commander*mask,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_commander_mask.pdf")
hp.gnomview(nilc*mask,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_nilc_mask.pdf")
hp.gnomview(sevem*mask,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_sevem_mask.pdf")
hp.gnomview(smica*mask,rot=[10.68, 41.27],coord=['G','C'])
plt.savefig(directory+"m31_smica_mask.pdf")

# EOF