#!/usr/bin/env python
# Read in the Planck high-frequency data,
# smooth to a common resolution,
# and subtract a model to look for residual emission
#
# History:
# Mike Peel   18-Jan-2016   Initial version.
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from spectra import *
from astroutils import *

###
# Configuration options
###
indir = 'PR2/'
outdir = "PR2_subtractmap/"
map_filenames = ['512_60.00smoothed_HFI_SkyMap_353_2048_PR2_full.fits', '512_60.00smoothed_HFI_SkyMap_217_2048_PR2_full.fits', '512_60.00smoothed_HFI_SkyMap_143_2048_PR2_full.fits', '512_60.00smoothed_HFI_SkyMap_100_2048_PR2_full.fits', '512_60.00smoothed_LFI_SkyMap_70_256_PR2_full.fits', '512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits', '512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits']
freqs = [353.0, 217.0, 143.0, 100.0, 70.0, 44.0, 30.0]
nside = 512
beta_dust = 1.6
t_dust = 20.0
beta_sync = -3.0
const = get_spectrum_constants()
solid_angle = 1.0

# Read in the maps
num_maps = len(map_filenames)
npix = hp.nside2npix(nside)
maps = np.zeros((num_maps, 3, npix))
for i,filenames in enumerate(map_filenames):
	map = hp.read_map(indir+map_filenames[i],(0,1,2))
	test = np.array(map)
	maps[i][:][:] = test[:]

# Check the output directory exists
ensure_dir(outdir)

# Calculate the amplitudes for the components to subtract
amp_353 = thermaldust(const, freqs[0], 1.0e-5, beta_dust, t_dust, freqs[0], solid_angle)
amp_30 = synchrotron(const, freqs[6], freqs[6], 1.0, beta_sync)

map_353 = maps[0]
map_30 = maps[6]
for i in range(0,num_maps):
	# Calculate the spectrum of the components at this frequency
	amp = thermaldust(const, freqs[i], 1.0e-5, beta_dust, t_dust, freqs[0], solid_angle)
	sync_amp = synchrotron(const, freqs[i], freqs[6], 1.0, beta_sync)
	# Subtract them off using the amplitudes calculated before.
	resid = maps[i] - map_353 * (amp/amp_353) - map_30 * (sync_amp/amp_30)
	# Write out the map
	hp.write_map(outdir+"subtract_dust_sync_planck_"+str(freqs[i])+".fits", (resid[0], resid[1], resid[2]))

