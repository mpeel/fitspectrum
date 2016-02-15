import os
import numpy as np
import healpy as hp
from math import pi

# Ensure that a directory exists
# Function from http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
def ensure_dir(f):
	d = os.path.dirname(f)
	if not os.path.exists(d):
		os.makedirs(d)

# Create a healpix mask given an nside, min/max longitudes and latitudes, and a coordinate system
def healpixmask(nside, long_min, long_max, lat_min, lat_max, coordsystem='G'):

	npix = hp.nside2npix(nside)
	masked_map = np.ones(npix)

	for i in range(0,npix):
		pos = hp.pixelfunc.pix2ang(nside, i)
		phi = pos[1]*(180.0/pi)
		theta = 90.0-(pos[0]*180.0/pi)

		# Check for negative longitude ranges
		if (phi <= long_max and phi >= long_min and theta <= lat_max and theta >= lat_min):
			masked_map[i] = 0
		if (long_max < 0): 
			# Assuming that this means that long_min is also less than zero
			if (phi-360.0 <= long_max and phi-360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 0
		if (long_min < 0): 
			# Longitudes higher than 0 will have been dealt with above, deal with the negative ones here.
			if (phi-360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 0


	return masked_map
