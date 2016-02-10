#!/usr/bin/env python
# Read in a map, smooth it, and write it out
#
# History:
# Mike Peel   10-Jan-2016   Initial version.
#
# Requirements:
# Numpy, healpy, matplotlib

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def smoothmap(input, output, fwhm_arcmin, pol=False,nside_out=0):
	if (pol):
		field = (0,1,2)
	else:
		field = 0
	map = hp.read_map(input, field=field)
	nside = hp.get_nside(map)
	smoothed_map = hp.sphtfunc.smoothing(map, fwhm=np.radians(fwhm_arcmin/60.0))

	if (nside_out == 0):
		nside_out = nside

	if nside_out != nside:
		map = hp.ud_grade(map, nside_out)

	hp.write_map(map)
