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

###
# Configuration options
###

resolution = 5.0 # arcmin
data = [[]

]




smoothed_map = hp.sphtfunc.smoothalm(map, fwhm=60, arcmin=True)