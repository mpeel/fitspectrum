#!/usr/bin/env python
# Test program for template fitting
# Currently set up for magnetic dust models
# MP, 22 September 2016
# MP, 9 November 2016, add noise, masks

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy import special
from spectra import *
from astroutils import *
from tplfit import *
import scipy.optimize as op

def chisq_minimize_magdust(params, npix, planck857, iris5000, freqs, polfrac, covar, maps, covar_inv):

	# Make new templates to compare to the data
	template1 = np.zeros((len(freqs), 2, npix))
	template2 = np.zeros((len(freqs), 2, npix))
	for i in range(0,len(freqs)):
		template1[i][0] = planck857*(freqs[i]/857.0)**params[1]
		template1[i][1] = params[0]*planck857*(freqs[i]/857.0)**params[1]
		template2[i][0] = iris5000*(freqs[i]/857.0)**params[1]
		template2[i][1] = polfrac[i]*iris5000*(freqs[i]/857.0)**params[1]

	template1.shape = len(freqs) * 2 * npix
	template2.shape = len(freqs) * 2 * npix
	templates = np.array([template1,template2])

	# Do the template fitting
	a, a_err, chisq = templatefit_givencovarinv(covar, templates, maps, covar_inv)

	# Return the chisq
	return chisq

# Output directory
outdir = 'magdust/'
ensure_dir(outdir)

# Read in magnetic dust model
nu, magdust_polfrac1, magdust_polfrac2, magdust_polfrac3, magdust_polfrac4, magdust_polfrac5 = np.loadtxt('data/magdust/Fe_pol.dat',dtype=float,usecols=(0,1,2,3,4,5),comments="#",unpack=True)

# Plot it
plt.xscale('log')
# plt.yscale('log')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Polarisation fraction')
plt.plot(nu, magdust_polfrac1,label='0')
plt.plot(nu, magdust_polfrac2,label='0.25')
plt.plot(nu, magdust_polfrac3,label='0.5')
plt.plot(nu, magdust_polfrac4,label='0.75')
plt.plot(nu, magdust_polfrac5,label='1.0')
plt.legend()
plt.savefig(outdir+"test_plotmagneticdust.png")

# Read in a couple of dust maps, and ud_grade them
nside = 64
npix = 12*nside**2
pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
planck857 = hp.read_map('data/512_60.00smoothed_HFI_SkyMap_857_2048_R2.00_full.fits')
planck857 = hp.ud_grade(planck857, nside)
iris5000 = hp.read_map('data/512_60.00smoothed_iris_4995.0_512.fits')
iris5000 = hp.ud_grade(iris5000, nside)

# Do a quick rescaling of the iras 5000 map so it has roughly the same amplitude as the planck 857 map
iris5000 = iris5000 * max(planck857)/max(iris5000)
hp.write_map(outdir+"iris5000_map.fits",iris5000)
hp.write_map(outdir+"planck857_map.fits",planck857)

# Get a Galactic plane mask
galmasksize = 10.0
# mask = galacticmask(nside, galmasksize)
mask = healpixmask(nside, 299-(15.0/2.0), 299+(15.0/2.0), -16.5-6, -16.5+6)
hp.write_map(outdir+"mask.fits",mask)
nregion = int(np.sum(mask))
# iris5000 = planck857

# Now let's create some temperature and polarisation maps, assuming a percentage pol for normal dust, and an amplitude rescaling for magnetic dust
dustpolfrac = 0.2
magdustamp = 0.05
thermaldustamp = 1.0
freqs = [143, 217, 353]
powerlaw = 1.5
tempnoise = [0.55, 0.78, 2.56] / pixarea # May be a factor of 1000 too high!
polnoise = [1.17, 1.75, 7.31] / pixarea # May be a factor of 1000 too high!
polnoise /= 1000.0
tempnoise /= 1000.0

maps = np.zeros((len(freqs), 2, npix))
maskmap = np.zeros((len(freqs), 2, nregion))
noise = np.zeros((len(freqs), 2, nregion))
polfrac = np.zeros(len(freqs))
for i in range(0,len(freqs)):
	polfrac[i] = np.interp(freqs[i], nu, magdust_polfrac2)
	maps[i][0] = thermaldustamp * planck857*(freqs[i]/857.0)**powerlaw + magdustamp*iris5000*(freqs[i]/857.0)**powerlaw + np.random.normal(0.0, tempnoise[i], len(planck857))
	maps[i][1] = thermaldustamp * dustpolfrac*planck857*(freqs[i]/857.0)**powerlaw + polfrac[i]*magdustamp*iris5000*(freqs[i]/857.0)**powerlaw + np.random.normal(0.0, polnoise[i], len(planck857))
	hp.write_map(outdir+"map_"+str(freqs[i])+".fits",maps[i][0])
	hp.write_map(outdir+"map_"+str(freqs[i])+"_pol.fits",maps[i][1])

	noise[i][0] = tempnoise[i] * np.ones(nregion)
	noise[i][1] = polnoise[i] * np.ones(nregion)
	maskmap[i][0] = maps[i][0][mask == 1]
	maskmap[i][1] = maps[i][1][mask == 1]


# Make some templates to go with the maps
template1 = np.zeros((len(freqs), 2, nregion))
template2 = np.zeros((len(freqs), 2, nregion))
for i in range(0,len(freqs)):
	template1[i][0] = planck857[mask == 1]*(freqs[i]/857.0)**powerlaw
	template1[i][1] = dustpolfrac*planck857[mask == 1]*(freqs[i]/857.0)**powerlaw
	template2[i][0] = iris5000[mask == 1]*(freqs[i]/857.0)**powerlaw
	template2[i][1] = polfrac[i]*iris5000[mask == 1]*(freqs[i]/857.0)**powerlaw


# Turn the maps into a dataset for template fitting
maps.shape = len(freqs) * 2 * npix
maskmap.shape = len(freqs) * 2 * nregion
noise.shape = len(freqs) * 2 * nregion
# And the same with the templates
template1.shape = len(freqs) * 2 * nregion
template2.shape = len(freqs) * 2 * nregion
templates = np.array([template1,template2])
# print templates.shape
# print maps.shape

# Create a covariance array - just an identity matrix
covar = np.identity(len(maskmap)) * noise
covar_inv = np.linalg.inv(covar)

# Do the template fitting
# a, a_err, chisq = templatefit_givencovarinv(covar, templates, maps,covar_inv)
# print a
# print a_err
# print chisq

# Parameters - dustpolfrac, powerlaw
params = [00, 1.0]
limits = ((0.0, 1.0), (0.0, 2.0))

chisq_minimize_magdust_fun = lambda *args: chisq_minimize_magdust(*args)
result = op.minimize(chisq_minimize_magdust_fun, params, args=(nregion, planck857[mask==1], iris5000[mask==1], freqs, polfrac, covar, maskmap, covar_inv), bounds=limits, method='L-BFGS-B', options={'maxiter': 40})

print result
print result.x
print result.jac
params = result.x

# And run again to get the template coefficients
last_template1 = np.zeros((len(freqs), 2, nregion))
last_template2 = np.zeros((len(freqs), 2, nregion))
for i in range(0,len(freqs)):
	last_template1[i][0] = planck857[mask==1]*(freqs[i]/857.0)**params[1]
	last_template1[i][1] = params[0]*planck857[mask==1]*(freqs[i]/857.0)**params[1]
	last_template2[i][0] = iris5000[mask==1]*(freqs[i]/857.0)**params[1]
	last_template2[i][1] = polfrac[i]*iris5000[mask==1]*(freqs[i]/857.0)**params[1]

last_template1.shape = len(freqs) * 2 * nregion
last_template2.shape = len(freqs) * 2 * nregion
last_templates = np.array([last_template1,last_template2])
a, a_err, chisq = templatefit_givencovarinv(covar, last_templates, maskmap,covar_inv)

print 'Conclusion: parameter values were:'
print 'Amplitude of thermal dust map: ' + str(a[0]) + ' +- ' + str(a_err[0]) + ' - input was ' + str(thermaldustamp)
print 'Amplitude of magnetic dust map: ' + str(a[1]) + ' +- ' + str(a_err[1]) + ' - input was ' + str(magdustamp)
print 'Thermal dust polarisation fraction: ' + str(params[0]) + ' +- ' + str(result.jac[0]) + ' - input was ' + str(dustpolfrac)
print 'Index of thermal dust emission: ' + str(params[1]) + ' +- ' + str(result.jac[1]) + ' - input was ' + str(powerlaw)
print 'Noise level input was: ' + str(tempnoise) + ', ' + str(polnoise)
print 'Chisq of fit: ' + str(chisq)
print 'Galactic mask size: ' + str(galmasksize) + ' degrees'
print 'Fitting ' + str(nregion) + ' pixels (nside = ' + str(nside) + ')'
print 'Plots are in the "'+outdir+'" folder.'
