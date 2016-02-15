# Fit an SED to given multifrequency data
# 
# Version history:
# IDL batch script for running haperflux on several maps
# and plotting/fitting spectra for Planck and ancillary data
#
#
# 28-Aug-2010  C. Dickinson  Added model fitting  
# 02-Sep-2010  C. Dickinson  Added inputs
# 20-Oct-2010  M. Peel		 Add different input map unit options (via filename),
#							 write out data file of results
# 09-Sep-2011  M. Peel/C. Dickinson Updated from Mike's version.
#
#                         Differences are
# 1. Commented out new parameters to haperflux
# 2. Added noise model to haperflux (not as input yet)
# 3. Added /nocmb option
# 4. Changed range of LFI frequencies (28.5 to 28.4 so it works with
# this frequency!)
# 5. Increased error bar thickness and plotsymbol size
# 6. Added extra oplot so plotsymbol can be seen
# 7. Addded extra IF statements to WMAP colour corrections in case no
# good data there.
#
# 22-Sep-2011 M. Peel  Output grep-able summary of fitted parameters; epstopdf call; input for noisemodel
# 10-Nov-2011 M. Peel  Adding code to deal with halpha and RRL maps
# 22-Dec-2011 C. Dickinson  Added minfreq and maxfreq options
# 04-Jan-2012 C. Dickinson  Added amemodels option. Residual plots now
#                           include modelling errors
# 08-Jan-2012 C. Dickinson  Added modelling errors when fitting for Te
# 12-Jan-2012 C. Dickinson  Added dust_optical_depth_freq for thermal dust
# 22-Jan-2012 C. Dickinson  Added usecochannels option
# 24-Jan-2012 M. Peel       Add option to read in flux densities from a text file.
# 01-Mar-2012 M. Peel       Also output the flux densities prior to colour correction.
# 18-Mar-2012 C. Dickinson  Added extra columns to .dat file
# 19-Mar-2012 M. Peel       Additional formatting in .dat file
# 16-Apr-2012 C. Dickinson  Added /shiftspdust option
# 18-Jul-2012 C. Dickinson  Added /nothermaldust option + minor bugs
# 04-Dec-2012 M. Peel		Added /nosyserr option + set s25000,s12000 to 0 when no map at that freq exists
# 05-Dec-2012 M. Peel       Fixed bug with increasing the uncertainty on CO-affected 100 and 217GHz channels
# 06-Dec-2012 M. Peel       Added /nocoerr option, change default noise_model to 1
# 11-Jan-2013 M. Peel       New wmap colour correction code, and updated LFI colour correction coefficients
# 14-Jan-2013 M. Peel       Updating HFI colour correction code
# 18-Jan-2013 M. Peel       Updated LFI frequency ranges
# 23-Jan-2013 M. Peel       Updated LFI colour corrections, creating planckcc function
# 19-Apr-2013 M. Peel       Update to published colour correction code ('LFI_fastcc' and 'hfi_colour_correction')
# 
# Mike Peel   02-Feb-2016   v1.0 Migrating from IDL runaperflux.pro


import numpy as np
#import scipy as sp
from mpfit import mpfit
from spectra import *
import copy
import matplotlib.pyplot as plt

# Define some constants, used later in the SED functions
const = get_spectrum_constants()

solid_angle = 1.0

def spectrum(p, fjac=None, x=None, y=None, err=None):
	print p
	print 'Data:'
	print y
	print synchrotron(const, x, 1.0, p[7], p[8])
	model = synchrotron(const, x, 1.0, p[7], p[8]) + freefree(const, x, p[0], p[1], solid_angle) + thermaldust(const, x, p[2], p[3], p[4], const['dust_optical_depth_freq'], solid_angle)
	print 'Model:'
	print model
	plt.plot(x, model, 'r', x, y, 'g')
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	status = 0
	return [status, (y-model)/err]

def fitspectrum(filename, indir='', outdir='', spd_file='amemodels/spdust2_wim.dat'):

	solid_angle = 1.0

	# Read in the spinning dust curve
	spd = np.loadtxt(spd_file, dtype=float,comments=';')
	# print spd

	# Read in the file containing the SED data
	# Should be of the format 'freq flux err'
	inputspectrum = np.loadtxt(filename, dtype=float)
	freqs = inputspectrum[:,0] / 1e9
	fd = inputspectrum[:,1]
	fd_err = inputspectrum[:,2]

	###
	# Use MPFit
	###

	# Set up the initial parameters
	nparams = 9
	parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
	parinfo=[]
	for i in range(0,nparams):
		parinfo.append(copy.deepcopy(parbase))

	parinfo[0]['limited'] = [1,1]          # EM
	parinfo[0]['limits'] = [0,1e10]
	parinfo[1]['limited'] = [1,1]          # Te
	parinfo[1]['limits'] = [1000.,30000.]
	parinfo[1]['fixed'] = 1  # fix Te
	parinfo[2]['limited'] = [1,1]          # thermal dust amplitude (optical depth)
	parinfo[2]['limits'] = [0.0,100.] 
	parinfo[3]['limited'] = [1,1]          # thermal dust index
	parinfo[3]['limits'] = [0.5,3.5] 
	parinfo[3]['fixed'] = 0  # fix dust emissivity index
	parinfo[4]['limited'] = [1,1]          # thermal dust temperature
	parinfo[4]['limits'] = [5.,200.] 
	parinfo[5]['limited'] = [1,1]          # spinning dust amplitude
	parinfo[5]['limits'] = [0.,1e25] 
	parinfo[5]['fixed'] = 1  # fix spinningdust_amp
	parinfo[6]['limited'] = [0,0]          # delta T (CMB) in K
	parinfo[6]['limits'] = [-150e-6,150e-6]
	parinfo[7]['limited'] = [1,1]		# sync amplitude
	parinfo[7]['limits'] = [0.0,1e10]
	parinfo[8]['limited'] = [1,1]		# sync index
	parinfo[8]['limits'] = [-3.0,1.0]
	# parinfo[9]['fixed'] = 1  # spdust shift

	p0 = [1.0,8000.,1e-6,1.7,20.,1.,0.,1.0 ,-1.0]  # standard start values   
	fa = {'x':freqs, 'y':fd, 'err':fd_err}
	m = mpfit(spectrum, p0, parinfo=parinfo,functkw=fa)
	if (m.status <= 0): 
		print 'error message = ', m.errmsg
	










