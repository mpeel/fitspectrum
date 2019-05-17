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
# Mike Peel   07-Apr-2016   v1.0.1 Expanding, adding emcee


import numpy as np
import scipy.optimize as op
from mpfit import mpfit
from spectra import *
from astroutils import *
import copy
import matplotlib.pyplot as plt
import emcee
import corner

def calc_beam_area(beam1,beam2):
	return (np.pi * (beam1*np.pi/180.0)*(beam2*np.pi/180.0))/(4.0*np.log(2.0))


# Define some constants, used later in the SED functions
const = get_spectrum_constants()
solid_angle = calc_beam_area(100.0/60.0,0.7*100.0/60.0)#1.0e-5 # For now
spd_freq = []
spd_amp = []

def spectrum(p, fjac=None, x=None):
	model = synchrotron(const, x, 1.0, p[7], p[8]) + freefree(const, x, p[0], p[1], solid_angle) + spinningdust(spd_freq, spd_amp, solid_angle, x, p[5], p[9]) + cmb(const, x, p[6]*1e-6, solid_angle) + thermaldust(const, x, p[2]*1e-5, p[3], p[4], const['dust_optical_depth_freq'], solid_angle)# + thermaldust(const, x, p[10]*1e-5, p[11], p[12], const['dust_optical_depth_freq'], solid_angle)

	return model

def mpfitfunction(p, fjac=None, x=None, y=None, err=None):

	model = spectrum(p, fjac, x)

	status = 0
	return [status, (y-model)/err]

def lnlike(p, x, y, yerr):
	model = spectrum(p, x=x)
	inv_sigma2 = 1.0/(yerr**2) # + model**2*np.exp(2*lnf))
	value = -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	return value

def lnprior(p, bounds):
	for i in range (0,len(p)):
		if p[i] < bounds[i][0] or p[i] > bounds[i][1]:
			return -np.inf
	# If we've got to here, we're within the bounds.
	return 0.0

def lnprob(p, bounds, x, y, yerr):
    lp = lnprior(p, bounds)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, x, y, yerr)

def plot_results(outfile,srcname,minfreq,maxfreq,params,const,solid_angle,spd_freq,spd_amp,freqs,fd,fd_err,goodvals,minflux,maxflux):
	
	x = np.arange(minfreq,maxfreq,(maxfreq-minfreq)/1000.0)

	# Generate the model and plot it
	# if nosync == False:
	model_sync = synchrotron(const, x, 1.0, params[7], params[8])
	plt.plot(x, model_sync, 'r')
	# if nofreefree == False:
	model_freefree = freefree(const, x, params[0], params[1], solid_angle)
	plt.plot(x, model_freefree, 'b')
	# if noame == False:
	model_spd = spinningdust(spd_freq, spd_amp, solid_angle, x, params[5], params[9])
	plt.plot(x, model_spd, 'black')
	# if nodust == False:
	model_dust1 = thermaldust(const, x, params[2]*1e-5, params[3], params[4], const['dust_optical_depth_freq'], solid_angle)
	plt.plot(x, model_dust1, 'g')
	# if nodust2 == False:
	# model_dust2 = thermaldust(const, x, params[10]*1e-5, params[11], params[12], const['dust_optical_depth_freq'], solid_angle)
	# plt.plot(x, model_dust2, 'y')

	model_overall = spectrum(params, x=x)
	plt.plot(x, model_overall, 'orange')
	# Add the data to the plot
	plt.errorbar(freqs[goodvals == 1], fd[goodvals == 1], fd_err[goodvals == 1],fmt='+')
	plt.errorbar(freqs[goodvals == 0], fd[goodvals == 0], fd_err[goodvals == 0])

	# Formatting, and output
	plt.title(srcname)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Frequency (GHz)')
	plt.ylabel('Flux density (Jy)')
	plt.ylim(ymin=minflux*0.1,ymax=maxflux*10)
	plt.savefig(outfile)
	plt.close()
	return

def write_results(outfile,srcname,freqs,fd,fd_err,params,perror):
	outputfile = open(outfile, "w")
	outputfile.write("# " + srcname + "\n")
	for i in range(0,len(freqs)):
		outputfile.write(str(freqs[i]) + "	" + str(fd[i]) + "	" + str(fd_err[i]) + "\n")

	outputfile.write('# Sync amplitude: ' + str(params[7]) + " +- " + str(perror[7]) + '\n')
	outputfile.write('# Sync index: ' + str(params[8]) + " +- " + str(perror[8]) + '\n')
	outputfile.write('# Freefree amplitude: ' + str(params[0]) + " +- " + str(perror[0]) + '\n')
	outputfile.write('# Freefree temp: ' + str(params[1]) + " +- " + str(perror[1]) + '\n')
	outputfile.write('# AME amplitude: ' + str(params[5]) + " +- " + str(perror[5]) + '\n')
	outputfile.write('# AME shift: ' + str(params[9]) + " +- " + str(perror[9]) + '\n')
	outputfile.write('# CMB amplitude: ' + str(params[6]*1e-6) + " +- " + str(perror[6]*1e-6) + '\n')
	# outputfile.write('# AME at 25GHz: ' + str(spinningdust(spd_freq, spd_amp, solid_angle, 25.0, m.params[5], m.params[9])) + '\n')
	# outputfile.write('# AME at 30GHz: ' + str(spinningdust(spd_freq, spd_amp, solid_angle, 30.0, m.params[5], m.params[9])) + '\n')
	outputfile.write('# Dust amplitude: ' + str(params[2]*1e-5) + " +- " + str(perror[2]*1e-5) + '\n')
	outputfile.write('# Dust index: ' + str(params[3]) + " +- " + str(perror[3]) + '\n')
	outputfile.write('# Dust temperature: ' + str(params[4]) + " +- " + str(perror[4]) + '\n')
	# outputfile.write('# Dust2 amplitude: ' + str(params[10]) + " +- " + str(perror[10]) + '\n')
	# outputfile.write('# Dust2 index: ' + str(params[11]) + " +- " + str(perror[11]) + '\n')
	# outputfile.write('# Dust2 temperature: ' + str(params[12]) + " +- " + str(perror[12]) + '\n')

	outputfile.close()
	return

def fitspectrum(filename, srcname='',indir='', outdir='', spd_file='amemodels/spdust2_wim.dat',
	format=1, inunits='Jy',fitunits='Jy',minfitfreq=0,maxfitfreq=0,
	nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,
	startparams=0,mcmc=True,quiet=False):#,nodust2=True, fixdust2temp=0

	print(outdir)
	ensure_dir(outdir)

	# Read in the spinning dust curve
	global spd_freq
	global spd_amp
	spd_freq, spd_amp = np.loadtxt(spd_file, dtype=float,usecols=(0,1),comments=';',unpack=True)

	# solid_angle = calc_beam_area(100.0/60.0,0.7*100.0/60.0)
	# print(solid_angle)

	###
	# Read in the file containing the SED data
	# Different formats are supported, namely:
	# 1: freq [GHz] flux [Jy] err [Jy]'
	# 2: freq [Hz] flux [Jy] err [Jy]'
	# 3: wavelength [um]  flux [W/m^2/um] err [W/m^2/um] 
	# Use # to mark comments in the input file, and these will be ignored.
	# For the actual fitting, we'll use freq [GHz], flux [Jy], err [Jy] - other formats will be converted to this.
	###
	inputspectrum = np.loadtxt(filename)
	if format == 1:
		freqs = inputspectrum[:,0]
		fd = inputspectrum[:,1]
		fd_err = inputspectrum[:,2]
	elif format == 2:
		freqs = inputspectrum[:,0] / 1e9
		fd = inputspectrum[:,1]
		fd_err = inputspectrum[:,2]
	elif format == 3:
		wavelength = inputspectrum[:,0]
		freqs = const['c'] / (wavelength*1e3)
		fd = inputspectrum[:,1] / ((1e-26 * (const['c']*1e6) / (wavelength**2))) # Convert from W/m^2/um to W/m^2/Hz then to Jy
		fd_err = inputspectrum[:,2] / ((1e-26 * (const['c']*1e6) / (wavelength**2)))


	num_datapoints = len(freqs)
	print('Fitting ' + str(num_datapoints) + ' data points')

	minfreq = min(freqs)
	maxfreq = max(freqs)
	minflux = min(fd)
	maxflux = max(fd)
	# Assume 10% uncertainties across the board for now
	fd_err = np.sqrt((0.1*fd)**2+fd_err**2)

	# Define which values we want to use in the fit
	goodvals = np.ones(num_datapoints)
	if maxfitfreq != 0:
		goodvals[freqs > maxfitfreq] = 0
	if minfitfreq != 0:
		goodvals[freqs < minfitfreq] = 0
	num_goodvals = np.sum(goodvals)
	badvals = True
	if (num_goodvals == num_datapoints):
		badvals = False

	###
	# Use MPFit
	###

	# Set up the initial parameters
	num_params = 10
	parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
	parinfo=[]
	for i in range(0,num_params):
		parinfo.append(copy.deepcopy(parbase))

	# Starting parameters
	p0 = [10.0,8000.,1.0,1.7,20.,1.,0.,10.0,-1.0, 0]#, 1.0,1.7,1000.]  # standard starting values
	if startparams != 0:
		p0 = startparams

	# Free-free, EM then Te
	parinfo[0]['limited'] = [1,1]
	parinfo[0]['limits'] = [0,1e10]
	parinfo[1]['limited'] = [1,1]
	parinfo[1]['limits'] = [1000.,30000.]
	parinfo[1]['fixed'] = 1  # fix Te
	if nofreefree == True:
		p0[0] = 0
		parinfo[0]['fixed'] = 1
		# parinfo[0]['limits'] = [p0[0],p0[0]]
		parinfo[1]['fixed'] = 1
		# parinfo[1]['limits'] = [p0[1],p0[1]]

	# Thermal dust: optical depth, then index, then temperature
	parinfo[2]['limited'] = [1,1]
	parinfo[2]['limits'] = [0.0,100.]
	parinfo[3]['limited'] = [1,1]
	parinfo[3]['limits'] = [0.5,3.5]
	parinfo[3]['fixed'] = 0  # fix dust emissivity index
	parinfo[4]['limited'] = [1,1]
	parinfo[4]['limits'] = [5.,600.]
	if nodust == True:
		p0[2] = 0
		parinfo[2]['fixed'] = 1
		parinfo[3]['fixed'] = 1
		parinfo[4]['fixed'] = 1

	# Spinning dust amplitude
	parinfo[5]['limited'] = [1,1]
	parinfo[5]['limits'] = [0.,1e25]
	parinfo[9]['fixed'] = 1  # spdust shift
	parinfo[9]['limited'] = [1,1]
	parinfo[9]['limits'] = [0.0,1e10]
	if noame == True:
		p0[5] = 0
		parinfo[5]['fixed'] = 1
		parinfo[9]['fixed'] = 1

	# CMB delta T (CMB) in K
	parinfo[6]['limited'] = [1,1]
	parinfo[6]['limits'] = [-150.0,150.0]#*1e-6
	if nocmb == True:
		p0[6] = 0
		parinfo[6]['fixed'] = 1

	# Synchrotron, amplitude then index
	parinfo[7]['limited'] = [1,1]
	parinfo[7]['limits'] = [0.0,1e10]
	parinfo[8]['limited'] = [1,1]
	parinfo[8]['limits'] = [-3.0,1.0]
	if nosync == True:
		p0[7] = 0
		parinfo[7]['fixed'] = 1
		parinfo[8]['fixed'] = 1
	p0[8] = -1.0
	# parinfo[8]['fixed'] = 1

	# Second thermal dust component: optical depth, then index, then temperature. Disabled by default.
	# parinfo[10]['limited'] = [1,1]
	# parinfo[10]['limits'] = [0.0,100.]
	# parinfo[11]['limited'] = [1,1]
	# parinfo[11]['limits'] = [0.0,3.5]
	# parinfo[12]['limited'] = [1,1]
	# parinfo[12]['limits'] = [5.,10000.] 
	# if nodust2 == True:
	# 	p0[10] = 0
	# 	parinfo[10]['fixed'] = 1
	# 	parinfo[11]['fixed'] = 1
	# 	parinfo[12]['fixed'] = 1
	# if fixdust2temp != 0:
	# 	p0[12] = fixdust2temp
	# 	parinfo[12]['fixed'] = 1

	print('Full info:')
	print(parinfo)

	# Do the fit
	fa = {'x':freqs[goodvals == 1], 'y':fd[goodvals == 1], 'err':fd_err[goodvals == 1]}
	m = mpfit(mpfitfunction, p0, parinfo=parinfo,functkw=fa,xtol=1e-30,quiet=True)

	print('status = ', m.status)
	if (m.status <= 0): 
		print('error message = ', m.errmsg)
	
	print('Number of iterations: ', m.niter)
	dof = len(freqs[goodvals == 1]) - len(m.params)
	print('Degrees of freedom: ', dof)
	print('Chisq: ', m.fnorm / dof)
	print('Parameters: ', m.params)


	plot_results(outdir+srcname+'.pdf',srcname,minfreq,maxfreq,m.params,const,solid_angle,spd_freq,spd_amp,freqs,fd,fd_err,goodvals,minflux,maxflux)

	write_results(outdir+srcname+".txt",srcname,freqs,fd,fd_err,m.params,m.perror)

	###
	# MCMC fitting
	###
	if mcmc == True:

		nll = lambda *args: -lnlike(*args)
		p = p0
		# p = m.params
		# x = freqs[goodvals == 1]
		# y = fd[goodvals == 1]
		# yerr = fd_err[goodvals == 1]
		bounds = np.zeros((num_params,2))
		for i in range(0,num_params):
			if parinfo[i]['fixed'] == 1:
				if p0[i] == 0:
					bounds[i][0] = -1e-5
				else:
					bounds[i][0] = p0[i]*0.999
			elif parinfo[i]['limited'][0] == 1:
				bounds[i][0] = parinfo[i]['limits'][0]
			else:
				bounds[i][0] = None
			if parinfo[i]['fixed'] == 1:
				if p0[i] == 0:
					bounds[i][1] = 1e-5
				else:
					bounds[i][1] = p0[i]*1.0001
			elif parinfo[i]['limited'][1] == 1:
				bounds[i][1] = parinfo[i]['limits'][1]
			else:
				bounds[i][1] = None
		print(bounds)
		result = op.minimize(nll, p, args=(freqs[goodvals == 1], fd[goodvals == 1], fd_err[goodvals == 1]), bounds=bounds, method='L-BFGS-B')#, options={'gtol': 1e-30, 'disp': False, 'maxiter': 4000})
		maxlikelihood = result["x"]
		print("Done")
		print(m.params)
		print(maxlikelihood)
		print(result['success'])
		print(result['message'])

		plot_results(outdir+srcname+'_likelihood.pdf',srcname,minfreq,maxfreq,maxlikelihood,const,solid_angle,spd_freq,spd_amp,freqs,fd,fd_err,goodvals,minflux,maxflux)
		write_results(outdir+srcname+"_likelihood.txt",srcname,freqs,fd,fd_err,maxlikelihood,np.zeros(len(maxlikelihood)))

		# Let's do some MCMC fitting!
		ndim = num_params
		nwalkers = 100
		# not_fixed *= 1e-4 # Use a random distribution only for values that aren't fixed.
		pos = [result["x"] + 1e-4*result["x"]*np.random.randn(ndim) for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=8,args=(bounds, freqs[goodvals == 1], fd[goodvals == 1], fd_err[goodvals == 1]))
		sampler.run_mcmc(pos, 5000)

		samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
		plt.plot(samples[:,0])
		plt.savefig(outdir+srcname+'_test.png')
		plt.close()
		plt.plot(samples[:,1])
		plt.savefig(outdir+srcname+'_test1.png')
		plt.close()
		plt.plot(samples[:,2])
		plt.savefig(outdir+srcname+'_test2.png')
		plt.close()
		plt.plot(samples[:,3])
		plt.savefig(outdir+srcname+'_test3.png')
		plt.close()


		# samples[:, 2] = np.exp(samples[:, 2])
		p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
		vals = np.zeros(num_params)
		vals[0] = p1[0]
		vals[1] = p2[0]
		vals[2] = p3[0]
		vals[3] = p4[0]
		vals[4] = p5[0]
		vals[5] = p6[0]
		vals[6] = p7[0]
		vals[7] = p8[0]
		vals[8] = p9[0]
		vals[9] = p10[0]
		# vals[10] = p11[0]
		# vals[11] = p12[0]
		# vals[12] = p13[0]
		errs = np.zeros(num_params)
		errs[0] = p1[1]
		errs[1] = p2[1]
		errs[2] = p3[1]
		errs[3] = p4[1]
		errs[4] = p5[1]
		errs[5] = p6[1]
		errs[6] = p7[1]
		errs[7] = p8[1]
		errs[8] = p9[1]
		errs[9] = p10[1]
		# errs[10] = p11[1]
		# errs[11] = p12[1]
		# errs[12] = p13[1]
		print(errs)
		# print(m_mcmc)
		# print(b_mcmc)
		# print(f_mcmc)
		# print(vals)
		# print(np.shape(vals))
		# exit()
		plot_range = []
		for i in range(0,num_params):
			if parinfo[i]['fixed'] == 1:
				if vals[i] == 0.0:
					plot_range.append((-0.01, 0.01))
				else:
					plot_range.append((0.99*vals[i], 1.01*vals[i]))
			elif vals[i] == 0.0:
				plot_range.append((-0.01, 0.01))
			else:
				plot_range.append((np.min(samples[:,i]),np.max(samples[:,i])))
				# plot_range.append((vals[i]-3.0*errs[i],vals[i]+3.0*errs[i]))
				# plot_range[i][1] = max(samples[:,i])
		print(plot_range)
		fig = corner.corner(samples, labels=['EM', 'Te', 'Damp', 'Dalph', 'Dtemp', 'SDamp', 'dTCMB','Samp','Salph','SDshift'],range=plot_range,show_titles=True)#,'Damp', 'Dalph', 'Dtemp'
                      # truths=[m_true, b_true, np.log(f_true)])
		fig.savefig(outdir+srcname+"_triangle.png")
		plt.close()

		plot_results(outdir+srcname+'_mcmc.pdf',srcname,minfreq,maxfreq,vals,const,solid_angle,spd_freq,spd_amp,freqs,fd,fd_err,goodvals,minflux,maxflux)
		write_results(outdir+srcname+"_mcmc.txt",srcname,freqs,fd,fd_err,vals,errs)


# That's all, folks!
