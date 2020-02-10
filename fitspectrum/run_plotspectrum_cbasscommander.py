from plotspectrum import plotspectrum
from smoothmap import smoothmap
import numpy as np
import healpy as hp
from astroutils import *

# Plot the spectra for the C-BASS comparison Commander runs of January 2020
# Overall parameters
basedir = '/Users/mpeel/Desktop/Low-Frequency_Survey2/'

run_rmsspectra = True
run_aperflux = False
run_extramapplots = True

sources = [['Perseus', 160.26,-18.62], ['California',160.26,-12.5], ['M31', 121.174322, -21.573311], ['M42', 209.0137, -19.3816], ['TauA', 184.55745, -05.78436]]

# 0  = 'GAL020  '           / 20% sky coverage                               
# 1  = 'GAL040  '           / 40% sky coverage                               
# 2  = 'GAL060  '           / 60% sky coverage                               
# 3  = 'GAL070  '           / 70% sky coverage                               
# 4  = 'GAL080  '           / 80% sky coverage                               
# 5  = 'GAL090  '           / 90% sky coverage                               
# 6  = 'GAL097  '           / 97% sky coverage                               
# 7  = 'GAL099  '           / 99% sky coverage   
mask_min=[['/Users/mpeel/Documents/maps/commander2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits',2048,4]]
mask_max=[['/Users/mpeel/Documents/maps/commander2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits',2048,5]]
nside = 64
res = 60.0
aper_inner_radius = res
aper_outer_radius1 = res
aper_outer_radius2 = 1.5*res
units = 'JyPix'
noise_model = 1

freqbands_I = [[0.408,0.408,'Haslam','grey', 0.5,'-'],
[0.820,0.820,'Dwingeloo','grey', 0.5,'-.'],
[1.42,1.420,'Stockert','grey', 0.5,'--'],
[2.3,2.30,'HartRAO','grey', 0.5,':'],
[4.5,5.5,'C-BASS','g', 0.3],
[19.5,25.0,'WMAP','b', 0.1],
[28.0,37.0,'','b', 0.1],
[35,46.0,'','b', 0.1],
[53.0,69.0,'','b', 0.1],
[82.0,106.0,'','b', 0.1],
[30.0-(0.2*30.0)/2,30.0+(0.2*30.0)/2,'Planck','r', 0.1], # Planck
[44.0-(0.2*44.0)/2,44.0+(0.2*44.0)/2,'','r', 0.1], # Planck
[70.0-(0.2*70.0)/2,70.0+(0.2*70.0)/2,'','r', 0.1], # Planck
[100.0-(0.33*100.0)/2,100.0+(0.33*100.0)/2,'','r', 0.1], # Planck
[143.0-(0.33*143.0)/2,143.0+(0.33*143.0)/2,'','r', 0.1], # Planck
[217.0-(0.33*217.0)/2,217.0+(0.33*217.0)/2,'','r', 0.1], # Planck
[353.0-(0.33*353.0)/2,353.0+(0.33*353.0)/2,'','r', 0.1], # Planck
[545.0-(0.33*545.0)/2,545.0+(0.33*545.0)/2,'','r', 0.1], # Planck
[857.0-(0.33*857.0)/2,857.0+(0.33*857.0)/2,'','r', 0.1]] # Planck

# vals: 0=sync_amp, 1=ff_amp, 2=ff_EM, 3=CMB_amp, 4=AME1_amp, 5=AME1_freq, 6=AME2_amp, 7=td_amp, 8=td, 9=td, 10=sync_beta
# params, 0=filename, 1=nside, 2=res, 3=hdu, 4=colnum, 5=valsnum, 6=rescale

# C-BASS-only

runs = ['haslam', 'cbass','chipass', 'spass','cbass_chipass', 'cbass_spass', 'haslam_cbass', 'haslam_chipass', 'haslam_spass','spass_chipass']

for run in runs:
	if '_' not in run:
		subdir = 'one/'+run+'/'
	else:
		subdir = 'two/'+run+'/'

	maps_I = [[basedir+subdir+'synch_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,0,1.0,'A_Sync'],
	[basedir+subdir+'synch_beta_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,10,1.0,'Beta_Sync'],
	[basedir+subdir+'ff_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,1,1.0,'A_ff'],
	[basedir+subdir+'ff_T_e_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,2,1.0,'Te'],
	[basedir+subdir+'cmb_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,3,1.0,'CMB'],
	[basedir+subdir+'ame1_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,4,1.0,'A_AME1'],
	[basedir+subdir+'ame1_nup_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,5,1.0,'nu_AME1']]

	if '_' not in run:
		null = 0
	else:
		for i in range(len(maps_I)):
			maps_I[i][0] = maps_I[i][0].replace('base_','')

	othermask = [basedir+'masks/mask_'+run+'_n0064.fits',64,0]
	print(mask_min)
	mask_min2 = mask_min.copy()
	mask_max2 = mask_max.copy()

	# if run == 'haslam':
	# 	othermask = []
	# else:
	mask_min2.append(othermask)
	mask_max2.append(othermask)
	print(maps_I)
	print(mask_min2)

	if run_rmsspectra:
		try:
			plotspectrum(outdir=basedir+subdir, name='spectra_'+run, maps=maps_I, mask_min=mask_min2,mask_max=mask_max2,spd_file='amemodels/spdust2_cnm.dat',minfreq=0.2,maxfreq=2000.0,numxpoints=1000,ymin=1,ymax=5e7,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True,freqbands=freqbands_I,nside=nside,res=res,pol=False,othermask=othermask)
		except:
			continue

	if run_extramapplots:
		mask = hp.read_map(othermask[0])
		mask2 = hp.read_map(mask_min[0][0])
		mask2 = hp.ud_grade(mask2, nside)
		mask[mask2 == 0] = 0
		for line in maps_I:
			mapdata = hp.read_map(line[0])
			numsigma = 2.0
			if np.min(mapdata[mask==1]) == 0:
				minplot = 0
			else:
				minplot = -numsigma*np.std(mapdata[mask==1])
			hp.mollview(mapdata-np.mean(mapdata[mask==1]),min=minplot,max=numsigma*np.std(mapdata[mask==1]),coord=['G','C'])
			plt.savefig(line[0].replace('.fits','_dec.pdf'))


	if run_aperflux:
		logfile = open(basedir+subdir+"/_aperflux.txt","w")
		logfile.write('Source')
		for line in maps_I:
			logfile.write(" & " + str(line[7].replace('_','')))
		logfile.write('\\\\\n')
		for source in sources:
			try:
				logfile.write('\n'+str(source[0]))
				# logfile.write('\n'+str(source)+'\n\\begin{itemize}\n')
				for line in maps_I:
					mapdata = hp.read_map(line[0])
					vals = haperflux(mapdata, 0.0, res, source[1], source[2], aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=noise_model, silent=True,quickplot=line[0].replace('.fits','_'+str(source[0])+'_aper.png'))
	#				logfile.write('\\item '+str(line[7].replace('_','-')) + "  {0:.2f}".format(vals[0]) + ' +- ' + "{0:.2f}".format(vals[1])+'\n')
					logfile.write(" & ${0:.2f}".format(vals[0]) + ' \\pm ' + "{0:.2f}$".format(vals[1]))
				# logfile.write("\\end{itemize}")
				logfile.write("\\\\\n")
			except:
				continue
		logfile.close()


