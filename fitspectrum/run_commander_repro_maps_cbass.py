from commander_repro_maps import commander_repro_maps
from smoothmap import smoothmap
import numpy as np

basedir = '/Users/mpeel/Desktop/Low-Frequency_Survey2/'

runs = ['haslam', 'cbass','chipass', 'spass', 'cbass_chipass', 'cbass_spass', 'haslam_cbass', 'haslam_chipass', 'haslam_spass','spass_chipass']
nside = 64
res = 60.0

for run in runs:
	if '_' not in run:
		subdir = 'one/'+run+'/'
	else:
		subdir = 'two/'+run+'/'

	# vals: 0=sync_amp, 1=ff_amp, 2=ff_EM, 3=CMB_amp, 4=AME1_amp, 5=AME1_freq, 6=AME2_amp, 7=td_amp, 8=td, 9=td, 10=sync_beta
	# params, 0=filename, 1=nside, 2=res, 3=hdu, 4=colnum, 5=valsnum, 6=rescale
	maps = [[basedir+subdir+'synch_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,0,1.0,'A_Sync'],
	[basedir+subdir+'synch_beta_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,10,1.0,'Beta_Sync'],
	[basedir+subdir+'ff_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,1,1.0,'A_ff'],
	[basedir+subdir+'ff_T_e_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,2,1.0,'Te'],
	[basedir+subdir+'cmb_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,3,1.0,'CMB'],
	[basedir+subdir+'ame1_amplitude_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,4,1.0,'A_AME1'],
	[basedir+subdir+'ame1_nup_base_'+run+'_60arcmin_n0064.fits',nside,res,1,0,5,1.0,'nu_AME1']]

	if '_' not in run:
		null = 0
	else:
		for i in range(len(maps)):
			maps[i][0] = maps[i][0].replace('base_','')

	mask = 	[basedir+'masks/mask_'+run+'_n0064.fits',64,0]

	frequencies = [4.76, 22.8, 28.4]
	# names = ['cbass', 'wmapk', 'planck30']

	for i in range(0,len(frequencies)):
		commander_repro_maps(outdir=basedir+'commander_freq_maps/', name=run, maps=maps, mask=mask,spd_file='amemodels/spdust2_cnm.dat',nside=nside,res=60.0,freq=frequencies[i],use_ame2=False,use_thermaldust=False,syncmodel=2,amemodel=3)

