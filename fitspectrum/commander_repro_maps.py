#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Plot the spectrum of the Commander maps
# 
# Version history:
#
# 26-Jan-2018  M. Peel       Started, forked from plotspectrum.py

import numpy as np
import scipy.optimize as op
import healpy as hp
from spectra import *
from astroutils import *
import copy
import matplotlib.pyplot as plt
from smoothmap import smoothmap

# Define some constants, used later in the SED functions
const = get_spectrum_constants()
solid_angle = 1.0e-10 # For now


def commander_repro_maps(outdir='', name='plot', maps=[''],spd_file='amemodels/spdust2_wim.dat', nside=256,res=60.0,freqbands=[[0,0,'nofreq','k']],pol=False,legend=True,freq=30.0,use_ame1=True,use_ame2=True,use_ff=True,use_sync=True,use_cmb=True,use_thermaldust=True,syncmodel=1,mask=[]):

	# ensure_dir(outdir)

	# Read in the spinning dust curve
	spd_freq, spd_amp = np.loadtxt(spd_file, dtype=float,usecols=(0,1),comments=';',unpack=True)
	if syncmodel==1:
		galprop_freq, galprop_amp = np.loadtxt('commander2015/Synchrotron_template_GHz_extended.txt',dtype=float,usecols=(0,1),comments=";",unpack=True)

	# Read in the maps, and smooth them if needed
	# nummaps = len(maps)
	nummaps = 11
	npix = hp.nside2npix(nside)
	mapdata = np.zeros((nummaps, npix))

	for i in range(0,len(maps)):
		if maps[i][2] != res:
			print("Hello!")
			newfilename = maps[i][0][:-5]+'_'+str(nside)+"_"+str(res)+".fits"
			smoothmap(maps[i][0],newfilename, np.sqrt((float(res)/60.0)**2-(float(maps[i][2]))**2),pol=pol,nside_out=nside)
			maps[i][0] = newfilename

		print(maps[i][0])
		print(maps[i][3])
		print(maps[i][4])
		mapdata_temp = hp.read_map(maps[i][0], maps[i][4],hdu=maps[i][3])
		if maps[i][1] != nside:
			mapdata_temp = hp.ud_grade(mapdata_temp, nside)
		mapdata[maps[i][5]] = mapdata_temp

	mapmask = np.ones(hp.nside2npix(nside))
	if len(mask) > 0:
		mapmask = hp.read_map(mask[0], field=mask[2])
		if mask[1] != nside:
			print('udgrading to ' + str(nside))
			mapmask = hp.ud_grade(mapmask, nside)
		hp.mollview(mapmask)
		plt.savefig(outdir+name+'_mask.pdf')
		plt.clf()


	sync = np.zeros(npix)
	ff = np.zeros(npix)
	ame = np.zeros(npix)
	ame1 = np.zeros(npix)
	ame2 = np.zeros(npix)
	cmb = np.zeros(npix)
	thermaldust = np.zeros(npix)
	total_fg = np.zeros(npix)
	for i in range(0,npix-1):
		if mapmask[i] == 1:
			if use_sync:
				if syncmodel == 1:
					sync[i] = syncshifted_comm(freq, mapdata[0][i], 4e-3, galprop_freq, galprop_amp)
				else:
					sync[i] = synchrotron(const, freq, 22.8, mapdata[0][i], mapdata[10][i])

			if use_ff:
				ff[i] = freefree(const, freq, mapdata[1][i], mapdata[2][i], solid_angle, equation=1,comm=1)
			if use_ame1:
				ame1[i] = spinningdust_comm(freq, mapdata[4][i], mapdata[5][i], spd_amp, spd_freq, 1)
			if use_ame2:
				ame2[i] = spinningdust_comm(freq, mapdata[6][i], 33.35, spd_amp, spd_freq, 2)
			ame[i] =  ame1[i] + ame2[i]
			if use_cmb:
				cmb[i] = cmb_comm(const, freq, mapdata[3][i])
			if use_thermaldust:
				thermaldust[i] = thermaldust_comm(const, freq, mapdata[7][i], mapdata[9][i], mapdata[8][i])
			total_fg[i] = sync[i] + ff[i] + ame[i] + thermaldust[i]

	numstd = 3.0
	if use_sync:
		hp.write_map(outdir+name+"_commander_sync_"+str(freq)+".fits", sync*1e-3,overwrite=True)
		median = np.median(sync*1e-3)
		std = np.std(sync*1e-3-median)
		if np.min(sync*1e-3) > median-numstd*std:
			minval = np.min(sync*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(sync*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_sync_"+str(freq)+".png")
	if use_ff:
		hp.write_map(outdir+name+"_commander_freefree_"+str(freq)+".fits", ff*1e-3,overwrite=True)
		median = np.median(ff*1e-3)
		std = np.std(ff*1e-3-median)
		if np.min(ff*1e-3) > median-numstd*std:
			minval = np.min(ff*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(ff*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_freefree_"+str(freq)+".png")
	if use_ame1 and use_ame2:
		hp.write_map(outdir+name+"_commander_ame_"+str(freq)+".fits", ame*1e-3,overwrite=True)
		median = np.median(ame*1e-3)
		std = np.std(ame*1e-3-median)
		if np.min(ame*1e-3) > median-numstd*std:
			minval = np.min(ame*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(ame*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_ame_"+str(freq)+".png")
	if use_ame1:
		hp.write_map(outdir+name+"_commander_ame1_"+str(freq)+".fits", ame1*1e-3,overwrite=True)
		median = np.median(ame1*1e-3)
		std = np.std(ame1*1e-3-median)
		if np.min(ame1*1e-3) > median-numstd*std:
			minval = np.min(ame1*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(ame1*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_ame1_"+str(freq)+".png")
	if use_ame2:
		hp.write_map(outdir+name+"_commander_ame2_"+str(freq)+".fits", ame2*1e-3,overwrite=True)
		median = np.median(ame2*1e-3)
		std = np.std(ame2*1e-3-median)
		if np.min(ame2*1e-3) > median-numstd*std:
			minval = np.min(ame2*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(ame2*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_ame2_"+str(freq)+".png")
	if use_cmb:
		hp.write_map(outdir+name+"_commander_cmb_"+str(freq)+".fits", cmb*1e-3,overwrite=True)
		median = np.median(cmb*1e-3)
		std = np.std(cmb*1e-3-median)
		if np.min(cmb*1e-3) > median-numstd*std:
			minval = np.min(cmb*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(cmb*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_cmb_"+str(freq)+".png")
	if use_thermaldust:
		hp.write_map(outdir+name+"_commander_thermaldust_"+str(freq)+".fits", thermaldust*1e-3,overwrite=True)
		median = np.median(thermaldust*1e-3)
		std = np.std(thermaldust*1e-3-median)
		if np.min(thermaldust*1e-3) > median-numstd*std:
			minval = np.min(thermaldust*1e-3)
		else:
			minval = median-numstd*std
		hp.mollview(thermaldust*1e-3, min=minval, max=median+numstd*std)
		plt.savefig(outdir+name+"_commander_thermaldust_"+str(freq)+".png")

	hp.write_map(outdir+name+"_commander_total_fg_"+str(freq)+".fits", total_fg*1e-3,overwrite=True)
	median = np.median(total_fg*1e-3)
	std = np.std(total_fg*1e-3-median)
	if np.min(total_fg*1e-3) > median-numstd*std:
		minval = np.min(total_fg*1e-3)
	else:
		minval = median-numstd*std
	hp.mollview(total_fg*1e-3, min=minval, max=median+numstd*std)
	plt.savefig(outdir+name+"_commander_total_fg_"+str(freq)+".png")


# That's all, folks!
