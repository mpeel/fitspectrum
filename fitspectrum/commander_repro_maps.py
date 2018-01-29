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


def commander_repro_maps(outdir='', name='plot', maps=[''],spd_file='amemodels/spdust2_wim.dat', nside=256,res=60.0,freqbands=[[0,0,'nofreq','k']],pol=False,legend=True,freq=30.0):

	# ensure_dir(outdir)

	# Read in the spinning dust curve
	spd_freq, spd_amp = np.loadtxt(spd_file, dtype=float,usecols=(0,1),comments=';',unpack=True)
	galprop_freq, galprop_amp = np.loadtxt('commander2015/Synchrotron_template_GHz_extended.txt',dtype=float,usecols=(0,1),comments=";",unpack=True)

	# Read in the maps, and smooth them if needed
	nummaps = len(maps)
	npix = hp.nside2npix(nside)
	mapdata = np.zeros((nummaps, npix))

	for i in range(0,nummaps):
		if maps[i][2] != res:
			print "Hello!"
			newfilename = maps[i][0][:-5]+'_'+str(nside)+"_"+str(res)+".fits"
			smoothmap(maps[i][0],newfilename, np.sqrt((float(res)/60.0)**2-(float(maps[i][2]))**2),pol=pol,nside_out=nside)
			maps[i][0] = newfilename

		print maps[i][0]
		print maps[i][3]
		print maps[i][4]
		mapdata_temp = hp.read_map(maps[i][0], maps[i][4],hdu=maps[i][3])
		if maps[i][1] != nside:
			mapdata_temp = hp.ud_grade(mapdata_temp, nside)
		mapdata[i] = mapdata_temp

	sync = np.zeros(npix)
	ff = np.zeros(npix)
	ame = np.zeros(npix)
	cmb = np.zeros(npix)
	thermaldust = np.zeros(npix)
	total_fg = np.zeros(npix)
	for i in range(0,npix-1):
		sync[i] = syncshifted_comm(freq, mapdata[0][i], 4e-3, galprop_freq, galprop_amp)
		ff[i] = freefree(const, freq, mapdata[1][i], mapdata[2][i], solid_angle, equation=1,comm=1)
		ame[i] = spinningdust_comm(freq, mapdata[4][i], mapdata[5][i], spd_amp, spd_freq, 1) + spinningdust_comm(freq, mapdata[6][i], 33.35, spd_amp, spd_freq, 2)
		cmb[i] = cmb_comm(const, freq, mapdata[3][i])
		thermaldust[i] = thermaldust_comm(const, freq, mapdata[7][i], mapdata[9][i], mapdata[8][i])
		total_fg[i] = sync[i] + ff[i] + ame[i] + thermaldust[i]

	hp.write_map(outdir+"commander_sync_"+str(freq)+".fits", sync*1e-3)
	hp.write_map(outdir+"commander_freefree_"+str(freq)+".fits", ff*1e-3)
	hp.write_map(outdir+"commander_ame_"+str(freq)+".fits", ame*1e-3)
	hp.write_map(outdir+"commander_cmb_"+str(freq)+".fits", cmb*1e-3)
	hp.write_map(outdir+"commander_thermaldust_"+str(freq)+".fits", thermaldust*1e-3)
	hp.write_map(outdir+"commander_total_fg_"+str(freq)+".fits", total_fg*1e-3)


# That's all, folks!
