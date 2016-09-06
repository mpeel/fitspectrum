#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Plot the spectrum of the Commander maps
# 
# Version history:
#
# 15-Aug-2016  M. Peel       Started, for C-BASS project paper

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


def plotspectrum(outdir='', name='plot', maps=[''],mask_min=[''],mask_max=[''],spd_file='amemodels/spdust2_wim.dat', minfreq=0, maxfreq=0, numxpoints=1000, ymin=0, ymax=0, nosync=False, nofreefree=False, noame=False, nocmb=False, nodust=False, nodust2=True, nside=256,res=60.0,freqbands=[[0,0,'nofreq','k']],pol=False,legend=True):

	ensure_dir(outdir)

	# Read in the spinning dust curve
	spd_freq, spd_amp = np.loadtxt(spd_file, dtype=float,usecols=(0,1),comments=';',unpack=True)
	galprop_freq, galprop_amp = np.loadtxt('commander2015/Synchrotron_template_GHz_extended.txt',dtype=float,usecols=(0,1),comments=";",unpack=True)

	# Figure out the X range
	x = np.arange(minfreq,maxfreq,(maxfreq-minfreq)/float(numxpoints))

	# Read in the mask for the minimum values, and ud_grade it if needed
	nummasks_min = len(mask_min)
	maskmap_min = np.ones(hp.nside2npix(nside))
	for i in range(0,nummasks_min):
		maskmap_min_temp = hp.read_map(mask_min[i][0], field=mask_min[i][2])
		if mask_min[i][1] != nside:
			maskmap_min_temp = hp.ud_grade(maskmap_min_temp, nside)
		maskmap_min *= maskmap_min_temp
	print 'Minimum mask: ' + str(100.0*(sum(maskmap_min[maskmap_min == 1]) / len(maskmap_min)))


	# Read in the mask for the maximum values, and ud_grade it if needed
	nummasks_max = len(mask_max)
	maskmap_max = np.ones(hp.nside2npix(nside))
	for i in range(0,nummasks_max):
		maskmap_max_temp = hp.read_map(mask_max[i][0], field=mask_max[i][2])
		if mask_max[i][1] != nside:
			maskmap_max_temp = hp.ud_grade(maskmap_max_temp, nside)
		maskmap_max *= maskmap_max_temp
	print 'Maximum mask: ' + str(100.0*(sum(maskmap_max[maskmap_max == 1]) / len(maskmap_max)))

	# Read in the maps, and smooth them if needed
	nummaps = len(maps)
	minvals = np.zeros(nummaps)
	maxvals = np.zeros(nummaps)
	for i in range(0,nummaps):
		if maps[i][2] != res:
			print "Hello!"
			newfilename = maps[i][0][:-5]+'_'+str(nside)+"_"+str(res)+".fits"
			smoothmap(maps[i][0],newfilename, np.sqrt((float(res)/60.0)**2-(float(maps[i][2]))**2),pol=pol,nside_out=nside)
			maps[i][0] = newfilename

		print maps[i][0]
		print maps[i][3]
		print maps[i][4]
		mapdata = hp.read_map(maps[i][0], maps[i][4],hdu=maps[i][3])
		# if pol:
		# 	mapdata = np.sqrt(mapdata[0]**2+mapdata[1]**2)
		if maps[i][1] != nside:
			mapdata = hp.ud_grade(mapdata, nside)


		minvals[maps[i][5]] = np.sqrt(np.mean(np.square(mapdata[maskmap_min == 1]))) / maps[i][6]
		maxvals[maps[i][5]] = np.sqrt(np.mean(np.square(mapdata[maskmap_max == 1]))) / maps[i][6]

	print minvals
	print maxvals

	# # Generate the model and plot it
	# if nodust == False:
	# 	model_dust1 = thermaldust(const, x, m.params[2], m.params[3], m.params[4], const['dust_optical_depth_freq'], solid_angle)
	# 	
	# if nodust2 == False:
	# 	model_dust2 = thermaldust(const, x, m.params[10], m.params[11], m.params[12], const['dust_optical_depth_freq'], solid_angle)
	# 	plt.plot(x, model_dust2, 'g')#, freqs, fd, 'g')

	# model_overall = spectrum(m.params, x=x)
	# plt.plot(x, model_overall, 'black')
	# # Add the data to the plot
	# plt.errorbar(freqs[goodvals == 1], fd[goodvals == 1], fd_err[goodvals == 1])
	# if badvals:
	# 	plt.errorbar(freqs[goodvals == 0], fd[goodvals == 0], fd_err[goodvals == 0])

	# Formatting, and output
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Frequency (GHz)')
	plt.ylabel('Brightness temperature r.m.s. (uK)')
	plt.xlim(xmin=minfreq,xmax=maxfreq)
	plt.ylim(ymin=ymin,ymax=ymax)

	if pol == False:
		sync_spectrum_min = syncshifted_comm(x, minvals[0], 4e-3, galprop_freq, galprop_amp)
		sync_spectrum_max = syncshifted_comm(x, maxvals[0], 4e-3, galprop_freq, galprop_amp)
		plt.fill_between(x, sync_spectrum_min, sync_spectrum_max,facecolor='magenta',lw=0,zorder=10,label='Synchrotron')

		freefree_spectrum_min = freefree(const, x, minvals[1], minvals[2], solid_angle, equation=1,comm=1)
		freefree_spectrum_max = freefree(const, x, maxvals[1], maxvals[2], solid_angle, equation=1,comm=1)
		plt.fill_between(x, freefree_spectrum_min, freefree_spectrum_max,facecolor='blue',lw=0,zorder=9,label='Free-free')

		ame_spectrum_min = spinningdust_comm(x, minvals[4], minvals[5], spd_amp, spd_freq, 1) + spinningdust_comm(x, minvals[6], 33.35, spd_amp, spd_freq, 2)
		ame_spectrum_max = spinningdust_comm(x, maxvals[4], maxvals[5], spd_amp, spd_freq, 1) + spinningdust_comm(x, maxvals[6], 33.35, spd_amp, spd_freq, 2)
		plt.fill_between(x, ame_spectrum_min, ame_spectrum_max,facecolor='gold',lw=0,zorder=9,label="AME")

		cmb_spectrum_min = cmb_comm(const, x, minvals[3])
		cmb_spectrum_max = cmb_comm(const, x, maxvals[3])
		plt.fill_between(x, cmb_spectrum_min, cmb_spectrum_max,facecolor='k',lw=1.0,zorder=11,label='CMB')

		thermaldust_spectrum_min = thermaldust_comm(const, x, minvals[7], minvals[9], minvals[8])
		thermaldust_spectrum_max = thermaldust_comm(const, x, maxvals[7], maxvals[9], maxvals[8])
		plt.fill_between(x, thermaldust_spectrum_min, thermaldust_spectrum_max,facecolor='r',lw=0,zorder=9,label='Thermal dust')

		fd_min = sync_spectrum_min + freefree_spectrum_min + ame_spectrum_min + thermaldust_spectrum_min
		fd_max = sync_spectrum_max + freefree_spectrum_max + ame_spectrum_max + thermaldust_spectrum_max
		plt.plot(x, fd_min, 'k--',zorder=12,label='Total-CMB')
		plt.plot(x, fd_max, 'k--',zorder=12)
	else:
		sync_spectrum_min = syncshifted_pol_comm(x, np.sqrt(minvals[0]**2+minvals[1]**2), 4e-3, galprop_freq, galprop_amp)
		sync_spectrum_max = syncshifted_pol_comm(x, np.sqrt(maxvals[0]**2+maxvals[1]**2), 4e-3, galprop_freq, galprop_amp)
		plt.fill_between(x, sync_spectrum_min, sync_spectrum_max,facecolor='magenta',lw=0,zorder=10,label="Synchrotron")
		
		cmb_spectrum_min = cmb_comm(const, x, np.sqrt(minvals[2]**2+minvals[3]**2))
		cmb_spectrum_max = cmb_comm(const, x, np.sqrt(maxvals[2]**2+maxvals[3]**2))
		plt.fill_between(x, cmb_spectrum_min, cmb_spectrum_max,facecolor='k',lw=1.0,zorder=11,label='CMB')

		thermaldust_spectrum_min = thermaldust_comm(const, x, np.sqrt(minvals[4]**2+minvals[5]**2), minvals[7], minvals[6])
		thermaldust_spectrum_max = thermaldust_comm(const, x, np.sqrt(maxvals[4]**2+maxvals[5]**2), maxvals[7], maxvals[6])
		plt.fill_between(x, thermaldust_spectrum_min, thermaldust_spectrum_max,facecolor='r',lw=0,zorder=9,label='Thermal dust')

		fd_min = sync_spectrum_min + thermaldust_spectrum_min
		fd_max = sync_spectrum_max + thermaldust_spectrum_max
		plt.plot(x, fd_min, 'k--',zorder=12,label='Total-CMB')
		plt.plot(x, fd_max, 'k--',zorder=12)

	numfreqbands = len(freqbands)
	for i in range (0,numfreqbands):
		if (freqbands[i][0] != freqbands[i][1]):
			plt.axvspan(freqbands[i][0], freqbands[i][1], color=freqbands[i][3], alpha=0.1, lw=1.0,zorder=0,label=freqbands[i][2])
		else:
			plt.axvspan(freqbands[i][0], freqbands[i][1], color=freqbands[i][3], alpha=0.3, lw=1.0,zorder=0,label=freqbands[i][2])
		# plt.plot((freqbands[i][0], freqbands[i][1]), (ymin, ymax), freqbands[i][3])

	if legend == True:
		l = plt.legend(prop={'size':8})
		l.set_zorder(20)
	plt.savefig(outdir+name+'.pdf')
	plt.close()

	# outputfile = open(outdir+srcname+".dat", "w")
	# outputfile.write("# " + srcname + "\n")
	# for i in range(0,num_datapoints):
	# 	outputfile.write(str(freqs[i]) + "	" + str(fd[i]) + "	" + str(fd_err[i]) + "\n")

	# if nodust == False:
	# 	outputfile.write('# Dust amplitude: ' + str(m.params[2]) + " +- " + str(m.perror[2]) + '\n')
	# 	outputfile.write('# Dust index: ' + str(m.params[3]) + " +- " + str(m.perror[3]) + '\n')
	# 	outputfile.write('# Dust temperature: ' + str(m.params[4]) + " +- " + str(m.perror[4]) + '\n')
	# if nodust2 == False:
	# 	outputfile.write('# Dust2 amplitude: ' + str(m.params[10]) + " +- " + str(m.perror[10]) + '\n')
	# 	outputfile.write('# Dust2 index: ' + str(m.params[11]) + " +- " + str(m.perror[11]) + '\n')
	# 	outputfile.write('# Dust2 temperature: ' + str(m.params[12]) + " +- " + str(m.perror[12]) + '\n')



	# outputfile.close()


# That's all, folks!
