#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Plot the spectrum of the Commander maps
# 
# Version history:
#
# 15-Aug-2016  M. Peel       Started, for C-BASS project paper
# 28-Jan-2020  M. Peel       Upgrade to handle C-BASS commander maps

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


def plotspectrum(outdir='', name='plot', maps=[''],mask_min=[''],mask_max=[''],spd_file='amemodels/spdust2_wim.dat', galprop_file='', minfreq=0, maxfreq=0, numxpoints=1000, ymin=0, ymax=0, nosync=False, nofreefree=False, noame=False, nocmb=False, nodust=False, nodust2=True, nside=256,res=60.0,freqbands=[[0,0,'nofreq','k']],pol=False,legend=True,othermask=[]):

	# ensure_dir(outdir)
	plt.style.use('classic')
	plt.xticks(fontsize=17)
	plt.yticks(fontsize=17)

	const = get_spectrum_constants()
	# Read in the spinning dust curve
	if spd_file != '':
		spd_freq, spd_amp = np.loadtxt(spd_file, dtype=float,usecols=(0,1),comments=';',unpack=True)
	if galprop_file != '':
		galprop_freq, galprop_amp = np.loadtxt('commander2015/Synchrotron_template_GHz_extended.txt',dtype=float,usecols=(0,1),comments=";",unpack=True)

	# Figure out the X range
	x = np.arange(minfreq,maxfreq,(maxfreq-minfreq)/float(numxpoints))

	# Read in the mask for the minimum values, and ud_grade it if needed
	nummasks_min = len(mask_min)
	maskmap_min = np.ones(hp.nside2npix(nside))
	for i in range(0,nummasks_min):
		maskmap_min_temp = hp.read_map(mask_min[i][0], field=mask_min[i][2])
		if mask_min[i][1] != nside:
			print('udgrading to ' + str(nside))
			maskmap_min_temp = hp.ud_grade(maskmap_min_temp, nside)
		maskmap_min *= maskmap_min_temp
	print('Minimum mask: ' + str(100.0*(sum(maskmap_min[maskmap_min == 1]) / len(maskmap_min))))
	hp.mollview(maskmap_min)
	plt.savefig(outdir+'mask_min.pdf')
	plt.clf()


	# Read in the mask for the maximum values, and ud_grade it if needed
	nummasks_max = len(mask_max)
	maskmap_max = np.ones(hp.nside2npix(nside))
	for i in range(0,nummasks_max):
		maskmap_max_temp = hp.read_map(mask_max[i][0], field=mask_max[i][2])
		if mask_max[i][1] != nside:
			print('udgrading to ' + str(nside))
			maskmap_max_temp = hp.ud_grade(maskmap_max_temp, nside)
		maskmap_max *= maskmap_max_temp
	print('Maximum mask: ' + str(100.0*(sum(maskmap_max[maskmap_max == 1]) / len(maskmap_max))))
	hp.mollview(maskmap_max)
	plt.savefig(outdir+'mask_max.pdf')
	plt.clf()

	maskmap_other = np.ones(hp.nside2npix(nside))
	if len(othermask) > 0:
		num = 1
	else:
		num = 0
	for i in range(0,num):
		maskmap_other = hp.read_map(othermask[0], field=othermask[2])
		if othermask[1] != nside:
			print('udgrading to ' + str(nside))
			maskmap_other = hp.ud_grade(maskmap_other, nside)
		hp.mollview(maskmap_other)
		plt.savefig(outdir+'mask_other.pdf')
		plt.clf()

	# Read in the maps, and smooth them if needed
	nummaps = 11#len(maps)
	minvals = np.zeros(nummaps)
	maxvals = np.zeros(nummaps)
	for i in range(0,len(maps)):

		if maps[i][2] != res:
			print("Hello!")
			newfilename = maps[i][0][:-5]+'_'+str(nside)+"_"+str(res)+".fits"
			smoothmap(maps[i][0],newfilename, np.sqrt((float(res)/60.0)**2-(float(maps[i][2]))**2),pol=pol,nside_out=nside)
			maps[i][0] = newfilename

		print(maps[i][0])
		print(maps[i][3])
		print(maps[i][4])
		mapdata = hp.read_map(maps[i][0], maps[i][4],hdu=maps[i][3])

		hp.mollview(mapdata,norm='hist')
		plt.savefig(maps[i][0].replace('.fits','_hist.pdf'))
		plt.clf()
		if np.min(mapdata[maskmap_min==1]) == 0:
			minplot = 0
		else:
			minplot = -3.0*np.std(mapdata[maskmap_min==1])
		hp.mollview(mapdata,min=minplot,max=3.0*np.std(mapdata[maskmap_min==1]))
		plt.savefig(maps[i][0].replace('.fits','_std.pdf'))
		plt.clf()

		# Create some quick histograms
		histrange_min = np.min(mapdata[maskmap_min==1])
		histrange_max = np.max(mapdata[maskmap_min==1])
		print(histrange_min)
		print(histrange_max)
		if histrange_min != histrange_max:
			if histrange_max - histrange_min > 100 and histrange_min >= 0.0:
				if histrange_min != 0:
					print(np.abs(np.log10(histrange_max)-np.log10(histrange_min))/100.0)
					bins = np.arange(np.log10(histrange_min),np.log10(histrange_max), np.abs(np.log10(histrange_max)-np.log10(histrange_min))/100.0)
				else:
					bins = np.arange(0,np.log10(histrange_max), np.abs(np.log10(histrange_max))/100.0)
				plt.hist(np.log10(mapdata[maskmap_other==1]), bins=bins,label='Mean is ' + ("{0:.2f}".format(np.mean(mapdata[maskmap_other==1]))) + ', median ' + ("{0:.2f}".format(np.median(mapdata[maskmap_other==1]))) + ', std ' + ("{0:.2f}".format(np.std(mapdata[maskmap_other==1]))))
				plt.xlabel('(log10)')
			else:
				bins = np.arange(histrange_min,histrange_max, (histrange_max-histrange_min)/100.0)
				plt.hist(mapdata[maskmap_other==1], bins=bins,label='Mean is ' + ("{0:.2f}".format(np.mean(mapdata[maskmap_other==1]))) + ', median ' + ("{0:.2f}".format(np.median(mapdata[maskmap_other==1]))) + ', std ' + ("{0:.2f}".format(np.std(mapdata[maskmap_other==1]))))
			plt.title(maps[i][0])
			l = plt.legend(prop={'size':9})
			# l = plt.legend(prop={'size':11})
			l.set_zorder(20)
			plt.savefig(maps[i][0].replace('.fits','_histogram.pdf'))
			plt.clf()

		# if pol:
		# 	mapdata = np.sqrt(mapdata[0]**2+mapdata[1]**2)
		if maps[i][1] != nside:
			mapdata = hp.ud_grade(mapdata, nside)

		print(maps[i][5])
		minvals[maps[i][5]] = np.sqrt(np.mean(np.square(mapdata[maskmap_min == 1]))) / maps[i][6]
		maxvals[maps[i][5]] = np.sqrt(np.mean(np.square(mapdata[maskmap_max == 1]))) / maps[i][6]

	print(minvals)
	print(maxvals)

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
	plt.style.use('classic')
	plt.figure(figsize=(7.2, 5.4), dpi=80)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Frequency (GHz)', fontsize=17)
	plt.ylabel('Brightness temperature rms ($\mu$K)', fontsize=17)
	plt.xlim(xmin=minfreq,xmax=maxfreq)
	plt.ylim(ymin=ymin,ymax=ymax)
	plt.xticks(fontsize=17)
	plt.yticks(fontsize=17)
	if pol == False:

		sync_spectrum_min = 0
		if galprop_file != '':
			sync_spectrum_min = syncshifted_comm(x, minvals[0], 4e-3, galprop_freq, galprop_amp)
			sync_spectrum_max = syncshifted_comm(x, maxvals[0], 4e-3, galprop_freq, galprop_amp)
			plt.fill_between(x, sync_spectrum_min, sync_spectrum_max,facecolor='magenta',lw=0,zorder=10,label='Synchrotron')
		elif minvals[10] != 0:
			print(minvals[0])
			sync_spectrum_min = synchrotron(const, x, 22.8, minvals[0], -minvals[10])
			sync_spectrum_max = synchrotron(const, x, 22.8, maxvals[0], -maxvals[10])
			print(sync_spectrum_min)
			# exit(0)
			plt.fill_between(x, sync_spectrum_min, sync_spectrum_max,facecolor='magenta',lw=0,zorder=10,label='Synchrotron')

		freefree_spectrum_min = freefree(const, x, minvals[1], minvals[2], solid_angle, equation=1,comm=1)
		freefree_spectrum_max = freefree(const, x, maxvals[1], maxvals[2], solid_angle, equation=1,comm=1)
		plt.fill_between(x, freefree_spectrum_min, freefree_spectrum_max,facecolor='blue',lw=0,zorder=9,label='Free-free')

		print(minvals[10])
		print(minvals[5])
		ame_spectrum_min = spinningdust_comm(x, minvals[4], minvals[5], spd_amp, spd_freq, 3)
		ame_spectrum_max = spinningdust_comm(x, maxvals[4], maxvals[5], spd_amp, spd_freq, 3)
		# ame_spectrum_min = spinningdust_comm(x, minvals[4], minvals[5], spd_amp, spd_freq, 1) + spinningdust_comm(x, minvals[6], 33.35, spd_amp, spd_freq, 2)
		# ame_spectrum_max = spinningdust_comm(x, maxvals[4], maxvals[5], spd_amp, spd_freq, 1) + spinningdust_comm(x, maxvals[6], 33.35, spd_amp, spd_freq, 2)
		plt.fill_between(x, ame_spectrum_min, ame_spectrum_max,facecolor='gold',lw=0,zorder=9,label="AME")

		cmb_spectrum_min = cmb_comm(const, x, minvals[3])
		cmb_spectrum_max = cmb_comm(const, x, maxvals[3])
		plt.fill_between(x, cmb_spectrum_min, cmb_spectrum_max,facecolor='k',lw=1.0,zorder=11,label='CMB')

		# thermaldust_spectrum_min = thermaldust_comm(const, x, minvals[7], minvals[9], minvals[8])
		# thermaldust_spectrum_max = thermaldust_comm(const, x, maxvals[7], maxvals[9], maxvals[8])
		# plt.fill_between(x, thermaldust_spectrum_min, thermaldust_spectrum_max,facecolor='r',lw=0,zorder=9,label='Thermal dust')

		fd_min = sync_spectrum_min + freefree_spectrum_min + ame_spectrum_min# + thermaldust_spectrum_min
		fd_max = sync_spectrum_max + freefree_spectrum_max + ame_spectrum_max# + thermaldust_spectrum_max
		plt.plot(x, fd_min, 'k--',zorder=12,label='Total foreground')
		plt.plot(x, fd_max, 'k--',zorder=12)
	else:
		sync_spectrum_min = syncshifted_pol_comm(x, np.sqrt(minvals[0]**2+minvals[1]**2), 4e-3, galprop_freq, galprop_amp)
		sync_spectrum_max = syncshifted_pol_comm(x, np.sqrt(maxvals[0]**2+maxvals[1]**2), 4e-3, galprop_freq, galprop_amp)
		plt.fill_between(x, sync_spectrum_min, sync_spectrum_max,facecolor='magenta',lw=0,zorder=10,label="Synchrotron")
		
		print(minvals[2])
		print(minvals[3])
		print(maxvals[2])
		print(maxvals[3])
		minvals[2] = 0.55/np.sqrt(2)
		minvals[3] = 0.55/np.sqrt(2)
		maxvals[2] = 0.7/np.sqrt(2)
		maxvals[3] = 0.7/np.sqrt(2)
		cmb_spectrum_min = cmb_comm(const, x, np.sqrt(minvals[2]**2+minvals[3]**2))
		cmb_spectrum_max = cmb_comm(const, x, np.sqrt(maxvals[2]**2+maxvals[3]**2))
		plt.fill_between(x, cmb_spectrum_min, cmb_spectrum_max,facecolor='k',lw=1.0,zorder=11,label='CMB')

		thermaldust_spectrum_min = thermaldust_comm(const, x, np.sqrt(minvals[4]**2+minvals[5]**2), minvals[7], minvals[6])
		thermaldust_spectrum_max = thermaldust_comm(const, x, np.sqrt(maxvals[4]**2+maxvals[5]**2), maxvals[7], maxvals[6])
		plt.fill_between(x, thermaldust_spectrum_min, thermaldust_spectrum_max,facecolor='r',lw=0,zorder=9,label='Thermal dust')

		fd_min = sync_spectrum_min + thermaldust_spectrum_min
		fd_max = sync_spectrum_max + thermaldust_spectrum_max
		plt.plot(x, fd_min, 'k--',zorder=12,label='Total foreground')
		plt.plot(x, fd_max, 'k--',zorder=12)

	numfreqbands = len(freqbands)
	for i in range (0,numfreqbands):
		if (freqbands[i][0] != freqbands[i][1]):
			plt.axvspan(freqbands[i][0], freqbands[i][1], color=freqbands[i][3], alpha=freqbands[i][4], lw=1.0,zorder=0,label=freqbands[i][2])
		else:
			plt.plot((freqbands[i][0], freqbands[i][1]), (ymin, ymax), freqbands[i][3], alpha=freqbands[i][4], lw=1.0,zorder=0,label=freqbands[i][2],linestyle=freqbands[i][5])

	if legend == True:
		l = plt.legend(prop={'size':9})
		# l = plt.legend(prop={'size':11})
		l.set_zorder(20)
	plt.savefig(outdir+name+'.pdf')
	plt.savefig(outdir+name+'.png')
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
