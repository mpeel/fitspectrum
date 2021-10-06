#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# For the NPIPE maps, remove the dipole based on the 2015 map they provide
# 
# Version history:
#
# 01-Sep-2020  M. Peel       Started
import os
import healpy as hp
import numpy as np

dipole_file = '/Users/mpeel/Documents/maps/planck_npipe/dipole_nside2048.fits'
mapindir = '/Users/mpeel/Documents/maps/wmap9_planck2020_tqu_v1.4/'
mapoutdir = '/Users/mpeel/Documents/maps/wmap9_planck2020_tqu_v1.4_nodipole/'

nsides = [32, 64, 128, 256, 512, 1024, 2048]

inputlist = os.listdir(mapindir)

dipole = hp.read_map(dipole_file)
dipole_32 = hp.ud_grade(dipole,nside_out=32)
dipole_64 = hp.ud_grade(dipole,nside_out=64)
dipole_128 = hp.ud_grade(dipole,nside_out=128)
dipole_256 = hp.ud_grade(dipole,nside_out=256)
dipole_512 = hp.ud_grade(dipole,nside_out=512)
dipole_1024 = hp.ud_grade(dipole,nside_out=1024)
dipoles = [dipole_32, dipole_64, dipole_128, dipole_256, dipole_512, dipole_1024, dipole]

for file in sorted(inputlist):
	if 'fits' in file[-4:]:
		nside = file.split('_')[0]
		# print(nside)
		if nside.isnumeric():
			if int(nside) in nsides:
				print(file)
				inmap = hp.read_map(mapindir+file,field=None)
				inmap[0] -= dipoles[nsides.index(int(nside))]
				hp.write_map(mapoutdir+file.replace('fullbeam','fullbeamnodp'),inmap)