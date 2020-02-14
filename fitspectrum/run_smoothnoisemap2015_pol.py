from smoothnoisemap import smoothnoisemap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header) # print header
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = [60.0]#,120.0,240.0]
output_nside = [512, 256, 128, 64]
numrealisations=1000
# directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
directory = '/net/nas/proyectos/quijote/mikepeel/'
outdirectory = directory+"wmap9_planck2015_tqu_noise/"

# beamtf_K = np.loadtxt(directory+'wmap9/wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Ka = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Q = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
# beamtf_V = np.loadtxt(directory+'wmap9/wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
# beamtf_W = np.loadtxt(directory+'wmap9/wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))


beamtf_p30 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',28)
beamtf_p44 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',29)
beamtf_p70 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',30)

HFIbeams = fits.open(directory+'planck2015/HFI_RIMO_Beams-100pc_R2.00.fits')
beamtf_p100 = HFIbeams[3].data[0][0]
beamtf_p143 = HFIbeams[4].data[0][0]
beamtf_p217 = HFIbeams[5].data[0][0]
beamtf_p353 = HFIbeams[6].data[0][0]
beamtf_p545 = HFIbeams[7].data[0][0]
beamtf_p857 = HFIbeams[8].data[0][0]

numres = len(output_resolution)

mapnumbers = [4,5,6,7,8,9]

for i in range(0,numres):
	for m in range(0,len(mapnumbers)):
		resolution = "%.2f" % output_resolution[i]

		smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_28.4_1024_2015_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p30,rescale=1000.0)
		smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_44.1_1024_2015_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p44,rescale=1000.0)
		smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_70.4_1024_2015_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_070_2048_R2.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p70,rescale=1000.0)
		# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_100_1024_2015_mKCMBunits', 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100,rescale=1000.0)
		# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_143_1024_2015_mKCMBunits', 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143,rescale=1000.0)
		# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_217_1024_2015_mKCMBunits', 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217,rescale=1000.0)
		# smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_353_1024_2015_mKCMBunits', 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353,rescale=1000.0)
		# try:
		# 	smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_545_1024_2015_MJySrunits', 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545)
		# 	smoothnoisemap(directory+'planck2015/', outdirectory, resolution+'smoothed_PlanckR2fullbeam_857_1024_2015_MJySrunits', 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857)
		# except:
		# 	null = 0
