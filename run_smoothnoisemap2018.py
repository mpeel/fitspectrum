from smoothnoisemap import smoothnoisemap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print header # print header
	print data.names # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = [60.0,120.0,240.0]
output_nside = [512, 256, 128, 64]
numrealisations=1000
directory = '/mirror/data/mpeel/smoothmaps/'

# beamtf_K = np.loadtxt(directory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Ka = np.loadtxt(directory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Q = np.loadtxt(directory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
# beamtf_V = np.loadtxt(directory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
# beamtf_W = np.loadtxt(directory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

beamtf_p30 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=28)
beamtf_p30 = beamtf_p30[0]
beamtf_p44 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=29)
beamtf_p44 = beamtf_p44[0]
beamtf_p70 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=30)
beamtf_p70 = beamtf_p70[0]

beamtf_p100 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
print len(beamtf_p100)
beamtf_p143 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

numres = len(output_resolution)
for i in range(0,numres):
	resolution = "%.2f" % output_resolution[i]

	# smoothnoisemap(directory, 'planck30_'+str(output_resolution[i]), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(33.16)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck44_'+str(output_resolution[i]), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(28.09)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck70_'+str(output_resolution[i]), 'LFI_SkyMap_070_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(13.08)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck100_'+str(output_resolution[i]), 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(9.58)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck143_'+str(output_resolution[i]), 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(7.18)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck217_'+str(output_resolution[i]), 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.87)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck353_'+str(output_resolution[i]), 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.7)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck545_'+str(output_resolution[i]), 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.73)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck857_'+str(output_resolution[i]), 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.51)**2),nside=output_nside)

	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits', 'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p30,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits', 'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p44,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_70.4_1024_2018_mKCMBunits', 'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p70,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_100_1024_2018_mKCMBunits', 'HFI_SkyMap_100_2048_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_143_1024_2018_mKCMBunits', 'HFI_SkyMap_143_2048_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_217_1024_2018_mKCMBunits', 'HFI_SkyMap_217_2048_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_353_1024_2018_mKCMBunits', 'HFI_SkyMap_353-psb_2048_R3.00_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353,rescale=1000.0)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_545_1024_2018_MJySrunits', 'HFI_SkyMap_545_2048_R3.00_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545)
	smoothnoisemap(directory, resolution+'smoothed_PlanckR3fullbeam_857_1024_2018_MJySrunits', 'HFI_SkyMap_857_2048_R3.00_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857)
