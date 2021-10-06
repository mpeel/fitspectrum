from smoothmap import *
from smoothnoisemap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header)
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = 60.0#20.0
output_nside = [2048]#, 512]#, 256, 128, 64]

directory = '/Users/mpeel/Documents/maps/'
# directory = '/scratch1/mpeel/maps/'
# outdirectory = directory+"planck2018_tqu/"
outdirectory = '/Volumes/Toshiba5TB2/maps/planck2018_tqu/'
# beamtf_p100 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
beamtf_p143 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')

numnside = len(output_nside)
mapnumbers = [5,6,7,8,9]
for i in range(0,numnside):

		# smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_100_2048_R3.01_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_100_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j])#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam_100_1024_2018_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2")
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_143_2048_R3.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam_143_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,units_out='mKCMB')
		for m in range(0,len(mapnumbers)):
			try:
				smoothnoisemap(directory+'planck2018/', outdirectory, str(output_resolution)+'smoothed_PlanckR3fullbeam_143_2048_2018_mKCMBunits_noisenum2_'+str(mapnumbers[m]), 'HFI_SkyMap_143_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=2,fwhm=output_resolution,nside=[output_nside[i]],windowfunction=beamtf_p143,rescale=1000.0)
			except:
				null = 0
		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_217_2048_R3.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam_217_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,units_out='mKCMB')#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam_217_1024_2018_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
		for m in range(0,len(mapnumbers)):
			try:
				smoothnoisemap(directory+'planck2018/', outdirectory, str(output_resolution)+'smoothed_PlanckR3fullbeam_217_2048_2018_mKCMBunits_noisenum2_'+str(mapnumbers[m]), 'HFI_SkyMap_217_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=2,fwhm=output_resolution,nside=[output_nside[i]],windowfunction=beamtf_p143,rescale=1000.0)
			except:
				null = 0

		smoothmap(directory+'planck2018/',outdirectory,'HFI_SkyMap_353-psb_2048_R3.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR3fullbeam_353_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,units_out='mKCMB')

		for m in range(0,len(mapnumbers)):
			try:
				smoothnoisemap(directory+'planck2018/', outdirectory, str(output_resolution)+'smoothed_PlanckR3fullbeam_353_2048_2018_mKCMBunits_noisenum2_'+str(mapnumbers[m]), 'HFI_SkyMap_353-psb_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=2,fwhm=output_resolution,nside=[output_nside[i]],windowfunction=beamtf_p143,rescale=1000.0)
			except:
				null = 0
# EOF