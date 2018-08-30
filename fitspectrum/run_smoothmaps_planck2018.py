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
	return data[0][0]

output_resolution = 60.0
output_nside = 256

#directory = '/scratch/nas_cbass/scratch/mpeel/smoothmaps/'
directory = '/mirror/data/mpeel/smoothmaps/'

beamtf_p30 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=28)
beamtf_p30 = beamtf_p30[0]
beamtf_p44 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=29)
beamtf_p44 = beamtf_p44[0]
beamtf_p70 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R3.31.fits',hdu=30)
beamtf_p70 = beamtf_p70[0]


# HFIbeams = fits.open(directory+'HFI_RIMO_Beams_R2.00.fits')
# beamtf_p100 = HFIbeams[3].data[0][0]
# beamtf_p143 = HFIbeams[4].data[0][0]
# beamtf_p217 = HFIbeams[5].data[0][0]
# beamtf_p353 = HFIbeams[6].data[0][0]
# beamtf_p545 = HFIbeams[7].data[0][0]
# beamtf_p857 = HFIbeams[8].data[0][0]

beamtf_p100 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
beamtf_p143 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

smoothmap(directory,directory,'LFI_SkyMap_030_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3full_28.4_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(33.16)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'LFI_SkyMap_044_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3full_44.1_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(28.09)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'LFI_SkyMap_070_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3full_70.4_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(13.08)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_100_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_100_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-9.59**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_143_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_143_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-7.18**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_217_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_217_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-4.87**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_353_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_353_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-4.7**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_545_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_545_2048_2018_MJySrunits.fits', np.sqrt(output_resolution**2-4.73**2),nside_out=output_nside,outputmaps=[0,2])
smoothmap(directory,directory,'HFI_SkyMap_857_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3full_857_2048_2018_MJySrunits.fits', np.sqrt(output_resolution**2-4.51**2),nside_out=output_nside,outputmaps=[0,2])

smoothmap(directory,directory,'LFI_SkyMap_030_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p30)
smoothmap(directory,directory,'LFI_SkyMap_044_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_44.1_1024_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p44)
smoothmap(directory,directory,'LFI_SkyMap_070_1024_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_70.4_1024_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p70)
smoothmap(directory,directory,'HFI_SkyMap_100_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_100_2048_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p100)
smoothmap(directory,directory,'HFI_SkyMap_143_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_143_2048_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p143)
smoothmap(directory,directory,'HFI_SkyMap_217_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_217_2048_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p217)
smoothmap(directory,directory,'HFI_SkyMap_353_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_353_2048_2018_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p353)
smoothmap(directory,directory,'HFI_SkyMap_545_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_545_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,2],windowfunction=beamtf_p545)
smoothmap(directory,directory,'HFI_SkyMap_857_2048_R3.00_full.fits','256_60.00smoothed_PlanckR3fullbeam_857_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,2],windowfunction=beamtf_p857)


exit()