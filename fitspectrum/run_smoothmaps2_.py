from smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

output_resolution = 60.0
output_nside = 256

#directory = '/scratch/nas_cbass/scratch/mpeel/smoothmaps/'
directory = '/mirror/data/mpeel/smoothmaps/'

beamtf_p30 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=28)
beamtf_p30 = beamtf_p30[0]
beamtf_p44 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=29)
beamtf_p44 = beamtf_p44[0]
beamtf_p70 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=30)
beamtf_p70 = beamtf_p70[0]

HFIbeams = fits.open(directory+'HFI_RIMO_Beams-100pc_R2.00.fits')
beamtf_p100 = HFIbeams[3].data[0][0]
beamtf_p143 = HFIbeams[4].data[0][0]
beamtf_p217 = HFIbeams[5].data[0][0]
beamtf_p353 = HFIbeams[6].data[0][0]
beamtf_p545 = HFIbeams[7].data[0][0]
beamtf_p857 = HFIbeams[8].data[0][0]

smoothmap(directory,directory,'LFI_SkyMap_030_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2full_28.4_1024_2015_KCMBunits.fits', np.sqrt(output_resolution**2-(33.16)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'LFI_SkyMap_044_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2full_44.1_1024_2015_KCMBunits.fits', np.sqrt(output_resolution**2-(28.09)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'LFI_SkyMap_070_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2full_70.4_1024_2015_KCMBunits.fits', np.sqrt(output_resolution**2-(13.08)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_100_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_100_2048_2015_KCMBunits.fits', np.sqrt(output_resolution**2-9.59**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_143_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_143_2048_2015_KCMBunits.fits', np.sqrt(output_resolution**2-7.18**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_217_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_217_2048_2015_KCMBunits.fits', np.sqrt(output_resolution**2-4.87**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_353_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_353_2048_2015_KCMBunits.fits', np.sqrt(output_resolution**2-4.7**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory,directory,'HFI_SkyMap_545_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_545_2048_2015_MJySrunits.fits', np.sqrt(output_resolution**2-4.73**2),nside_out=output_nside,outputmaps=[0,2])
smoothmap(directory,directory,'HFI_SkyMap_857_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2full_857_2048_2015_MJySrunits.fits', np.sqrt(output_resolution**2-4.51**2),nside_out=output_nside,outputmaps=[0,2])

smoothmap(directory,directory,'LFI_SkyMap_030_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2fullbeam_28.4_1024_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p30)
smoothmap(directory,directory,'LFI_SkyMap_044_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2fullbeam_44.1_1024_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p44)
smoothmap(directory,directory,'LFI_SkyMap_070_1024_R2.01_full.fits','256_60.00smoothed_PlanckR2fullbeam_70.4_1024_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p70)
smoothmap(directory,directory,'HFI_SkyMap_100_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_100_2048_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p100)
smoothmap(directory,directory,'HFI_SkyMap_143_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_143_2048_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p143)
smoothmap(directory,directory,'HFI_SkyMap_217_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_217_2048_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p217)
smoothmap(directory,directory,'HFI_SkyMap_353_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_353_2048_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p353)
smoothmap(directory,directory,'HFI_SkyMap_545_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_545_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,2],windowfunction=beamtf_p545)
smoothmap(directory,directory,'HFI_SkyMap_857_2048_R2.02_full.fits','256_60.00smoothed_PlanckR2fullbeam_857_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,2],windowfunction=beamtf_p857)


exit()