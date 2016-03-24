from smoothmap import smoothmap
import numpy as np


output_resolution = 60.0
output_nside = 512

# WMAP 9-year maps
directory = 'wmap/9yr/'
smoothmap(directory+'wmap_band_iqumap_r9_9yr_K_v5.fits',directory+'512_60.00smoothed_wmap_band_iqumap_r9_9yr_K_v5.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),pol=True,nside_out=output_nside,nobsmap=3)

smoothmap(directory+'wmap_band_iqumap_r9_9yr_Ka_v5.fits',directory+'512_60.00smoothed_wmap_band_iqumap_r9_9yr_Ka_v5.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),pol=True,nside_out=output_nside,nobsmap=3)

smoothmap(directory+'wmap_band_iqumap_r9_9yr_Q_v5.fits',directory+'512_60.00smoothed_wmap_band_iqumap_r9_9yr_Q_v5.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),pol=True,nside_out=output_nside,nobsmap=3)

smoothmap(directory+'wmap_band_iqumap_r9_9yr_V_v5.fits',directory+'512_60.00smoothed_wmap_band_iqumap_r9_9yr_V_v5.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),pol=True,nside_out=output_nside,nobsmap=3)

smoothmap(directory+'wmap_band_iqumap_r9_9yr_W_v5.fits',directory+'512_60.00smoothed_wmap_band_iqumap_r9_9yr_W_v5.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),pol=True,nside_out=output_nside,nobsmap=3)

exit()

# Planck PR2 maps
directory = 'PR2/'
smoothmap(directory+'LFI_SkyMap_030-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
smoothmap(directory+'LFI_SkyMap_044-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
smoothmap(directory+'LFI_SkyMap_070-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_70_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_100-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_100_2048_PR2_full.fits', np.sqrt(output_resolution**2-9.68**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_143-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_143_2048_PR2_full.fits', np.sqrt(output_resolution**2-7.30**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_217-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_217_2048_PR2_full.fits', np.sqrt(output_resolution**2-5.02**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_353-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_353_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.94**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_545-field-Int_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_545_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.83**2),pol=False,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_857-field-Int_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_857_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.64**2),pol=False,nside_out=output_nside)
