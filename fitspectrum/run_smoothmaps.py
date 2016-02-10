from smoothmap import smoothmap
import numpy as np
# Planck PR2 maps
directory = 'PR2/'
output_resolution = 60.0
output_nside = 512

smoothmap(directory+'HFI_SkyMap_100-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_100_2048_PR2_full.fits', np.sqrt(output_resolution**2-9.68**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_143-field-IQU_2048_R2.02_full.fits',directory+'/512_60.00smoothed_HFI_SkyMap_143_2048_PR2_full.fits', np.sqrt(output_resolution**2-7.30**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_217-field-IQU_2048_R2.02_full.fits',directory+'/512_60.00smoothed_HFI_SkyMap_217_2048_PR2_full.fits', np.sqrt(output_resolution**2-5.02**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_353-field-IQU_2048_R2.02_full.fits',directory+'/512_60.00smoothed_HFI_SkyMap_353_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.94**2),pol=True,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_545-field-IQU_2048_R2.02_full.fits',directory+'/512_60.00smoothed_HFI_SkyMap_545_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.83**2),pol=False,nside_out=output_nside)
smoothmap(directory+'HFI_SkyMap_857-field-IQU_2048_R2.02_full.fits',directory+'/512_60.00smoothed_HFI_SkyMap_857_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.64**2),pol=False,nside_out=output_nside)
