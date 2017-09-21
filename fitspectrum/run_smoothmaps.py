from smoothmap import smoothmap
import numpy as np


output_resolution = 60.0
output_nside = 256

#directory = '/scratch/nas_cbass/scratch/mpeel/smoothmaps/'
directory = '/mirror/data/mpeel/smoothmaps/'

smoothmap(directory+'LFI_SkyMap_030_1024_R2.01_full.fits',directory+'256_60.00smoothed_LFI_SkyMap_30_256_PR2.01_full.fits',np.sqrt(output_resolution**2-34.2**2),nside_out=output_nside,units_out='mK_RJ')

# WMAP 9-year maps
smoothmap(directory+'wmap_band_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9K.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429)

# smoothmap(directory+'wmap_band_imap_r9_9yr_Ka_v5.fits',directory+'256_60.00smoothed_wmap9Ka.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside,sigma_0=1.466)
exit()


smoothmap(directory+'wmap_band_imap_r9_9yr_Q_v5.fits',directory+'256_60.00smoothed_wmap9Q.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside,sigma_0=2.188)

smoothmap(directory+'wmap_band_imap_r9_9yr_V_v5.fits',directory+'256_60.00smoothed_wmap9V.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside,sigma_0=3.131)

smoothmap(directory+'wmap_band_imap_r9_9yr_W_v5.fits',directory+'256_60.00smoothed_wmap9W.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside,sigma_0=6.544)

exit()


# smoothmap(directory+'LFI_SkyMap_030-field-IQU_1024_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2.01_full.fits',nside_out=output_nside,units_out='mK_CMB')
#smoothmap("/Users/mpeel/Desktop/wmap_band_smth_imap_r9_7yr_K_v4.fits","/Users/mpeel/Desktop/512_240.00smoothed_wmap_band_smth_imap_r9_7yr_K_v4.fits",np.sqrt(output_resolution**2-51.3**2),nside_out=output_nside)
# smoothmap(directory+'LFI_SkyMap_030-field-IQU_1024_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2.01_full.fits', np.sqrt(output_resolution**2-32.29**2),nside_out=output_nside)

exit()

# WMAP 1-year maps
directory = 'wmap/1yr/'
smoothmap(directory+'map_k_imap_yr1_v1.fits',directory+'512_60.00smoothed_map_k_imap_yr1_v1.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),pol=False,nside_out=output_nside,nobsmap=1)

smoothmap(directory+'map_ka_imap_yr1_v1.fits',directory+'512_60.00smoothed_map_ka_imap_yr1_v1.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),pol=False,nside_out=output_nside,nobsmap=1)

smoothmap(directory+'map_q_imap_yr1_v1.fits',directory+'512_60.00smoothed_map_q_imap_yr1_v1.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),pol=False,nside_out=output_nside,nobsmap=1)

smoothmap(directory+'map_v_imap_yr1_v1.fits',directory+'512_60.00smoothed_map_v_imap_yr1_v1.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),pol=False,nside_out=output_nside,nobsmap=1)

smoothmap(directory+'map_w_imap_yr1_v1.fits',directory+'512_60.00smoothed_map_w_imap_yr1_v1.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),pol=False,nside_out=output_nside,nobsmap=1)

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
