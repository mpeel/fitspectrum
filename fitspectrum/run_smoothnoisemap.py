from smoothnoisemap import smoothnoisemap
import numpy as np

output_resolution = 0.88*60.0
output_nside = 512
numrealisations=100
# WMAP 9-year maps
# directory = '/Users/mpeel/Desktop/noisetest/'
directory = '/mirror/data/mpeel/smoothmaps'

smoothnoisemap(directory, 'planck30_'+str(output_resolution), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(33.16)**2))

exit()

smoothnoisemap(directory, 'wmap_K_'+str(output_resolution), 'wmap_band_imap_r9_9yr_K_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.88*60.0)**2),sigma_0=1.429)
smoothnoisemap(directory, 'wmap_Ka_'+str(output_resolution), 'wmap_band_imap_r9_9yr_Ka_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.66*60.0)**2),sigma_0=1.466)
smoothnoisemap(directory, 'wmap_Q_'+str(output_resolution), 'wmap_band_imap_r9_9yr_Q_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.51*60.0)**2),sigma_0=2.188)
smoothnoisemap(directory, 'wmap_V_'+str(output_resolution), 'wmap_band_imap_r9_9yr_V_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.35*60.0)**2),sigma_0=3.131)
smoothnoisemap(directory, 'wmap_W_'+str(output_resolution), 'wmap_band_imap_r9_9yr_W_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.22*60.0)**2),sigma_0=6.544)

smoothnoisemap(directory, 'planck30_'+str(output_resolution), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(33.16)**2))
smoothnoisemap(directory, 'planck44_'+str(output_resolution), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(28.09)**2))
smoothnoisemap(directory, 'planck70_'+str(output_resolution), 'LFI_SkyMap_070_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(13.08)**2))
smoothnoisemap(directory, 'planck100_'+str(output_resolution), 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(9.58)**2))
smoothnoisemap(directory, 'planck143_'+str(output_resolution), 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(7.18)**2))
smoothnoisemap(directory, 'planck217_'+str(output_resolution), 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.87)**2))
smoothnoisemap(directory, 'planck353_'+str(output_resolution), 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.7)**2))
smoothnoisemap(directory, 'planck545_'+str(output_resolution), 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.73)**2))
smoothnoisemap(directory, 'planck857_'+str(output_resolution), 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.51)**2))



# exit()

# # Planck PR2 maps
# directory = 'PR2/'
# smoothmap(directory+'LFI_SkyMap_030-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'LFI_SkyMap_044-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_44_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'LFI_SkyMap_070-BPassCorrected-field-IQU_0256_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_70_256_PR2_full.fits', np.sqrt(output_resolution**2-60.0**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_100-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_100_2048_PR2_full.fits', np.sqrt(output_resolution**2-9.68**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_143-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_143_2048_PR2_full.fits', np.sqrt(output_resolution**2-7.30**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_217-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_217_2048_PR2_full.fits', np.sqrt(output_resolution**2-5.02**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_353-field-IQU_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_353_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.94**2),pol=True,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_545-field-Int_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_545_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.83**2),pol=False,nside_out=output_nside)
# smoothmap(directory+'HFI_SkyMap_857-field-Int_2048_R2.02_full.fits',directory+'512_60.00smoothed_HFI_SkyMap_857_2048_PR2_full.fits', np.sqrt(output_resolution**2-4.64**2),pol=False,nside_out=output_nside)
