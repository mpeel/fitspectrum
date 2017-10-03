from smoothnoisemap import smoothnoisemap
import numpy as np

output_resolution = 52.8
output_nside = [512, 256, 128, 64]
numrealisations=1000
# directory = '/Users/mpeel/Desktop/noisetest/'
directory = '/mirror/data/mpeel/smoothmaps'

smoothnoisemap(directory, 'test_wmap_K_'+str(output_resolution), 'wmap_band_imap_r9_9yr_K_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.88*60.0)**2),sigma_0=1.429,nside=output_nside)

exit()

smoothnoisemap(directory, 'wmap_Ka_'+str(output_resolution), 'wmap_band_imap_r9_9yr_Ka_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.66*60.0)**2),sigma_0=1.466,nside=output_nside)
smoothnoisemap(directory, 'wmap_Q_'+str(output_resolution), 'wmap_band_imap_r9_9yr_Q_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.51*60.0)**2),sigma_0=2.188,nside=output_nside)
smoothnoisemap(directory, 'wmap_V_'+str(output_resolution), 'wmap_band_imap_r9_9yr_V_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.35*60.0)**2),sigma_0=3.131,nside=output_nside)
smoothnoisemap(directory, 'wmap_W_'+str(output_resolution), 'wmap_band_imap_r9_9yr_W_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(0.22*60.0)**2),sigma_0=6.544,nside=output_nside)

smoothnoisemap(directory, 'planck30_'+str(output_resolution), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(33.16)**2),nside=output_nside)
smoothnoisemap(directory, 'planck44_'+str(output_resolution), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(28.09)**2),nside=output_nside)
smoothnoisemap(directory, 'planck70_'+str(output_resolution), 'LFI_SkyMap_070_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(13.08)**2),nside=output_nside)
smoothnoisemap(directory, 'planck100_'+str(output_resolution), 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(9.58)**2),nside=output_nside)
smoothnoisemap(directory, 'planck143_'+str(output_resolution), 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(7.18)**2),nside=output_nside)
smoothnoisemap(directory, 'planck217_'+str(output_resolution), 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.87)**2),nside=output_nside)
smoothnoisemap(directory, 'planck353_'+str(output_resolution), 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.7)**2),nside=output_nside)
smoothnoisemap(directory, 'planck545_'+str(output_resolution), 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.73)**2),nside=output_nside)
smoothnoisemap(directory, 'planck857_'+str(output_resolution), 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution**2-(4.51)**2),nside=output_nside)
