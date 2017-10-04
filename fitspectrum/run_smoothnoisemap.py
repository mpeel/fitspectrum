from smoothnoisemap import smoothnoisemap
import numpy as np

output_resolution = [60.0,120.0,240.0]
output_nside = [512, 256, 128, 64]
numrealisations=1000
directory = '/mirror/data/mpeel/smoothmaps'

beamtf_K = np.loadtxt(directory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

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

numres = len(output_resolution[i])
for i in range(0,numres):
	smoothnoisemap(directory, 'wmap_K_'+str(output_resolution[i][]), 'wmap_band_imap_r9_9yr_K_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(0.88*60.0)**2),sigma_0=1.429,nside=output_nside)
	# smoothnoisemap(directory, 'wmap_Ka_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_Ka_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(0.66*60.0)**2),sigma_0=1.466,nside=output_nside)
	# smoothnoisemap(directory, 'wmap_Q_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_Q_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(0.51*60.0)**2),sigma_0=2.188,nside=output_nside)
	# smoothnoisemap(directory, 'wmap_V_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_V_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(0.35*60.0)**2),sigma_0=3.131,nside=output_nside)
	# smoothnoisemap(directory, 'wmap_W_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_W_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(0.22*60.0)**2),sigma_0=6.544,nside=output_nside)

	# smoothnoisemap(directory, 'planck30_'+str(output_resolution[i]), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(33.16)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck44_'+str(output_resolution[i]), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(28.09)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck70_'+str(output_resolution[i]), 'LFI_SkyMap_070_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(13.08)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck100_'+str(output_resolution[i]), 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(9.58)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck143_'+str(output_resolution[i]), 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(7.18)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck217_'+str(output_resolution[i]), 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.87)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck353_'+str(output_resolution[i]), 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.7)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck545_'+str(output_resolution[i]), 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.73)**2),nside=output_nside)
	# smoothnoisemap(directory, 'planck857_'+str(output_resolution[i]), 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=np.sqrt(output_resolution[i]**2-(4.51)**2),nside=output_nside)

	smoothnoisemap(directory, 'wmap_K_beam_'+str(output_resolution[i][]), 'wmap_band_imap_r9_9yr_K_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.429,nside=output_nside,windowfunction=beamtf_K)
	# smoothnoisemap(directory, 'wmap_Ka_beam_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_Ka_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.466,nside=output_nside,windowfunction=beamtf_Ka)
	# smoothnoisemap(directory, 'wmap_Q_beam_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_Q_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=2.188,nside=output_nside,windowfunction=beamtf_Q)
	# smoothnoisemap(directory, 'wmap_V_beam_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_V_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=3.131,nside=output_nside,windowfunction=beamtf_V)
	# smoothnoisemap(directory, 'wmap_W_beam_'+str(output_resolution[i]), 'wmap_band_imap_r9_9yr_W_v5.fits',mapnumber=1,numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=6.544,nside=output_nside,windowfunction=beamtf_W)

	# smoothnoisemap(directory, 'planck30_beam_'+str(output_resolution[i]), 'LFI_SkyMap_030_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,output_resolution[i],nside=output_nside,windowfunction=beamtf_p30)
	# smoothnoisemap(directory, 'planck44_beam_'+str(output_resolution[i]), 'LFI_SkyMap_044_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,output_resolution[i],nside=output_nside,windowfunction=beamtf_p44)
	# smoothnoisemap(directory, 'planck70_beam_'+str(output_resolution[i]), 'LFI_SkyMap_070_1024_R2.01_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p70)
	# smoothnoisemap(directory, 'planck100_beam_'+str(output_resolution[i]), 'HFI_SkyMap_100_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100)
	# smoothnoisemap(directory, 'planck143_beam_'+str(output_resolution[i]), 'HFI_SkyMap_143_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143)
	# smoothnoisemap(directory, 'planck217_beam_'+str(output_resolution[i]), 'HFI_SkyMap_217_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217)
	# smoothnoisemap(directory, 'planck353_beam_'+str(output_resolution[i]), 'HFI_SkyMap_353_2048_R2.02_full.fits',mapnumber=4,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353)
	# smoothnoisemap(directory, 'planck545_beam_'+str(output_resolution[i]), 'HFI_SkyMap_545_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545)
	# smoothnoisemap(directory, 'planck857_beam_'+str(output_resolution[i]), 'HFI_SkyMap_857_2048_R2.02_full.fits',mapnumber=2,numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857)
