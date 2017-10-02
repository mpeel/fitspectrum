from smoothmap import smoothmap
import numpy as np


output_resolution = 60.0
output_nside = 256

#directory = '/scratch/nas_cbass/scratch/mpeel/smoothmaps/'
directory = '/mirror/data/mpeel/smoothmaps/'


beamtf_K = np.loadtxt(directory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9dec_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0])
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits',directory+'256_60.00smoothed_wmap9dec_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0])
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',directory+'256_60.00smoothed_wmap9dec_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0])
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',directory+'256_60.00smoothed_wmap9dec_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0])
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',directory+'256_60.00smoothed_wmap9dec_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0])

smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9decbeam_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_K)
smoothmap(directory+'wmap_band_deconv_imap_r9_9yr_Ka_v5.fits',directory+'256_60.00smoothed_wmap9decbeam_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Ka)
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',directory+'256_60.00smoothed_wmap9decbeam_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Q)
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',directory+'256_60.00smoothed_wmap9decbeam_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_V)
smoothmap(directory+'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',directory+'256_60.00smoothed_wmap9decbeam_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_W)

smoothmap(directory+'wmap_band_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK')
smoothmap(directory+'wmap_band_imap_r9_9yr_Ka_v5.fits',directory+'256_60.00smoothed_wmap9_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK')
smoothmap(directory+'wmap_band_imap_r9_9yr_Q_v5.fits',directory+'256_60.00smoothed_wmap9_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK')
smoothmap(directory+'wmap_band_imap_r9_9yr_V_v5.fits',directory+'256_60.00smoothed_wmap9_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK')
smoothmap(directory+'wmap_band_imap_r9_9yr_W_v5.fits',directory+'256_60.00smoothed_wmap9_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK')

smoothmap(directory+'wmap_band_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',windowfunction=beamtf_K)
smoothmap(directory+'wmap_band_imap_r9_9yr_Ka_v5.fits',directory+'256_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',windowfunction=beamtf_Ka)
smoothmap(directory+'wmap_band_imap_r9_9yr_Q_v5.fits',directory+'256_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',windowfunction=beamtf_Q)
smoothmap(directory+'wmap_band_imap_r9_9yr_V_v5.fits',directory+'256_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',windowfunction=beamtf_V)
smoothmap(directory+'wmap_band_imap_r9_9yr_W_v5.fits',directory+'256_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',windowfunction=beamtf_W)

beamtf_p30 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=28)
beamtf_p30 = beamtf_p30[0]
beamtf_p44 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=29)
beamtf_p44 = beamtf_p44[0]
beamtf_p70 = hp.fitsfunc.mrdfits(directory+'LFI_RIMO_R2.50.fits',hdu=30)
beamtf_p70 = beamtf_p70[0]

HFIbeams = fits.open('HFI_RIMO_Beams-100pc_R2.00.fits',)
beamtf_p100 = inputfits[3].data[0][0]
beamtf_p143 = inputfits[4].data[0][0]
beamtf_p217 = inputfits[5].data[0][0]
beamtf_p353 = inputfits[6].data[0][0]
beamtf_p545 = inputfits[7].data[0][0]
beamtf_p857 = inputfits[8].data[0][0]

smoothmap(directory+'LFI_SkyMap_030_1024_R2.01_full.fits',directory+'256_60.00smoothed_PlanckR2full_28.4_1024_2015_KCMBunits.fits', np.sqrt(output_resolution**2-(33.16)**2),nside_out=output_nside,outputmaps=[0,4])
smoothmap(directory+'LFI_SkyMap_030_1024_R2.01_full.fits',directory+'256_60.00smoothed_PlanckR2fullbeam_28.4_1024_2015_KCMBunits.fits', output_resolution,nside_out=output_nside,outputmaps=[0,4],windowfunction=beamtf_p30)


exit()

# smoothmap(directory+'LFI_SkyMap_030_1024_R2.01_full.fits',directory+'256_60.00smoothed_LFI_30_256_PR2.01_full_test1.fits',np.sqrt(output_resolution**2-34.2**2),nside_out=output_nside)#,units_out='mK_RJ')
smoothmap(directory+'LFI_SkyMap_030_1024_R2.01_full.fits',directory+'256_60.00smoothed_LFI_30_256_PR2.01_full_test2a.fits',np.sqrt(output_resolution**2-34.2**2),nside_out=output_nside,outputmaps=[0,4])

# smoothmap(directory+'wmap_band_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9_K_test1.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK')
# smoothmap(directory+'wmap_band_imap_r9_9yr_K_v5.fits',directory+'256_60.00smoothed_wmap9_K_test2.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0])
exit()

# smoothmap(directory+'LFI_SkyMap_030-field-IQU_1024_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2.01_full.fits',nside_out=output_nside,units_out='mK_CMB')
#smoothmap("/Users/mpeel/Desktop/wmap_band_smth_imap_r9_7yr_K_v4.fits","/Users/mpeel/Desktop/512_240.00smoothed_wmap_band_smth_imap_r9_7yr_K_v4.fits",np.sqrt(output_resolution**2-51.3**2),nside_out=output_nside)
# smoothmap(directory+'LFI_SkyMap_030-field-IQU_1024_R2.01_full.fits',directory+'512_60.00smoothed_LFI_SkyMap_30_256_PR2.01_full.fits', np.sqrt(output_resolution**2-32.29**2),nside_out=output_nside)

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
