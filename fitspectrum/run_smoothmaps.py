from smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits


output_resolution = 60.0
output_nside = [512, 256, 128, 64]

directory = '/mirror/data/mpeel/smoothmaps/'
outdirectory = directory+"to_clive/"
beamtf_K = np.loadtxt(directory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

numnside = len(output_nside)
for i in range(0,numnside):
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits','256_60.00smoothed_wmap9dec_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0])
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits','256_60.00smoothed_wmap9dec_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0])
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits','256_60.00smoothed_wmap9dec_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0])
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits','256_60.00smoothed_wmap9dec_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0])
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits','256_60.00smoothed_wmap9dec_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0])

	smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_K,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits')
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits','256_60.00smoothed_wmap9decbeam_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Ka)
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits','256_60.00smoothed_wmap9decbeam_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Q)
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits','256_60.00smoothed_wmap9decbeam_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_V)
	# smoothmap(directory,directory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits','256_60.00smoothed_wmap9decbeam_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_W)

	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_K_v5.fits','256_60.00smoothed_wmap9_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK')
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_Ka_v5.fits','256_60.00smoothed_wmap9_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK')
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_Q_v5.fits','256_60.00smoothed_wmap9_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK')
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_V_v5.fits','256_60.00smoothed_wmap9_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK')
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_W_v5.fits','256_60.00smoothed_wmap9_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK')

	smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.429,sigma_0_unit='mK',windowfunction=beamtf_K,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits')
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_Ka_v5.fits','256_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=1.466,sigma_0_unit='mK',windowfunction=beamtf_Ka)
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_Q_v5.fits','256_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=2.188,sigma_0_unit='mK',windowfunction=beamtf_Q)
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_V_v5.fits','256_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=3.131,sigma_0_unit='mK',windowfunction=beamtf_V)
	# smoothmap(directory,directory,'wmap_band_imap_r9_9yr_W_v5.fits','256_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside,sigma_0=6.544,sigma_0_unit='mK',windowfunction=beamtf_W)

# For Planck, see run_smoothmaps2.py
exit()