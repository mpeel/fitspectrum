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
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = 60.0
output_nside = [512, 256, 128, 64]

#directory = '/scratch/nas_cbass/scratch/mpeel/smoothmaps/'
directory = '/mirror/data/mpeel/smoothmaps/'
outdirectory = directory+"to_clive_R3_CMB/"

# beamtf_K = np.loadtxt(directory+'wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Ka = np.loadtxt(directory+'wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
# beamtf_Q = np.loadtxt(directory+'wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
# beamtf_V = np.loadtxt(directory+'wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
# beamtf_W = np.loadtxt(directory+'wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

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
print len(beamtf_p100)
beamtf_p143 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

numnside = len(output_nside)
for i in range(0,numnside):
	# Smooth CMB maps
	smoothmap(directory,directory,'COM_CMB_IQU-commander_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBCommander_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	smoothmap(directory,directory,'COM_CMB_IQU-nilc_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBNILC_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	smoothmap(directory,directory,'COM_CMB_IQU-sevem_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSevem_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	smoothmap(directory,directory,'COM_CMB_IQU-smica_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	smoothmap(directory,directory,'COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica-nosz_0.0_2048_2018_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')

	subtractmaps = [str(output_nside[i])+'_60.00smoothed_PlanckCMBCommander_0.0_2048_2018_mKCMBunits.fits', str(output_nside[i])+'_60.00smoothed_PlanckCMBNILC_0.0_2048_2018_mKCMBunits.fits', str(output_nside[i])+'_60.00smoothed_PlanckCMBSevem_0.0_2048_2018_mKCMBunits.fits', str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica_0.0_2048_2018_mKCMBunits.fits', str(output_nside[i])+'_60.00smoothed_PlanckCMBSmica-nosz_0.0_2048_2018_mKCMBunits.fits','']
	subtractmaps_name = ['CMBcommandersub', 'CMBNILCsub', 'CMBSevemsub', 'CMBSmicasub', 'CMBSmicanoszsub','']
	numsubtract = 6
	for j in range(0,numsubtract):
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0])

		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_K,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Ka,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Q,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_V,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_W,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])

		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK')

		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam'+subtractmaps_name[j]+'_22.8_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',windowfunction=beamtf_K,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam'+subtractmaps_name[j]+'_33.0_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',windowfunction=beamtf_Ka,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam'+subtractmaps_name[j]+'_40.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',windowfunction=beamtf_Q,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam'+subtractmaps_name[j]+'_60.7_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',windowfunction=beamtf_V,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9beam'+subtractmaps_name[j]+'_93.5_512_2013_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',windowfunction=beamtf_W,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])

		# smoothmap(directory,outdirectory,'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_28.4_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(33.16)**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_44.1_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(28.09)**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_70.4_1024_2018_KCMBunits.fits', np.sqrt(output_resolution**2-(13.08)**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_100_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_100_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-9.59**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_143_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_143_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-7.18**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_217_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_217_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-4.87**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_353_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_353_2048_2018_KCMBunits.fits', np.sqrt(output_resolution**2-4.7**2),nside_out=output_nside[i],outputmaps=[0,4])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_545_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_545_2048_2018_MJySrunits.fits', np.sqrt(output_resolution**2-4.73**2),nside_out=output_nside[i],outputmaps=[0,2])
		# smoothmap(directory,outdirectory,'HFI_SkyMap_857_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3full_857_2048_2018_MJySrunits.fits', np.sqrt(output_resolution**2-4.51**2),nside_out=output_nside[i],outputmaps=[0,2])

		smoothmap(directory,outdirectory,'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_28.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_28.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_44.1_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_44.1_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_70.4_1024_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_70.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'HFI_SkyMap_100_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_100_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_100_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'HFI_SkyMap_143_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_143_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_143_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'HFI_SkyMap_217_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_217_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_217_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		smoothmap(directory,outdirectory,'HFI_SkyMap_353-psb_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_353_2048_2018_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_353_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",units_out='mKCMB',subtractmap=subtractmaps[j])
		if j == 4:
			smoothmap(directory,outdirectory,'HFI_SkyMap_545_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_545_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p545,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_545_1024_2015_MJySrunits_variance.fits',appendmapname="II_cov",appendmapunit="(MJySr)^2")
			smoothmap(directory,outdirectory,'HFI_SkyMap_857_2048_R3.00_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR3fullbeam'+subtractmaps_name[j]+'_857_2048_2018_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p857,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_857_1024_2015_MJySrunits_variance.fits',appendmapname="II_cov",appendmapunit="(MJySr)^2")
exit()