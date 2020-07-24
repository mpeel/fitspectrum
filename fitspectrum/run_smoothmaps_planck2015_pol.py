from smoothmap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header)
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

output_resolution = 60.0
output_nside = [2048, 1024, 512, 256, 128, 64, 32]
smoothvariance = False

# directory = '/Users/mpeel/Documents/maps/'
directory = '/scratch1/mpeel/maps/'
outdirectory = directory+"wmap9_planck2015_tqu_v1.4/"

beamtf_K = np.loadtxt(directory+'wmap9/wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap9/wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap9/wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

beamtf_p30 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',28)
beamtf_p44 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',29)
beamtf_p70 = get_beam(directory+'planck2015/LFI_RIMO_R2.50.fits',30)

HFIbeams = fits.open(directory+'planck2015/HFI_RIMO_Beams-100pc_R2.00.fits')
beamtf_p100 = HFIbeams[3].data[0][0]
beamtf_p143 = HFIbeams[4].data[0][0]
beamtf_p217 = HFIbeams[5].data[0][0]
beamtf_p353 = HFIbeams[6].data[0][0]
beamtf_p545 = HFIbeams[7].data[0][0]
beamtf_p857 = HFIbeams[8].data[0][0]

numnside = len(output_nside)
for i in range(0,numnside):
	# Smooth CMB maps
	# smoothmap(directory+'planck2015/',directory+'planck2015/','COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBCommander_0.0_2048_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	# smoothmap(directory+'planck2015/',directory+'planck2015/','COM_CMB_IQU-nilc-field-Int_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBNILC_0.0_2048_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	# smoothmap(directory+'planck2015/',directory+'planck2015/','COM_CMB_IQU-sevem-field-Int_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBSevem_0.0_2048_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')
	# smoothmap(directory+'planck2015/',directory+'planck2015/','COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBSmica_0.0_2048_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-5.0**2),nside_out=output_nside[i],units_out='mKCMB')

	# subtractmaps = [str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBCommander_0.0_2048_2015_mKCMBunits.fits', str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBNILC_0.0_2048_2015_mKCMBunits.fits', str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBSevem_0.0_2048_2015_mKCMBunits.fits', str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckCMBSmica_0.0_2048_2015_mKCMBunits.fits','']
	# subtractmaps_name = ['CMBcommandersub', 'CMBNILCsub', 'CMBSevemsub', 'CMBSmicasub','']
	subtractmaps = ['']
	subtractmaps_name = ['']
	numsubtract = len(subtractmaps)
	for j in range(0,numsubtract):
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9dec_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0])

		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_22.8_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_K,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_33.0_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Ka,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_40.7_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_Q,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_60.7_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_V,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])
		# smoothmap(directory,outdirectory,'wmap_band_smth_deconv_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9decbeam'+subtractmaps_name[j]+'_93.5_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',nosmooth=[0],windowfunction=beamtf_W,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2",subtractmap=subtractmaps[j])

		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_K_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_22.8_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.88*60.0)**2),nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_33.0_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.66*60.0)**2),nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_40.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.51*60.0)**2),nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_V_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_60.7_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.35*60.0)**2),nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK')
		# smoothmap(directory,outdirectory,'wmap_band_imap_r9_9yr_W_v5.fits',str(output_nside[i])+'_60.00smoothed_wmap9_93.5_512_2013_mKCMBunits.fits', np.sqrt(output_resolution**2-(0.22*60.0)**2),nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK')

		if output_nside[i] <= 512:
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_K_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_22.8_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.429,sigma_0_unit='mK',windowfunction=beamtf_K,subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_22.8_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_Ka_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_33.0_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=1.466,sigma_0_unit='mK',windowfunction=beamtf_Ka,subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_33.0_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_Q_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_40.7_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=2.188,sigma_0_unit='mK',windowfunction=beamtf_Q,subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_40.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_V_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_60.7_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=3.131,sigma_0_unit='mK',windowfunction=beamtf_V,subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_60.7_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
			smoothmap(directory+'wmap9/',outdirectory,'wmap_band_iqumap_r9_9yr_W_v5.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_wmap9beam'+subtractmaps_name[j]+'_93.5_512_20132015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],sigma_0=6.544,sigma_0_unit='mK',windowfunction=beamtf_W,subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_wmap9beam_93.5_512_2013_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"


		# NB: These are 1° smoothed already!!
		if output_nside[i] <= 256:
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_030-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_28.4_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_28.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"

			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_044-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_44.1_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_44.1_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_070-BPassCorrected_0256_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeambpcorr'+subtractmaps_name[j]+'_70.4_256_2015_mKCMBunits.fits', np.sqrt(output_resolution**2-(60.0)**2),nside_out=output_nside[i],units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_70.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"

		if output_nside[i] <= 1024:
			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_030_1024_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_28.4_1024_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p30,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_28.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"

			smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_044_1024_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_44.1_1024_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p44,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_44.1_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
		
		smoothmap(directory+'planck2015/',outdirectory,'LFI_SkyMap_070_2048_R2.01_full.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_70.4_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p70,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_70.4_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"

		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_100_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_100_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p100,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_100_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2")
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_143_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_143_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p143,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_143_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_217_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_217_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p217,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_217_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"
		smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_353_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_353_2048_2015_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p353,units_out='mKCMB',subtractmap=subtractmaps[j],smoothvariance=smoothvariance)#,maxnummaps=1,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_353_1024_2015_mKCMBunits_variance.fits',appendmapname="II_cov",appendmapunit="(mK)^2"

		# Only want these for the non-CMB subtracted sets.
		if subtractmaps_name[j] == '':
			smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_545_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_545_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p545,smoothvariance=smoothvariance)#,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_545_1024_2015_MJySrunits_variance.fits',appendmapname="II_cov",appendmapunit="(MJySr)^2")
			smoothmap(directory+'planck2015/',outdirectory,'HFI_SkyMap_857_2048_R2.02_full.fits',str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam'+subtractmaps_name[j]+'_857_2048_2015_MJySrunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=beamtf_p857,smoothvariance=smoothvariance)#,appendmap=directory+str(output_nside[i])+'_60.00smoothed_PlanckR2fullbeam_857_1024_2015_MJySrunits_variance.fits',appendmapname="II_cov",appendmapunit="(MJySr)^2")

# EOF