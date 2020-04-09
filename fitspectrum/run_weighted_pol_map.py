from weighted_pol_map import *
from astrocode.colourcorrections.fastcc import *

# General settings
nside = 512
index = -3.0

# Settings for this run file
doing_quijote = False

# Settings to pass to the code
use_halfrings = False
use_weights = False
use_reweight_by_rms = True
use_reweight_by_rms_method = 2 # 1 = ricardo, 2 = alberto
use_planck = True
use_cbass = False
freqs = [17,19,11,13,17,19]#11,13,
# normfreq = 10.0
normfreq = 28.4

if doing_quijote:
	# QUIJOTE
	indirectory = '/Users/mpeel/Documents/maps/quijote_201907/smooth/'
	outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/analysetemp/'
	date='201907'
	prefix='mfi'

	# indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
	# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'
	# date='201905'

	# Set up QUIJOTE input
	prefix='half1mfi'
	maps_half1 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

	prefix='half2mfi'
	maps_half2 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

	prefix='mfi'
	maps = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

else:
	date = ''
	maps_half1=[]
	maps_half2=[]
	# Combine Planck+WMAP
	# indirectory = '/Users/mpeel/Documents/maps/quijote_201907/smoothcomb/'
	# outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/analysecomb/'
	# freqs = [17,19,11,13,17,19, 28.4, 44.1, 22.8, 33.0, 40.7]
	# maps.append('512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits')
	# maps.append('512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits')
	# maps.append('512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits')
	# maps.append('512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits')
	# maps.append('512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits')
	# normfreq = 28.4

	prefix='wmap9_planck2018_tqu_TEST3ONLY'
	freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
	rescale_amp = [1.0, 1.0, 1.0, 1.0, 1.0]
	maps = ['512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']
	indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
	outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_weight/'


	# prefix='wmap9_planck2015_tqu_10ghz'
	# freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
	rescale_amp = [1.0, 1.0, 1.0, 1.0, 1.0]
	rescale_variance = rescale_amp.copy()
	# maps = ['512_60.0smoothed_PlanckR2fullbeambpcorr_28.4_256_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeambpcorr_44.1_256_2015_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','../wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']
	# indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2015_tqu/'
	# outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2015_weight/'

	# Planck and WMAP colour corrections
	rescale_amp[0] *= fastcc('30',index+2.0)
	rescale_amp[1] *= fastcc('44',index+2.0)
	rescale_amp[2] *= fastcc('K',index+2.0)
	rescale_amp[3] *= fastcc('Ka',index+2.0)
	rescale_amp[4] *= fastcc('Q',index+2.0)
	print(rescale_amp)

	# Rescale the WMAP variances as the wrong sigma_0 was used to generate them.
	rescale_variance[2] *= (1.435/1.429)**2
	rescale_variance[3] *= (1.472/1.466)**2
	rescale_variance[4] *= (2.197/2.188)**2
	print(rescale_variance)
	# exit()

	apply_extra_mask=np.zeros(len(freqs))
	apply_extra_mask[0] = 1
	apply_extra_mask[1] = 1
	extra_mask='compare_wmap_planck_polmap_mask.fits'

	# prefix='wmap9'
	# freqs = [22.8, 33.0, 40.7]
	# maps = ['512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']

	# prefix='planck2018'
	# freqs = [28.4, 44.1]
	# maps = ['512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits']

	# normfreq = 28.4

weighted_pol_map(nside=nside,indirectory=indirectory,outdirectory=outdirectory,date=date,prefix=prefix,index=index,freqs=freqs,maps=maps,maps_half1=maps_half1,maps_half2=maps_half2,use_halfrings=use_halfrings,use_weights=use_weights,use_reweight_by_rms=use_reweight_by_rms,use_reweight_by_rms_method=use_reweight_by_rms_method,use_planck=use_planck,use_cbass=use_cbass,normfreq=normfreq,rescale_amp=rescale_amp,rescale_variance=rescale_variance,apply_extra_mask=apply_extra_mask,extra_mask=extra_mask)
