import healpy as hp
import numpy as np

basedir = '/Users/mpeel/Documents/maps/'
p353_file = basedir+'planck2018_tqu/20.0smoothed_PlanckR3fullbeam_353_2048_2018_mKCMBunits_noisenum2_7_actualvariance.fits'
wgt_file = basedir+'planck2018_tqu_weight/planck2018_hfi_tqu_v0_combine_u_unc.fits'
p353 = hp.read_map(p353_file)
wgt = hp.read_map(wgt_file)
print(np.median(np.sqrt(p353)))
print(np.median(np.sqrt(wgt)))
print(np.median(np.sqrt(p353))/np.median(np.sqrt(wgt)))
