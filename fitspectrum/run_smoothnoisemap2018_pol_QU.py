from smoothnoisemap import *
import numpy as np
import healpy as hp
import astropy.io.fits as fits

from smoothmap import *
import numba
import matplotlib.pyplot as plt
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix

def get_hfi_beam(FITSfile):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print(header) # print header
	print(data.names) # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][0]
	return newdata

n = 10000
x, y = np.random.multivariate_normal([0,0], [[1, 0.9], [0.9, 1]], n).T
plt.plot(x,y,'.')
plt.savefig('noisetest2.png')
# plt.clf()

# vals2 = np.asarray([[vals[0:100], vals[100:200]], [vals[200:300], vals[300:400]]]).T

# exit()
A = np.asarray([[1.0, 0.9], [0.9, 1.0]])
B = np.asarray([A]*n)
C = np.linalg.cholesky(B)
vals = np.random.normal(0,1, size = (n,2))
x,y = np.einsum('nij,nj->ni', C, vals).T
# vals = np.random.normal(scale=1.0, size=100)
# vals2 = np.asarray([[vals, vals], [vals, vals]]).T
# arr = vals2 * C
# print(np.shape(C))
# print(np.shape(vals))
# print(np.shape(vals2))
# arr = lsqr(C,vals2)
# arr = np.matmul(np.asarray(C[0]),np.asarray(vals))
# x = arr[:,0,0]
# y = arr[:,0,0]
plt.plot(x,y,'.')
plt.savefig('noisetest.png')
exit()

output_resolution = [60.0]#,120.0,240.0]
output_nside = [512, 256, 128, 64]
numrealisations=1000
# directory = '/Users/mpeel/Documents/maps/'
directory = '/scratch1/mpeel/maps/'
outdirectory = directory+"wmap9_planck2018_tqu_noise_QU/"

nside = 8
npix = 12*nside*nside

cov = [[np.ones(npix), np.ones(npix)*0.01], [-np.ones(npix)*0.01, np.ones(npix)]]
testmap = noiserealisation_QU(np.ones(npix), np.ones(npix), np.ones(npix)*0.01, npix)
print(testmap)
print(np.shape(testmap))
exit()

beamtf_K = np.loadtxt(directory+'wmap9/wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))
beamtf_Ka = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Ka1_9yr_v5p1.txt',usecols=(1,))
beamtf_Q = np.loadtxt(directory+'wmap9/wmap_ampl_bl_Q1_9yr_v5p1.txt',usecols=(1,))
beamtf_V = np.loadtxt(directory+'wmap9/wmap_ampl_bl_V1_9yr_v5p1.txt',usecols=(1,))
beamtf_W = np.loadtxt(directory+'wmap9/wmap_ampl_bl_W1_9yr_v5p1.txt',usecols=(1,))

beamtf_p30 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',28)
beamtf_p44 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',29)
beamtf_p70 = get_beam(directory+'planck2018/LFI_RIMO_R3.31.fits',30)

beamtf_p100 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_100x100.fits')
beamtf_p143 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_143x143.fits')
beamtf_p217 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_217x217.fits')
beamtf_p353 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_353x353.fits')
beamtf_p545 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_545x545.fits')
beamtf_p857 = get_hfi_beam(directory+'planck2018/BeamWf_HFI_R3.01/Bl_T_R3.01_fullsky_857x857.fits')

numres = len(output_resolution)
mapnumbers = [0,1,2,3]
hdu=2
for i in range(0,numres):
	for m in range(0,len(mapnumbers)):
		resolution = "%.2f" % output_resolution[i]

		# if m == 0:
		# 	sigma_0 = 1.429
		# else:
		# 	sigma_0 = 1.435
		# smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_22.8_512_20132018_mKCMBunits_'+str(mapnumbers[m])+'.fits', 'wmap_band_iqumap_r9_9yr_K_v5.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_K,sigma_0=sigma_0,nside=output_nside,hdu=hdu)
		# if m == 0:
		# 	sigma_0 = 1.466
		# else:
		# 	sigma_0 = 1.472
		# smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_33.0_512_20132018_mKCMBunits_'+str(mapnumbers[m])+'.fits', 'wmap_band_iqumap_r9_9yr_Ka_v5.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_Ka,sigma_0=sigma_0,nside=output_nside,hdu=hdu)
		# if m == 0:
		# 	sigma_0 = 2.188
		# else:
		# 	sigma_0 = 2.197
		# smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_40.7_512_20132018_mKCMBunits_'+str(mapnumbers[m])+'.fits', 'wmap_band_iqumap_r9_9yr_Q_v5.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_Q,sigma_0=sigma_0,nside=output_nside,hdu=hdu)
		# if m == 0:
		# 	sigma_0 = 3.131
		# else:
		# 	sigma_0 = 3.141
		# smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_60.7_512_20132018_mKCMBunits_'+str(mapnumbers[m])+'.fits', 'wmap_band_iqumap_r9_9yr_V_v5.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_V,sigma_0=sigma_0,nside=output_nside,hdu=hdu)
		# if m == 0:
		# 	sigma_0 = 6.544
		# else:
		# 	sigma_0 = 6.560
		# smoothnoisemap(directory+'/wmap9/', outdirectory, str(output_resolution[i])+'smoothed_wmap9beam_93.5_512_20132018_mKCMBunits_'+str(mapnumbers[m])+'.fits', 'wmap_band_iqumap_r9_9yr_W_v5.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],windowfunction=beamtf_W,sigma_0=sigma_0,nside=output_nside,hdu=hdu)

mapnumbers = [4,5,6,7,8,9]
mapnumbers = [5,6,7,8,9]

for i in range(0,numres):
	for m in range(0,len(mapnumbers)):
		resolution = "%.2f" % output_resolution[i]

		smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p30,rescale=1000.0)
		smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p44,rescale=1000.0)
		smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_70.4_1024_2018_mKCMBunits_'+str(mapnumbers[m]), 'LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p70,rescale=1000.0)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_100_1024_2018_mKCMBunits', 'HFI_SkyMap_100_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p100,rescale=1000.0)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_143_1024_2018_mKCMBunits', 'HFI_SkyMap_143_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p143,rescale=1000.0)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_217_1024_2018_mKCMBunits', 'HFI_SkyMap_217_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p217,rescale=1000.0)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_353_1024_2018_mKCMBunits', 'HFI_SkyMap_353-psb_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p353,rescale=1000.0)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_545_1024_2018_MJySrunits', 'HFI_SkyMap_545_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p545)
		# smoothnoisemap(directory+'planck2018/', outdirectory, resolution+'smoothed_PlanckR3fullbeam_857_1024_2018_MJySrunits', 'HFI_SkyMap_857_2048_R3.01_full.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=beamtf_p857)
