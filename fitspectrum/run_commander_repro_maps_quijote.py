from commander_repro_maps import commander_repro_maps
from smoothmap import smoothmap
import numpy as np

# smoothmap('commander2015/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits','commander2015/256_60.00smoothed_COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits', np.sqrt(60.0**2-40.0**2),pol=True,nside_out=256)
# smoothmap('commander2015/COM_CompMap_DustPol-commander_1024_R2.00.fits','commander2015/1024_60.00smoothed_COM_CompMap_DustPol-commander_1024_R2.00.fits', np.sqrt(60.0**2-10.0**2),pol=True,nside_out=1024,maxnummaps=2)
# smoothmap('commander2015/COM_CMB_IQU-commander_1024_R2.02_full.fits','commander2015/1024_60.00smoothed_COM_CMB_IQU-commander_1024_R2.02_full.fits', np.sqrt(60.0**2-10.0**2),pol=True,nside_out=1024)

maps_I = [['commander2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits',256,60.0,1,0,0,1.0],
['commander2015/COM_CompMap_freefree-commander_0256_R2.00.fits',256,60.0,1,0,1,1.0],
['commander2015/COM_CompMap_freefree-commander_0256_R2.00.fits',256,60.0,1,3,2,1.0],
['commander2015/COM_CompMap_CMB-commander_0256_R2.00.fits',256,60.0,1,0,3,1.0],
['commander2015/COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,1,0,4,1.0],
['commander2015/COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,1,3,5,1.0],
['commander2015/COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,2,0,6,1.0],
['commander2015/COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,0,7,1.0],
['commander2015/COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,3,8,1.0],
['commander2015/COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,6,9,1.0]]

maps_P = [['commander2015/256_60.00smoothed_COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits',256,60.0,1,0,0,1.0],
['commander2015/256_60.00smoothed_COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits',256,60.0,1,1,1,1.0],
['commander2015/1024_60.00smoothed_COM_CMB_IQU-commander_1024_R2.02_full.fits',1024,60.0,1,1,2,1e-6],
['commander2015/1024_60.00smoothed_COM_CMB_IQU-commander_1024_R2.02_full.fits',1024,60.0,1,2,3,1e-6],
['commander2015/1024_60.00smoothed_COM_CompMap_DustPol-commander_1024_R2.00.fits',1024,60.0,1,0,4,1.0],
['commander2015/1024_60.00smoothed_COM_CompMap_DustPol-commander_1024_R2.00.fits',1024,60.0,1,1,5,1.0],
['commander2015/COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,3,6,1.0],
['commander2015/COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,6,7,1.0]]

freqbands_I = [[0.408,0.408,'Haslam','grey', 0.5,'-'],
[0.820,0.820,'Dwingeloo','grey', 0.5,'-.'],
[1.42,1.420,'Stockert','grey', 0.5,'--'],
[2.3,2.30,'HartRAO','grey', 0.5,':'],
[4.5,5.5,'C-BASS','g', 0.3],
[19.5,25.0,'WMAP','b', 0.1],
[28.0,37.0,'','b', 0.1],
[35,46.0,'','b', 0.1],
[53.0,69.0,'','b', 0.1],
[82.0,106.0,'','b', 0.1],
[30.0-(0.2*30.0)/2,30.0+(0.2*30.0)/2,'Planck','r', 0.1], # Planck
[44.0-(0.2*44.0)/2,44.0+(0.2*44.0)/2,'','r', 0.1], # Planck
[70.0-(0.2*70.0)/2,70.0+(0.2*70.0)/2,'','r', 0.1], # Planck
[100.0-(0.33*100.0)/2,100.0+(0.33*100.0)/2,'','r', 0.1], # Planck
[143.0-(0.33*143.0)/2,143.0+(0.33*143.0)/2,'','r', 0.1], # Planck
[217.0-(0.33*217.0)/2,217.0+(0.33*217.0)/2,'','r', 0.1], # Planck
[353.0-(0.33*353.0)/2,353.0+(0.33*353.0)/2,'','r', 0.1], # Planck
[545.0-(0.33*545.0)/2,545.0+(0.33*545.0)/2,'','r', 0.1], # Planck
[857.0-(0.33*857.0)/2,857.0+(0.33*857.0)/2,'','r', 0.1]] # Planck

freqbands_P = [[4.5,5.5,'C-BASS','g', 0.3],
[19.5,25.0,'WMAP','b', 0.1],
[28.0,37.0,'','b', 0.1],
[35,46.0,'','b', 0.1],
[53.0,69.0,'','b', 0.1],
[82.0,106.0,'','b', 0.1],
[30.0-(0.2*30.0)/2,30.0+(0.2*30.0)/2,'Planck','r', 0.1],
[44.0-(0.2*44.0)/2,44.0+(0.2*44.0)/2,'','r', 0.1],
[70.0-(0.2*70.0)/2,70.0+(0.2*70.0)/2,'','r', 0.1],
[100.0-(0.33*100.0)/2,100.0+(0.33*100.0)/2,'','r', 0.1],
[143.0-(0.33*143.0)/2,143.0+(0.33*143.0)/2,'','r', 0.1],
[217.0-(0.33*217.0)/2,217.0+(0.33*217.0)/2,'','r', 0.1],
[353.0-(0.33*353.0)/2,353.0+(0.33*353.0)/2,'','r', 0.1]]

# 0  = 'GAL020  '           / 20% sky coverage                               
# 1  = 'GAL040  '           / 40% sky coverage                               
# 2  = 'GAL060  '           / 60% sky coverage                               
# 3  = 'GAL070  '           / 70% sky coverage                               
# 4  = 'GAL080  '           / 80% sky coverage                               
# 5  = 'GAL090  '           / 90% sky coverage                               
# 6  = 'GAL097  '           / 97% sky coverage                               
# 7  = 'GAL099  '           / 99% sky coverage   
mask_min=[['commander2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits',2048,4]]
mask_max=[['commander2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits',2048,5]]
# mask_min=[['commander2015/commander_dx11d2_mask_temp_n2048_fullres_v3_n0256.fits',256,0]]
# mask_max=[['commander2015/commander_mask_n256_likelihood_v1.fits',256,0]]

frequencies = [4.76]#, 11.0, 13.0, 17.0, 19.0, 22.8]
names = ['cbass']#'commander_sim_q11', 'commander_sim_q13', 'commander_sim_q17', 'commander_sim_q19', 'commander_sim_q23']

for i in range(0,len(frequencies)):
	commander_repro_maps(outdir='commander_freq_maps/', name=names[i], maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',nside=256,res=60.0,freq=frequencies[i])

