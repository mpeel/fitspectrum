from plotspectrum import plotspectrum
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

freqbands_I = [[0.408,0.408,'Haslam','maroon'],
[0.820,0.820,'Dwingaloo','purple'],
[1.42,1.420,'Reich','mediumblue'],
[2.3,2.30,'Jonas','darkblue'],
[4.5,5.5,'C-BASS','g'],
[19.5,25.0,'WMAP bands','b'],
[28.0,37.0,'','b'],
[35,46.0,'','b'],
[53.0,69.0,'','b'],
[82.0,106.0,'','b'],
[30.0-(0.2*30.0)/2,30.0+(0.2*30.0)/2,'Planck bands','r'], # Planck
[44.0-(0.2*44.0)/2,44.0+(0.2*44.0)/2,'','r'], # Planck
[70.0-(0.2*70.0)/2,70.0+(0.2*70.0)/2,'','r'], # Planck
[100.0-(0.33*100.0)/2,100.0+(0.33*100.0)/2,'','r'], # Planck
[143.0-(0.33*143.0)/2,143.0+(0.33*143.0)/2,'','r'], # Planck
[217.0-(0.33*217.0)/2,217.0+(0.33*217.0)/2,'','r'], # Planck
[353.0-(0.33*353.0)/2,353.0+(0.33*353.0)/2,'','r'], # Planck
[545.0-(0.33*545.0)/2,545.0+(0.33*545.0)/2,'','r'], # Planck
[857.0-(0.33*857.0)/2,857.0+(0.33*857.0)/2,'','r']] # Planck

freqbands_P = [[4.5,5.5,'C-BASS','g'],
[19.5,25.0,'WMAP bands','b'],
[28.0,37.0,'','b'],
[35,46.0,'','b'],
[53.0,69.0,'','b'],
[82.0,106.0,'','b'],
[30.0-(0.2*30.0)/2,30.0+(0.2*30.0)/2,'Planck bands','r'],
[44.0-(0.2*44.0)/2,44.0+(0.2*44.0)/2,'','r'],
[70.0-(0.2*70.0)/2,70.0+(0.2*70.0)/2,'','r'],
[100.0-(0.33*100.0)/2,100.0+(0.33*100.0)/2,'','r'],
[143.0-(0.33*143.0)/2,143.0+(0.33*143.0)/2,'','r'],
[217.0-(0.33*217.0)/2,217.0+(0.33*217.0)/2,'','r'],
[353.0-(0.33*353.0)/2,353.0+(0.33*353.0)/2,'','r']]

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

plotspectrum(outdir='', name='cbass_fig1', maps=maps_I, mask_min=mask_min,mask_max=mask_max,spd_file='amemodels/spdust2_wim.dat',minfreq=0.2,maxfreq=2000.0,numxpoints=1000,ymin=1,ymax=5e7,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True,freqbands=freqbands_I,nside=256,res=60.0,pol=False)

plotspectrum(outdir='', name='cbass_fig1_reproPlanck', maps=maps_I, mask_min=mask_min,mask_max=mask_max,spd_file='amemodels/spdust2_wim.dat',minfreq=10,maxfreq=1300.0,numxpoints=1000,ymin=2e-2,ymax=8e2,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True,freqbands=freqbands_I,nside=256,res=60.0,pol=False,legend=False)

# mask_min=[['commander2015/wmap_polarization_analysis_mask_r8_9yr_v5_trim.fits',250,0]]

plotspectrum(outdir='', name='cbass_fig1_pol', maps=maps_P, mask_min=mask_min,mask_max=mask_max,spd_file='amemodels/spdust2_wim.dat',minfreq=0.2,maxfreq=2000.0,numxpoints=1000,ymin=0.1,ymax=5e6,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True,freqbands=freqbands_P,nside=256,res=60.0,pol=True)

plotspectrum(outdir='', name='cbass_fig1_pol_reproPlanck', maps=maps_P, mask_min=mask_min,mask_max=mask_max,spd_file='amemodels/spdust2_wim.dat',minfreq=10,maxfreq=1300.0,numxpoints=1000,ymin=2e-2,ymax=8e2,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True,freqbands=freqbands_P,nside=256,res=60.0,pol=True,legend=False)