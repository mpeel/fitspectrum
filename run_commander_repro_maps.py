from commander_repro_maps import commander_repro_maps

mapdir = '/Users/mpeel/Documents/maps/planck_commander2015/'
outdir = mapdir+'commander_freq_maps/'

maps_I = [[mapdir+'COM_CompMap_Synchrotron-commander_0256_R2.00.fits',256,60.0,1,0,0,1.0],
[mapdir+'COM_CompMap_freefree-commander_0256_R2.00.fits',256,60.0,1,0,1,1.0],
[mapdir+'COM_CompMap_freefree-commander_0256_R2.00.fits',256,60.0,1,3,2,1.0],
[mapdir+'COM_CompMap_CMB-commander_0256_R2.00.fits',256,60.0,1,0,3,1.0],
[mapdir+'COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,1,0,4,1.0],
[mapdir+'COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,1,3,5,1.0],
[mapdir+'COM_CompMap_AME-commander_0256_R2.00.fits',256,60.0,2,0,6,1.0],
[mapdir+'COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,0,7,1.0],
[mapdir+'COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,3,8,1.0],
[mapdir+'COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,6,9,1.0]]

maps_P = [[mapdir+'256_60.00smoothed_COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits',256,60.0,1,0,0,1.0],
[mapdir+'256_60.00smoothed_COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits',256,60.0,1,1,1,1.0],
[mapdir+'1024_60.00smoothed_COM_CMB_IQU-commander_1024_R2.02_full.fits',1024,60.0,1,1,2,1e-6],
[mapdir+'1024_60.00smoothed_COM_CMB_IQU-commander_1024_R2.02_full.fits',1024,60.0,1,2,3,1e-6],
[mapdir+'1024_60.00smoothed_COM_CompMap_DustPol-commander_1024_R2.00.fits',1024,60.0,1,0,4,1.0],
[mapdir+'1024_60.00smoothed_COM_CompMap_DustPol-commander_1024_R2.00.fits',1024,60.0,1,1,5,1.0],
[mapdir+'COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,3,6,1.0],
[mapdir+'COM_CompMap_dust-commander_0256_R2.00.fits',256,60.0,1,6,7,1.0]]


# commander_repro_maps(outdir='commander_freq_maps/', name='cbass', maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',nside=256,res=60.0,freq=30.0)
# commander_repro_maps(outdir='commander_freq_maps/', name='cbass', maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',nside=256,res=60.0,freq=4.7)
# commander_repro_maps(outdir='commander_freq_maps/', name='wmap', maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',nside=256,res=60.0,freq=22.8)
# commander_repro_maps(outdir=outdir, name='', maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',galprop_file=mapdir+'Synchrotron_template_GHz_extended.txt',nside=256,res=60.0,freq=28.4)
commander_repro_maps(outdir=outdir, name='', maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',galprop_file=mapdir+'Synchrotron_template_GHz_extended.txt',nside=256,res=60.0,freq=10.0)

# frequencies = [4.76]#, 11.0, 13.0, 17.0, 19.0, 22.8]
# names = ['cbass']#'commander_sim_q11', 'commander_sim_q13', 'commander_sim_q17', 'commander_sim_q19', 'commander_sim_q23']
#
# for i in range(0,len(frequencies)):
# 	commander_repro_maps(outdir='commander_freq_maps/', name=names[i], maps=maps_I, spd_file='amemodels/spdust2_cnm.dat',nside=256,res=60.0,freq=frequencies[i])
