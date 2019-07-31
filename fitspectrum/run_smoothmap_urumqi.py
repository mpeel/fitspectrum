from smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

# indir = '/Users/mpeel/Documents/maps/urumqi5ghz/'

# smoothmap(indir,indir,'urumqi_hp.fits','512_60.00smoothed_urumqi_4.8_512_mKCMBunits.fits', fwhm_arcmin=np.sqrt(60.0**2-9.5**2),nside_out=512,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2])

indir = '/Users/mpeel/Documents/maps/cbass2019/'

smoothmap(indir,indir,'cbass_global8p8deg_swapQU_NIGHT_v28allelsNs_37_noiseCut_masked5pc_G_1024_ol500_lessTol_g_map_g_1deg_0256.fits','512_60.00smoothed_cbass_4.76_512_mKCMBunits.fits', fwhm_arcmin=np.sqrt(60.0**2-45**2),nside_out=512,usehealpixfits=True,minmapvalue=-1e10,maxmapvalue=1e10,minmaxmaps=[0,1,2])
