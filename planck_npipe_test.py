import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

indirectory='/Users/mpeel/Documents/maps/planck_npipe/Single-frequency/'
filenames = ['HFI_SkyMap_100_2048_R4.00_full.fits','HFI_SkyMap_143_2048_R4.00_full.fits','HFI_SkyMap_217_2048_R4.00_full.fits','HFI_SkyMap_353_2048_R4.00_full.fits','HFI_SkyMap_545_2048_R4.00_full.fits','HFI_SkyMap_857_2048_R4.00_full.fits','LFI_SkyMap_030_1024_R4.00_full.fits','LFI_SkyMap_044_1024_R4.00_full.fits','LFI_SkyMap_070_1024_R4.00_full.fits']
for filename in filenames:
	data,hdr = hp.read_map(indirectory+filename,field=None,h=True)
	print(hdr)
	hp.mollview(data[0],max=np.max(data[0])/1000.0)
	plt.savefig(indirectory+filename+'.png')
	plt.clf()
	hp.mollview(np.sqrt(data[1]**2+data[2]**2),max=np.max(data[0])/1000.0)
	plt.savefig(indirectory+filename+'_pol.png')
	plt.clf()
