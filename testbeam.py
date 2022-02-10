import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits

def get_hfi_beam(FITSfile, indexnum=0):
	fits.info(FITSfile) # print list of extensions found in FITSfile
	data, header = fits.getdata(FITSfile, 0, header=True) # read extension #10 (data and header)
	# data, header = fits.getdata(FITSfile, 'ABC', header=True) # read extension having EXTNAME='ABC' (data and header)
	print header # print header
	print data.names # print column names
	# pylab.plot( data.field(0).flatten() ) # plot 1st column of binary table
	newdata = np.zeros(len(data))
	for i in range(0,len(data)):
		newdata[i] = data[i][indexnum]
	return newdata


#beamtf_K = np.loadtxt('wmap_ampl_bl_K1_9yr_v5p1.txt',usecols=(1,))

# HFI 2015
inputfits = fits.open('HFI_RIMO_Beams-100pc_R2.00.fits',)
print inputfits[3].header
print inputfits[3].data[0][0]
print np.shape(inputfits[3].data[0][0])
# exit()
beam2015 = inputfits[3].data[0][0]

# HFI 2018
beam = get_hfi_beam('Bl_T_R3.01_fullsky_100x100.fits')
beam2 = get_hfi_beam('Bl_TEB_R3.01_fullsky_100x100.fits')
beam2E = get_hfi_beam('Bl_TEB_R3.01_fullsky_100x100.fits',indexnum=1)
beam2B = get_hfi_beam('Bl_TEB_R3.01_fullsky_100x100.fits',indexnum=2)

conv_windowfunction = hp.gauss_beam(np.radians(9.65/60.0),3*2048)

print '1 degree is l=' + str(180/1.0)
print '2 degree is l=' + str(180/2.0)
print '5 arcmin is l=' + str(180*(60/5.0))

# plt.xscale('log')
plt.yscale('log')
plt.xlabel('l')
plt.ylabel('B(l)')
plt.plot(conv_windowfunction,label='Gaussian window function')
plt.plot(beam,label='Beam transfer function')
plt.plot(beam2015,label='Beam transfer function (2015)')
plt.plot(beam2,label='Beam transfer function (2018 T)')
plt.plot(beam2E,label='Beam transfer function (2018 E)')
plt.plot(beam2B,label='Beam transfer function (2018 B)')
# plt.legend()
l = plt.legend(prop={'size':8})
l.set_zorder(20)
plt.savefig('test_plotwf_2018.pdf')
