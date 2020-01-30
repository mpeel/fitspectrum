import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
basedir = '/Users/mpeel/Desktop/Low-Frequency_Survey2/masks/'

cbass = hp.read_map(basedir+'mask_cbass_n0064.fits')
hp.mollview(cbass)
plt.savefig(basedir+'mask_cbass_n0064.png')
plt.clf()


chipass = hp.read_map(basedir+'mask_chipass_n0064.fits')
hp.mollview(chipass)
plt.savefig(basedir+'mask_chipass_n0064.png')
plt.clf()

spass = hp.read_map(basedir+'mask_spass_n0064.fits')
hp.mollview(spass)
plt.savefig(basedir+'mask_spass_n0064.png')
plt.clf()

numpix = len(cbass)

haslam = np.ones(numpix)
try:
	hp.write_map(basedir+'mask_haslam_n0064.fits',haslam)
except:
	null = 0
hp.mollview(haslam)
plt.savefig(basedir+'mask_haslam_n0064.png')
plt.clf()

cbass_chipass = cbass+chipass
cbass_chipass[cbass_chipass > 1] = 1
try:
	hp.write_map(basedir+'mask_cbass_chipass_n0064.fits',cbass_chipass)
except:
	null = 0
hp.mollview(cbass_chipass)
plt.savefig(basedir+'mask_cbass_chipass_n0064.png')
plt.clf()

cbass_spass = cbass+spass
cbass_spass[cbass_spass > 1] = 1
try:
	hp.write_map(basedir+'mask_cbass_spass_n0064.fits',cbass_spass)
except:
	null = 0
hp.mollview(cbass_spass)
plt.savefig(basedir+'mask_cbass_spass_n0064.png')
plt.clf()


haslam_cbass = haslam+cbass
haslam_cbass[haslam_cbass > 1] = 1
try:
	hp.write_map(basedir+'mask_haslam_cbass_n0064.fits',haslam_cbass)
except:
	null = 0
hp.mollview(haslam_cbass)
plt.savefig(basedir+'mask_haslam_cbass_n0064.png')
plt.clf()


haslam_chipass = haslam+chipass
haslam_chipass[haslam_chipass > 1] = 1
try:
	hp.write_map(basedir+'mask_haslam_chipass_n0064.fits',haslam_chipass)
except:
	null = 0
hp.mollview(haslam_chipass)
plt.savefig(basedir+'mask_haslam_chipass_n0064.png')
plt.clf()


haslam_spass = haslam+spass
haslam_spass[haslam_spass > 1] = 1
try:
	hp.write_map(basedir+'mask_haslam_spass_n0064.fits',haslam_spass)
except:
	null = 0
hp.mollview(haslam_spass)
plt.savefig(basedir+'mask_haslam_spass_n0064.png')
plt.clf()


spass_chipass = spass+chipass
spass_chipass[spass_chipass > 1] = 1
try:
	hp.write_map(basedir+'mask_spass_chipass_n0064.fits',spass_chipass)
except:
	null = 0
hp.mollview(spass_chipass)
plt.savefig(basedir+'mask_spass_chipass_n0064.png')
plt.clf()
