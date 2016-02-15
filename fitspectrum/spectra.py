import numpy as np
import scipy as sp
from math import pi

def get_spectrum_constants():
	const = {
		'h': 6.626e-34,
		'k': 1.381e-23,
		'c': 2.997e8,
		'pi': pi,
		'dust_optical_depth_freq': 1198.8,
		'tcmb': 2.7255
	}
	return const

def synchrotron(const, freq, basefreq, sync_amp, sync_spec):
	return sync_amp * (freq/basefreq)**sync_spec

def freefree(const, freq, EM, Te, solid_angle, equation=1):

	if equation == 1:
		# new equations from Draine (2011) book!
		g_ff = np.log(np.exp(5.960 - (np.sqrt(3.0)/const['pi'])*np.log(1.0*freq*(Te/10000.0)**(-3.0/2.0))) + 2.71828)
		tau_ff = 5.468e-2 * Te**(-1.5) * freq**(-2.0) * EM * g_ff
		# tau_ff = 1.772d-2 * (Te)^(-1.5) * (freq)^(-2d) * EM * g_ff
	elif equation == 2:
		# previous equations, from old/standard text books!
		g_ff = np.log(4.955e-2 / freq) + 1.5*np.log(Te)
		tau_ff = 3.014e-2 * Te**(-1.5) * freq**(-2.0) * EM * g_ff
	elif equation == 3:
		# old Altenhoff approximation, also in many text books including
		# Scheffler book
		tau_ff = 8.235e-2 * Te**(-1.35) * freq**(-2.1) * EM

	# finally get brightness temperature from the optical depth
	T_ff = Te * (1.0 - np.exp(-tau_ff))

	# deal with very small tau values in the exponential
	#smalltau = where(tau_ff < 1.0e-10)
	#if (smalltau[0] GE 0) THEN T_ff[smalltau] = Te * (tau_ff[smalltau] - (-tau_ff[smalltau])^2.0/ 2.0 - (-tau_ff[smalltau]^3.0)/6.0)

	S = 2.0 * 1381.0 * T_ff * (freq*1e9)**2 / const['c']**2  * solid_angle

	return S

def spinningdust(spd_freq, spd_em, solid_angle, freq, amp, shift):

	# shift and interpolate
	interp = sp.interpolate.interp1d(spd_freq, spd_em)
	spinningdust = interp(freq)

	# ; if outside range of spinning dust model, then set to 0
	# IF (spdustshift NE 0.) THEN BEGIN
	# badvals = where(nu LT min(spinningdust_freq)*2. OR nu GT (max(spinningdust_freq)-100.), nbadvals)
	# IF (nbadvals GT 0) THEN spinningdust[badvals] = 0.
	# ENDIF

	S = amp * spinningdust * solid_angle * 1.0e20
	return S

	# ; if outside range of spinning dust model, then set to 0
	# IF (spdustshift NE 0) THEN BEGIN
	# badvals = where(S LT 0.,nbadvals)
	# IF (nbadvals GT 0) THEN S[badvals]=0.
	# ENDIF


def cmb(const, freq, amp, solid_angle):

	xx = const['h'] * freq*1.0e9 / (const['k'] * const['tcmb'])
	xx2 = const['h'] * freq*1.0e9 / (const['k'] * (const['tcmb']+deltaT))
	S = (2.0e26 * const['h'] * (freq*1e9)**3 * solid_angle  / const['c']**2) * ((1.0/(np.exp(xx2)-1.0)) - (1.0/(np.exp(xx)-1.0)))
	return S

def thermaldust(const, freq, amp, index, temperature, optical_depth_freq, solid_angle):

	xx = const['h'] * freq*1.0e9 / (const['k'] * temperature)

	# fit optical depth at 250 um (1198.8 GHz) - dust_optical_depth_freq
	thermaldust = 2.0 * const['h'] * (freq*1e9)**3/(const['c'])**2 * (freq/optical_depth_freq)**index / (np.exp(xx)-1.0)
	S = amp * thermaldust * solid_angle *1e26
	return S

def paraobla(freq, coeff1, coeff2, coeff3):
	S = np.exp(coeff1 + coeff2*np.log(freq) + coeff3*np.log(freq)**2 )
	return S
 
def planckcorr(const, nu_ghz):
	nu = nu_ghz * 1.0e9
	x = const['h'] * nu / (const['k'] * const['tcmb'])
	value = (np.exp(x)-1.0)**2.0 / (x**2. * np.exp(x))
	return value

# Function to calculate the conversion factors between different units.
# 
# Mike Peel, 1 October 2014 - Forked from haperflux.pro
# Mike Peel, 15 February 2016 - converted from IDL to python
def convertunits(const, units_in, units_out, frequency, pix_area=1.0):
	unitslist = ['K','mK','uK','K_RJ','mK_RJ','uK_RJ','K_CMB','mK_CMB','uK_CMB', 'MJy/sr', 'Jy/pix' ]
	# get conversion from the input units to Jy/pix
	if (units_in == 'K' or units_in == 'K_RJ' or units_in == 'KRJ'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area
	elif (units_in == 'mK' or units_in == 'mK_RJ' or units_in == 'mKRJ'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e3
	elif (units_in == 'uK' or units_in == 'uK_RJ' or units_in == 'uKRJ'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e6
	elif (units_in == 'K_CMB' or units_in == 'KCMB'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / planckcorr(const, frequency)
	elif (units_in == 'mK_CMB' or units_in == 'mKCMB'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e3 / planckcorr(const, frequency)
	elif (units_in == 'uK_CMB' or units_in == 'uKCMB'):
		factor_in = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e6 / planckcorr(const, frequency)
	elif (units_in == 'MJy/sr' or units_in == "MJY/SR" or units_in == "MjySr"):
		factor_in = pix_area * 1.0e6
	elif (units_in == 'Jy/pixel' or units_in == 'JY/PIXEL' or units_in == 'JY/PIX' or units_in == 'JyPix' or units_in == 'Jy/Pix'):
		factor_in = 1.0
	else:
		print 'Invalid unit conversion specified for convertunits (units_out='+units_out+'). ***Returning 1.*** Please use one of the following available units:'
		print unitlist
		return 1.0

	if (units_out == 'K' or units_out == 'K_RJ' or units_out == 'KRJ'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area
	elif (units_out == 'mK' or units_out == 'mK_RJ' or units_out == 'mKRJ'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e3
	elif (units_out == 'uK' or units_out == 'uK_RJ' or units_out == 'uKRJ'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e6
	elif (units_out == 'K_CMB' or units_out == 'KCMB'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / planckcorr(const, frequency)
	elif (units_out == 'mK_CMB' or units_out == 'mKCMB'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e3 / planckcorr(const, frequency)
	elif (units_out == 'uK_CMB' or units_out == 'uKCMB'):
		factor_out = 2.*1381.*(frequency*1.0e9)**2/(const['c'])**2 * pix_area / 1.0e6 / planckcorr(const, frequency)
	elif (units_out == 'MJy/sr' or units_out == "MJY/SR" or units_out == "MjySr"):
		factor_out = pix_area * 1.0e6
	elif (units_out == 'Jy/pixel' or units_out == 'JY/PIXEL' or units_out == 'JY/PIX' or units_out == 'JyPix' or units_out == 'Jy/Pix'):
		factor_out = 1.0
	else:
		print 'Invalid unit conversion specified for convertunits (units_out='+units_out+'). ***Returning 1.*** Please use one of the following available units:'
		print unitlist
		return 1.0

	# Return the ratio of the factors
	return factor_in/factor_out

