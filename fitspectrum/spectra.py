import numpy as np
import scipy as sp

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
	thermaldust = 2.0 * const['h'] * (freq*1e9)**3/(2.997e8)**2 * (freq/optical_depth_freq)**index / (np.exp(xx)-1.0)
	S = amp * thermaldust * solid_angle *1e26
	return S

def paraobla (freq, coeff1, coeff2, coeff3):
	S = np.exp(coeff1 + coeff2*np.log(freq) + coeff3*np.log(freq)**2 )
	return S
 