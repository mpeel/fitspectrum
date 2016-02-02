import numpy as np
from math import pi

def synchrotron(freq, basefreq, sync_amp, sync_spec):
	return sync_amp * (freq/basefreq)**sync_spec

def freefree(freq, EM, Te):

	# previous equations, from old/standard text books!
	# g_ff = alog(4.955d-2 / freq) + 1.5*alog(Te)
	# tau_ff = 3.014d-2 * Te^(-1.5) * freq^(-2d) * EM * g_ff
	# T_ff = Te * (1d - exp(-tau_ff))

	# new equations from Draine (2011) book!
	g_ff = alog(exp(5.960 - (np.sqrt(3.0)/pi)*alog(1.0*freq*(Te/10000.0)^(-3.0/2.0))) + 2.71828)
	tau_ff = 5.468-2 * Te^(-1.5) * freq^(-2.0) * EM * g_ff
	# tau_ff = 1.772d-2 * (Te)^(-1.5) * (freq)^(-2d) * EM * g_ff

	# old Altenhoff approximation, also in many text books including
	# Scheffler book
	# tau_ff = 8.235d-2 * Te^(-1.35) * freq^(-2.1) * EM

	# finally get brightness temperature from the optical depthh
	T_ff = Te * (1.0 - exp(-tau_ff))

	# deal with very small tau values in the exponential
	#smalltau = where(tau_ff < 1.0e-10)
	#if (smalltau[0] GE 0) THEN T_ff[smalltau] = Te * (tau_ff[smalltau] - (-tau_ff[smalltau])^2.0/ 2.0 - (-tau_ff[smalltau]^3.0)/6.0)

	S = 2.0 * 1381.0 * T_ff * (freq*1e9)^2 / c^2  * solid_angle

	return S
