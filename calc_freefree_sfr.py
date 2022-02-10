import numpy as np

def sfr_murphy(Te, nu_ghz, flux_jy, m):
	return 4.6e-28 * np.power(Te/1e4, -0.45) * np.power(nu_ghz,0.1) * (flux_jy * 4*3.141592 * 100*m * 100*m / 1e23)

def sfr_condon(nu_ghz, flux_jy, m):
	return (1.0/5.5e20) * np.power(nu_ghz, 0.1) * (flux_jy * 4*3.141592 * m * m / 1e26)

def nc_mezger(Te, nu_ghz, flux_jy, D_kpc):
	return 4.761e48 * np.power(nu_ghz, 0.1) * np.power(Te, -0.45) * flux_jy * np.power(D_kpc, 2.0)

def sfr_kennicutt(Te, nu_ghz, flux_jy, D_kpc):
	return 7.5e-54 * nc_mezger(Te, nu_ghz, flux_jy, D_kpc)

distance = 0.79 # MPc
m = distance * 1000000 * 3.08e16 # Distance in metres
flux_mike = 0.7
freq_mike = 10.1
flux_srt = 0.35
freq_srt = 1.0
Te = 8000
print(sfr_murphy(Te, freq_srt, flux_srt, m))
print(sfr_murphy(Te, freq_mike, flux_mike, m))

print(sfr_condon(freq_srt, flux_srt,m))
print(sfr_condon(freq_mike, flux_mike,m))

print(sfr_kennicutt(Te, freq_srt, flux_srt,distance * 1000))
print(sfr_kennicutt(Te, freq_mike, flux_mike,distance * 1000))
print(sfr_kennicutt(Te, freq_srt, flux_srt*2.0,distance * 1000))
