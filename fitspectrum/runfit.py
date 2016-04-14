from fitspectrum import fitspectrum
# from spectra import *
# import numpy as np

# const = get_spectrum_constants()
# test1 = freefree(const, 100.0, 1.0, 4000.0, 1.0)
# test2 = freefree(const, 101.0, 1.0, 4000.0, 1.0)
# alpha = np.log(test1/test2) / np.log((100.0) / (101.0))
# print alpha

# fitspectrum('seds/m31_dx11d.txt',srcname="M31",format=1,maxfitfreq=2000)

# fitspectrum('seds/ngc253_spectrum.txt',format=2)
# fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis",format=3,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,mcmc=False)
# fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis Tdust2",format=3,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
fitspectrum('seds/l2pup_sed_band7.dat',srcname="L2 Puppis Tdust2_band7",format=1,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
fitspectrum('seds/l2pup_sed_noband7.dat',srcname="L2 Puppis Tdust2_noband7",format=1,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
