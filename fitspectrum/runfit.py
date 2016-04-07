from fitspectrum import fitspectrum
# from spectra import *
# import numpy as np

# const = get_spectrum_constants()
# test1 = freefree(const, 100.0, 1.0, 4000.0, 1.0)
# test2 = freefree(const, 101.0, 1.0, 4000.0, 1.0)
# alpha = np.log(test1/test2) / np.log((100.0) / (101.0))
# print alpha

# fitspectrum('seds/ngc253_spectrum.txt')
fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis",format=2,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5)
# fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis Tdust2",format=2,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True)
