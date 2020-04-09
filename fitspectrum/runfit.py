from fitspectrum import fitspectrum
# from spectra import *
# import numpy as np

# const = get_spectrum_constants()
# test1 = freefree(const, 100.0, 1.0, 4000.0, 1.0)
# test2 = freefree(const, 101.0, 1.0, 4000.0, 1.0)
# alpha = np.log(test1/test2) / np.log((100.0) / (101.0))
# print alpha

# fitspectrum('seds/ngc2403.txt',srcname="NGC2403",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nosyncbeta=True,nofreefree=False,noame=True,nocmb=True,nodust=False)
# fitspectrum('seds/m101.txt',srcname="M101",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=False,nocmb=True,nodust=False)
fitspectrum('seds/ngc6946_fake.txt',srcname="NGC 6946 Added AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=False,nocmb=True,nodust=False,nosyncbeta=True)
fitspectrum('seds/m51_fake.txt',srcname="M51 Added AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=False,nocmb=True,nodust=False,nosyncbeta=True)
fitspectrum('seds/ngc891_fake.txt',srcname="NGC 891 Added AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=False,nocmb=True,nodust=False,nosyncbeta=True)

fitspectrum('seds/ngc6946_fake.txt',srcname="NGC 6946 not fitting AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False,nosyncbeta=True)
fitspectrum('seds/m51_fake.txt',srcname="M51 not fitting AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False,nosyncbeta=True)
fitspectrum('seds/ngc891_fake.txt',srcname="NGC 891 not fitting AME",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False,nosyncbeta=True)
# exit()
fitspectrum('seds/ngc6946.txt',srcname="NGC 6946",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False)#,nosyncbeta=True)
fitspectrum('seds/ngc891.txt',srcname="NGC 891",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False)#,nosyncbeta=True)
fitspectrum('seds/m51.txt',srcname="M51",outdir='sedfits_srt/',format=1,maxfitfreq=2000,mcmc=False,nosync=False,nofreefree=False,noame=True,nocmb=True,nodust=False)#,nosyncbeta=True)
exit()


# fitspectrum('seds/m31_srt.txt',srcname="M31_srt",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)

fitspectrum('seds/m31_dx11d_cbassonly.txt',srcname="M31c_joint",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False)

fitspectrum('seds/m31_dx11d_quijotejoint_cbass.txt',srcname="M31qc_joint",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False)

# fitspectrum('seds/m31_dx11d_quijoterasterrescale_cbass.txt',srcname="M31qc_rescale",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False)
# fitspectrum('seds/m31_dx11d.txt',srcname="M31",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False)
exit()

fitspectrum('seds/m31_dx11d_quijote.txt',srcname="M31q",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)
fitspectrum('seds/m31_dx11d_cbass.txt',srcname="M31c",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)
fitspectrum('seds/m31_dx11d_quijote_cbass.txt',srcname="M31qc",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)
fitspectrum('seds/m31_dx11d_srt.txt',srcname="M31s",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)
fitspectrum('seds/m31_dx11d_quijote_srt.txt',srcname="M31qs",outdir='sedfits/',format=1,maxfitfreq=2000,mcmc=True,nosync=False,nofreefree=False,noame=False,nocmb=False,nodust=False,nodust2=True)

# fitspectrum('seds/ngc253_spectrum.txt',format=2)
# fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis",format=3,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,mcmc=False)
# fitspectrum('seds/l2pup_sed.dat',srcname="L2 Puppis Tdust2",format=3,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
# fitspectrum('seds/l2pup_sed_band7.dat',srcname="L2 Puppis Tdust2_band7",format=1,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
# fitspectrum('seds/l2pup_sed_noband7.dat',srcname="L2 Puppis Tdust2_noband7",format=1,nosync=True,nofreefree=True,noame=True,nodust=False,nodust2=False,maxfitfreq=3.0e5,fixdust2temp=3500,quiet=True,mcmc=False)
