PRO write_cvbl_planck3, fwhm, nsin, nsout, btag, annot
;
;  Finds beam window function amplitude CVBL to convolve the variance
;  maps. The variance beam is the square of the temperature beam, so
;  first we transform the temperature beam window function CONVBL to
;  get the beam profile, then square it and transform back to get CVBL.

; Input:
;        fwhm:  Required fwhm beamwidth (degrees)
;        nsin:  Healpix N_Side of input map
;        nsout: N_side of output map
;        btag:  Last part of input bl file name
;        annot: String, added to end of main part of file name. (eg."v1").
;
; Version 0.2 18 Oct 2013: add pixel window functions (in & out)
; Version 0.3  4 Nov 2013: try truncating before transforming
; Version 3   23 Mar 2014: Update to read new bl fits file format.
; Version 3.1 23 Mar 2014: Added annot, tidied a bit.
; Version 3.2  3 Jul 2014: Added btag input (saves hacking code each time!)
;
code   = ['30','44','70']
prog = 'write_cvbl_planck'
version = '3.1'

IF N_ELEMENTS(annot) NE 1 THEN annot = ''
IF N_ELEMENTS(btag) NE 1 THEN btag = '_T'
FOR i =0,2 DO BEGIN
                                ; Open beam window function file
   beamdir = 'bls/'
   fname = beamdir+'bls_GB_0'+code[i]+btag+'.fits'

   fits2cl, bl, fname
   nbl = n_elements(bl) - 1
   ll = lindgen(nbl + 1)

   pixw_in = HEALPIXWINDOW(nsin)
   sin = SIZE(pixw_in)
   nbli = sin[1]-1
   pixw_out = HEALPIXWINDOW(nsout)
   sout = SIZE(pixw_out)
   nlmax = 3*nsout-1
   nblo = sout[1]-1
   nbl = MIN([nbl,nbli,nblo, nlmax])
; Generate normalised radial profile of convolving beam:
   convbl = bl[0] * gaussbeam(fwhm*60.,nbl) * pixw_out[0:nbl,0] $
            / (bl[0:nbl] * pixw_in[0:nbl,0])

; Fudge to kill divide-by zeroes from inadequate-precision bl values:
   mask = WHERE(bl EQ 0.0,nz)
   IF mask[0] NE -1 THEN BEGIN
      PRINT, 'Suppressing', nz,'nulls in B_l'
      convbl[mask] = 0.0
   ENDIF

; Make Lagrange-polynomial transform of amp to get beam window
; function

; Choose scale to match size of beam. First find half-power point
   junk = WHERE(convbl GT 0.5, lcount)
   lhalf = junk[lcount-1]

; Calculate beam out to about 40 * half-power radius, roughly (note that
; this gives 100 points out to the half-power point).

   rad = FINDGEN(4001)*10.0*!pi/(lhalf*4000.)

   x = COS(rad)
   sinrad = SIN(rad)
   lgndr = FLTARR(4001,nbl+1)

   FOR l= 0,nbl DO lgndr[*,l] = LEGENDRE(x,l)

; Generate radial profile of convolving beam:
   conva = FLTARR(4001)
   FOR j=0,4000 DO conva[j] = TOTAL((ll+0.5)*convbl*lgndr[j,*],/DOUBLE)

   conva = conva / (2.0*!pi)

   PRINT, 'Peak of convolving beam is ', conva[0], MAX(conva)

; Define variance beam amplitude array of same size as convbl
   cvbl = convbl

; Square convolving beam and convert back to window function
   mult = sinrad*conva^2
   FOR l = 0,nbl DO cvbl[l] = INT_TABULATED(rad,mult*lgndr[*,l])

; Put in 2pi normalization factor:
   cvbl = 2.0*!pi*cvbl

   PRINT, 'CVBL[0] =', cvbl[0]

; Write out beam window function
;    filout = code[i]+STRING(fwhm,FORMAT="('_c',F3.1,'_var_window.txt')")
   filout = '0'+code[i]+'_c'+strtrim(fwhm,2)+'_var_window_' + $
            STRTRIM(nsin,2)+'_'+STRTRIM(nsout,2)+annot+'.txt'
   PRINT, 'Writing ', filout
   OPENW, 2, beamdir+filout
   PRINTF, 2, '# Beam window function for convolving variance maps'
;    PRINTF, 2, '# Output beam: ',STRING(fwhm*60.,FORMAT="(F5.1)"),' arcmin'
   PRINTF, 2, '# Output beam: ',STRtrim(fwhm*60.,2)+' arcmin'
   PRINTF, 2, '# Created by '+prog+' v'+version+' at ', SYSTIME()
   FOR l=0,nbl DO PRINTF, 2, l, cvbl[l]
   
   CLOSE, 2
    
ENDFOR

QUIT:

END
