;Planckcorr function for cmb conversion
FUNCTION planckcorr, nu_ghz

  k_b = 1.3806E-23              ; J/K
  h = 6.6262E-34                ; J*s
  T_cmb = 2.73                  ; K	
  
  nu = nu_ghz*1.0E9             ; Hz
  x = h * nu / (k_b * T_cmb)
  
  result = (exp(x)-1.)^2. / (x^2. * exp(x)) ; thermodynamic uK  

  return, result
END

PRO new_template_fitting,templatelist,regionlist,nside,savename,data=data,montecarlo=montecarlo,$
                         domask=domask,cmbsub=cmbsub,dographs=dographs,silent=silent
  
;Pro for template fitting both wmap and plank maps using both cmb
;covariance matrix as well as cmb subtraction depending on choice made
  
;templatelist -> text file containing a list of file names of templates to be fitted
;to the data. File names should be of the format:
;512_60.00smoothed_wmap9dec_22.8_mK  

;regionlist -> text file containg list of regions to be used in
;fitting process. This should be formatted as:
;region number,min longitude,max longitude,min lattitude,max lattitude

;nside -> The nside value to work at.
  
;savename -> The name for that all resulting products will be named
;followed by an appropriate add on ie:data_coeff for the coeficients

;data -> keyword to be set to a text file containing a list of file
;names of data to be fit e.g data='data.txt' (same format as template
;names within file). THIS PRO PERFORMS TEMPLATE FITTING FOR WMAP AND PLANCK

;montecarlo -> keyword to be set for the number of montecarlo
;simulations to be performed when testing the code. In this case the
;simulatio coefficients wanted must be input as a second column in the
;template list text file.

;domask -> keyword to be set to perform masking

;cmbsub -> keyword to be set to perform subtraction of the CMB from
;the data instead of using the CMB covariance matrix which is used by
;default

;silent -> Keyword set to reduce output to simple percentage of completion


;---------------------------
;Check Keywords
;---------------------------
run=0
;Check if montecarlo set other keywords that aren't used cannot be set
IF keyword_set(montecarlo) THEN BEGIN
   IF keyword_set(data) THEN BEGIN 
      print,'ERROR: Cannot set both data and montecarlo keywords at the same time!'
      run=1
      ENDIF ELSE IF keyword_set(domask) THEN BEGIN
      print,'ERROR: Cannot set both domask and montecarlo keywords at the same time!'
      run=1
      ENDIF ELSE IF keyword_set(cmbsub) THEN BEGIN
      print,'ERROR: Cannot set both cmbsub and montecarlo keywords at the same time!'
      run=1
      ENDIF ELSE IF keyword_set(dographs) THEN BEGIN
      print,'ERROR: Cannot set both dographs and montecarlo keywords at the same time!'
      run=1
   ENDIF 
ENDIF  
;Check either montecarlo or data set
IF ~keyword_set(data) AND ~keyword_set(montecarlo) THEN BEGIN
   print,'ERROR: You must set either the montecarlo or data keywords!'
   run=1
ENDIF
;Check number of parameters
IF n_params() NE 4 THEN BEGIN
   print,'ERROR: Incorrect number of parameters!'
   run=1
ENDIF
;Don't run and print syntax if any of the previous keyword errors occur
IF run EQ 1 THEN BEGIN
   print,''
   print,'SYNTAX: PRO joe_template_fitting,templatelist,regionlist,nside,savename [data=data,montecarlo=montecarlo,/domask,/domaps,/cmbsub,/dographs]'
   print,''
   stop
ENDIF 
;Set number of simulations to 1 if fitting data so code only runs once
IF keyword_set(montecarlo) THEN numsims = long(montecarlo) ELSE IF keyword_set(data) THEN numsims = 1L

;Set nside to long
nside = long(nside)
;---------------------------
;Perform setup
;---------------------------
;Set up for the save directory/base directories
baseDir = '/mirror/data/jKaemena/semester2/template_fitting/source/'
resultsDir = '/mirror/data/jKaemena/semester2/template_fitting/results/'
graphDir = '/mirror/data/jKaemena/semester2/template_fitting/graphs/'

;Load ct if needed
IF keyword_set(dographs) THEN loadct,39

;Print what file is being used for templates
IF ~keyword_set(silent) THEN print, 'Reading in templates from: '+templatelist $
ELSE print,string(13b),'Template Fitting: 0%',format='(A,A,$)'


;Read in templates from text file
IF keyword_set(montecarlo) THEN readcol,templatelist,templatenames,coeffs,format='A,D',comment='#',/silent $
ELSE readcol,templatelist,templatenames,format='A',comment='#',/silent
;Get number of templates
ntemplates=size(templatenames)
ntemplates=ntemplates[1]

;Make two dimensional array to hold maps once downgraded to desired
;nside
npix = (nside^2)*12L
maps = dblarr(ntemplates,npix)
;Make array for split template name strings
mapsinfo = strarr(ntemplates,5)

;Read in templates fits files and downgrade
FOR i=0L, ntemplates-1 DO BEGIN
   ;Split fits file name string
   mapsinfo[i,*]=templatenames[i].split('_')
   ;Read map in
   read_fits_map,baseDir+templatenames[i],map,ordering=ordering
   ;Downgrade and store map in maps array 
   ud_grade,map,map,nside=nside,order_in=ordering,order_out='ring'
   ;Store map in maps array
   ;maps[i,*]=map[*]
   maps[i,*] = map[*,0]
   ;Convert to uK to ensure inversions are accurate (leave rydberg as is)
   IF strmatch(mapsinfo[i,4],'*mK*') EQ 1 THEN maps[i,*] = maps[i,*]*1.e3 $
   ELSE IF strmatch(mapsinfo[i,4],'*K*') EQ 1 THEN maps[i,*] = maps[i,*]*1.e6 ;$
   ;ELSE IF strmatch(mapsinfo[i,6],'*R*') EQ 1 THEN BEGIN
   ;   maps[i,*] = maps[i,*]*(13.60569*1.602e-19)
   ;   maps[i,*] = maps[i,*]/1.38e-23
   ;   maps[i,*] = maps[i,*] * 1.e6
   ;ENDIF   
ENDFOR
;Check if fitting data from text file
IF keyword_set(data) THEN BEGIN
   ;Read in data maps form text file
   IF ~keyword_set(silent) THEN print,'Reading in data from: '+data
   readcol,data,datanames,format='A',comment='#',/silent
   ;Get number of data
   ndata=size(datanames)
   ndata=ndata[1]
   ;Split data names
   datainfo=strarr(ndata,5)
   FOR i=0L,ndata-1 DO datainfo[i,*]=datanames[i].split('_')
   ;Open text files to print out results to
   OPENW,unit1,resultsDir+savename+'_coeff',/get_lun
   OPENW,unit2,resultsDir+savename+'_spectral',/get_lun
   ;Print out header for the results files
   printf,unit1,'#Coefficient results for file:',savename
   printf,unit1,''
   printf,unit2,'#Spectral Index results for file:',savename
   printf,unit2,''
ENDIF ELSE ndata=1;Set ndata=1 if not fitting data so when looping through data in montecarlo mode loop only perfromed once

;Read in regions to be used
IF ~keyword_set(silent) THEN print,'Reading in regions from: '+regionlist
readcol,regionlist,regionnumber,regionminlong,regionmaxlong,regionminlat,regionmaxlat,/silent,comment='#'
;Get number of regions
nregions = size(regionnumber)
nregions=nregions[1]

;----------------------------------
;LOOP OVER ALL REGIONS STARTS HERE
;----------------------------------
   
FOR regioncount=0L,nregions-1 DO BEGIN
;Set region for this iteration of loop
region=dblarr(7)
region[0]=regionnumber[regioncount]
region[1]=regionminlong[regioncount]
region[2]=regionmaxlong[regioncount]
region[3]=regionminlat[regioncount]
region[4]=regionmaxlat[regioncount]

;Print divider and info on loop status
IF ~keyword_set(silent) THEN BEGIN
print,'-------------------------------------'
print,'Performing fitting for region: '+string(region[0])
print,''
ENDIF
;Create arrays for graphs if needed (one graph per region)
IF keyword_set(dographs) THEN BEGIN
graph_results = dblarr(ndata,ntemplates)
graph_errors = dblarr(ndata,ntemplates)
graph_freqs = dblarr(ndata);All templates fitted to same data frequency
ENDIF

;-------------------------------------
;LOOP OVER NUMBER OF DATA STARTS HERE
;-------------------------------------

FOR datacount=0L,ndata-1 DO BEGIN

;Read in data if required
IF keyword_set(data) THEN BEGIN
   IF ~keyword_set(silent) THEN print,'Fitting: '+datanames[datacount]
   ;Read in map
   read_fits_map,baseDir+datanames[datacount],y_map,ordering=ordering
   ;Downgrade to desired nside
   ud_grade,y_map,y_map,order_in=ordering,order_out='ring',nside=nside
   ;Convert to uK 
   IF strmatch(datainfo[datacount,4],'*mK*') EQ 1 THEN y_map=y_map*1.e3 $
   ELSE IF strmatch(datainfo[datacount,4],'*K*') EQ 1 THEN y_map=y_map*1.e6 
   ;Store frequency if needed for graphs
   IF keyword_set(dographs) THEN graph_freqs[datacount] = double(datainfo[datacount,3])
ENDIF

;----------------
;CMB subtraction
;----------------
IF keyword_set(cmbsub) THEN BEGIN
   IF ~keyword_set(silent) THEN print,'--->Subtracting cmb...'
   ;Read in cmb map
   IF strmatch(datainfo[datacount,2],'*wmap*') EQ 1 THEN read_fits_map,basedir+'wmap_cmb_9yr_smooth.fits',cmb,ordering=ordering $
   ELSE IF strmatch(datainfo[datacount,2],'*planck*') EQ 1 THEN read_fits_map,basedir+'COM_CompMap_CMB-smica_2048_R1.11.fits',cmb,ordering=ordering
   ud_grade,cmb,cmb,nside=nside,order_in=ordering,order_out='ring'
   ;Convert to uK
   cmb*=1.e3
   ;Subtract
   y_map[*,0]-=cmb
ENDIF

;--------
;Masking
;--------
;Check for masking and set up mask if needed (if no masking desired
;then mask is simply an array of 1's and hence has no effect on data)
;Standard mask (all values =1)
mask = dblarr(npix)
mask[*]=1.

IF keyword_set(domask) THEN BEGIN
   ;WMAP MASKING
   IF strmatch(datainfo[datacount,2],'*wmap*') EQ 1 THEN BEGIN
      IF ~keyword_set(silent) THEN print,'--->Using WMAP masking' 
      ;Read in masks 
      read_fits_map,basedir+'wmap_temperature_source_mask_r9_9yr_v5.fits',mask1,ordering=ordering
      ud_grade,mask1,mask1,nside=nside,order_in=ordering,order_out='ring'
      read_fits_map,basedir+'wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits',mask2,ordering=ordering
      ud_grade,mask2,mask2,nside=nside,order_in=ordering,order_out='ring'
      ;Combine masks
      FOR i=0L,npix-1 DO mask[i]=mask1[i]*mask2[i]
   ENDIF
   ;PLANCK MASKING
   IF strmatch(datainfo[datacount,2],'*planck*') EQ 1 THEN BEGIN
      IF ~keyword_set(silent) THEN print,'--->Using PLANCK masking'
      ;Read in masks
      read_fits_map,basedir+'LFI_Mask_PointSrc_2048_R2.00.fits',mask,ordering=ordering
      ud_grade,mask,mask,nside=nside,order_in=ordering,order_out='ring'
   ENDIF
ENDIF

;--------------------------
;Get relevant region of map
;--------------------------
;Get positions within selected area
pix2ang_ring, nside, [0:(npix-1)], thetarad, phirad

theta = thetarad * (180./!PI)
phi = phirad * (180./!PI)

longitude = phi
latitude = 90.-theta

;Find indicies of pixels in selected area (check for areas where
;region crosses from the <360 to >0 in longitude  
IF region[1] GT region[2] THEN indices = where(((longitude GT region[1]) OR (longitude LT region[2]))AND((latitude GT region[3])AND(latitude LT region[4])), nindices)$
ELSE indices = where(((longitude GT region[1])AND(longitude LT region[2]))AND((latitude GT region[3])AND(latitude LT region[4])), nvals)

;---------------------------------------
;Create box maps and temperature arrays
;---------------------------------------
;Create array for templates 
boxTemplates = dblarr(ntemplates,nvals)
FOR i=0L,ntemplates-1 DO boxTemplates[i,*]=maps[i,indices]*mask[indices]
;Angular positions
boxlat = latitude[indices]
boxlong = longitude[indices]
;Create arrays for fitting data if needed
IF keyword_set(data) THEN BEGIN
    boxData=y_map[indices,0]*mask[indices]
    ;Hits array for wmap or variance map for planck
    IF strmatch(datainfo[datacount,2],'*wmap*') EQ 1 THEN BEGIN 
       ;Define double values of nside_in and nside_out for wmap hits rescaling
       nside_in=double(datainfo[datacount,0])
       nside_out=double(nside)
       ;Input values for hits array
       boxHits=y_map[indices,1]*(nside_out/nside_in)
    ENDIF ELSE IF strmatch(datainfo[datacount,2],'*planck*') EQ 1 THEN boxvar=y_map[indices,1];THIS WILL NEED TO BE EDITED TO SUIT THE FORMAT OF PLANCK MAPS
ENDIF

;--------------
;Matrix Setup
;--------------
;Create Y array
;Produce linear combination simulation map if wanted
IF keyword_set(montecarlo) THEN BEGIN
   simY = dblarr(nvals)
   FOR i=0L,ntemplates-1 DO simY = simY+(coeffs[i]*boxTemplates[i,*])
  ENDIF ELSE IF keyword_set(data) THEN Y = boxData ;If fitting data           

;Input Values For X
X=dblarr(nvals,ntemplates)
FOR i=0L,ntemplates-1 DO X[*,i]=boxTemplates[i,*]
;Montecarlo arrays for results
IF keyword_set(montecarlo) THEN BEGIN
   a_results=dblarr(numsims,ntemplates)
   error_results=dblarr(numsims,ntemplates)
   chisq_results = dblarr(numsims)
ENDIF ELSE IF keyword_set(data) THEN BEGIN
   ;Create array to hold error results (sqrt(sigma))
   a_err=dblarr(ntemplates+1);+1 for offset
   ;Offset
   offset=dblarr(nvals)
   offset[*]=1
   ;Add offset 
   X=[[X],[offset]]
ENDIF
;Transpose X
transX=transpose(X)

;-----------------------
;Noise Covariance Matrix
;-----------------------
;Define the noise covariance matrix
noisecovar = identity(nvals,/double)
;Define noisescale to be used(mK)
noisescale=20.
;Set the covariance matrix up
IF keyword_set(data) THEN BEGIN
   IF (strmatch(datainfo[datacount,2],'*wmap*')EQ 1) THEN BEGIN
      ;Get appropriate noisescale for each wmap frequency(uK)
      IF datainfo[datacount,3] EQ 22.8 THEN wmapnoisescale=1424. $
      ELSE IF datainfo[datacount,3] EQ 33.0 THEN wmapnoiscale=1449. $
      ELSE IF datainfo[datacount,3] EQ 40.7 THEN wmapnoisescale=2267. $
      ELSE IF datainfo[datacount,3] EQ 60.7 THEN wmapnoisescale=3288. $
      ELSE IF datainfo[datacount,3] EQ 93.5 THEN wmapnoisescale=5852.
      ;WMAP noise scale
      FOR i=0L,nvals-1 DO noisecovar(i,i)=(wmapnoisescale^2)/boxHits[i] 
      ;Planck noise covariance matrix
   ENDIF ELSE IF strmatch(datainfo[datacount,2],'*planck*') EQ 1 THEN FOR i=0L,nvals-1 DO noisecovar(i,i)=boxvar[i]^0.5
ENDIF ELSE noisecovar = noisecovar * (noisescale^2);Simple noise covariance matrix (diagonal)

;---------------------
;CMB COVARIANCE MATRIX
;---------------------
IF ~keyword_set(cmbsub) THEN BEGIN
IF ~keyword_set(silent) THEN print,'--->Calculating CMB covariance matrix...'
;Define lmax to be used
lmax = 3*nside-1L

;Pixel Window Function
wpix=healpixwindow(nside)
;Beam function
B_l = gaussbeam(60.,lmax)

;Get CMB cl
;Text File version
readcol,basedir+'wmap_lcdm_pl_model_yr1_v1.txt',l,cmbpower,/silent
;Convert from power to cl
cl = 2.*!PI*cmbpower/(l*(l+1.))
;Convert to mK^2
;cl = cl/1.d6

;Set up vectors
cl=[0,0,cl];Insert first to values of cl as zero for l=0,1
l=[0,1,l];Insert l=0,1, both for text file version only
cl = cl[0:lmax]
l=l[0:lmax]
wpix=wpix[0:lmax]
B_l=B_l[0:lmax]

;Define CMB covariance matrix
cmbcovar = dblarr(nvals,nvals)
;Calculate matrix (via vectors)
FOR i=0L, nvals-1 DO BEGIN
   FOR j=0L, nvals-1 DO BEGIN
;Get Angle differences
      anglediff = sphdist(boxLong[i], boxLat[i], boxLong[j], boxLat[j],/degrees)
      anglediff = (anglediff*!PI)/180.
      cosine = cos(anglediff)
      ;Use vectors now
      pl = legendre(cosine, l, 0, /double)
      temp = (2*l+1.)* cl *pl* (wpix^2)*(B_l^2)      
      cmbcovar[i,j] = total(temp,/double)
   ENDFOR
ENDFOR
;Normalisation
cmbcovar = cmbcovar/(4.*!PI)

;Add noise and cmb together 
covar = noisecovar+cmbcovar
ENDIF ELSE BEGIN
;If using cmb subtraction add extra noise
FOR i=0L,nvals-1 DO noisecovar(i,i)+=10.^2
;Set covar
covar=noisecovar
ENDELSE

;--------------
;Invert Matrix
;--------------
IF ~keyword_set(silent) THEN print, '--->Inverting...'
invcovar = invert(covar,status,/double)
IF status NE 0 THEN BEGIN 
   print,''
   print,'*****WARNING COVARIANCE MATRIX INVERSION NOT ACCURATE*****'
   stop
ENDIF

;----------------------------
;MONTE CARLO LOOP STARTS HERE 
;----------------------------
;This is performed for both data and simulation but in the case of
;data only loops once for each data set
FOR q=0L,numsims-1 DO BEGIN
;Print status
IF keyword_set(montecarlo) AND ~keyword_set(silent) THEN print,string(13b),'--->Performing Montecarlo loop: ',q+1,'/',montecarlo,format='(A,A,I,A,I,$)'

;Perform any required functions for simulations
IF keyword_set(montecarlo) THEN BEGIN
   ;Create random noise 
   random=randomn(seed,nvals)
   noise = random*noisescale
   ;Add noise to data
   Y = simY+noise
   ;Get cl again to be safe
   readcol,basedir+'wmap_lcdm_pl_model_yr1_v1.txt',l,cmbpower2,/silent
   ;Convert from power to cl
   cl2 = 2.*!PI*cmbpower2/(l*(l+1L))
   ;Create CMB map to be added to simulated data
   myseed = randomn(seed,1,/long)
   isynfast,cl2,cmb,nside=nside,nlmax=lmax,fwhm_arcmin=60.,iseed=myseed,/silent
   ;Create temperature array of CMB for selected region
   boxCMB=cmb[indices,0]
   ;Add CMB to data
   Y = Y+boxCMB
ENDIF 

;--------------------
;Matrix Calculations
;--------------------
;Calculate sigma via invert
XTinvCX = transX # invcovar # X
sigma = invert(XTinvCX,/double,status)
IF status NE 0 THEN BEGIN
   print,''
   print,'*****WARNING MATRIX INVERSION NOT ACCURATE*****'
   stop
ENDIF
;Calculate a
a = sigma # transX # invcovar # Y

;Calculate chi squared
xA = a ## X
chisq = transpose(y-xA) # invcovar # (y-xA)

;Record montecarlo results
IF keyword_set(montecarlo) THEN BEGIN
chisq_results[q] = chisq
FOR i=0L,ntemplates-1 DO BEGIN
   a_results[q,i]=a[i]
   error_results[q,i]=sqrt(sigma[i,i])
ENDFOR ;Get errors from sigma when fitting data
ENDIF ELSE FOR i=0L,ntemplates DO a_err[i]=sqrt(sigma(i,i))

;END OF MONTE CARLO LOOP
ENDFOR

;Convert using planckorr from uKCMB/uK into uK/K or uK/R
IF keyword_set(data) THEN FOR i=0L,ntemplates-1 DO BEGIN 
   IF strmatch(mapsinfo[i,4],'*R*') EQ 1 THEN convert=1. ELSE convert=1.e6
   ;Planckorr conversion
   aa=double(datainfo[datacount,3])
   IF strmatch(datainfo[datacount,4],'*CMB*') EQ 1 THEN pc = planckcorr(aa) ELSE pc =1.;Only use conversion factor in CMB units case
   a[i]=convert*a[i]/pc
   a_err[i]=convert*a_err[i]/pc
ENDFOR

;Set reduced chi-squared
IF keyword_set(data) THEN red_chisq=chisq/(nvals-(ntemplates+1)) ELSE IF keyword_set(montecarlo) THEN red_chisq_results=chisq_results/(nvals-ntemplates)

;Store graph results if wanted
IF keyword_set(dographs) THEN BEGIN
   FOR i=0L,ntemplates-1 DO BEGIN
      graph_results[datacount,i]=a[i]
      graph_errors[datacount,i]=a_err[i]
   ENDFOR
ENDIF

;-------------------------
;Calculate spectral index
;-------------------------
IF keyword_set(data) THEN BEGIN
;Declare arrays
beta_results=dblarr(ntemplates,2);Two dimensional in order to hold info in whether or not it is an upper bound or not
beta_error=dblarr(ntemplates)

;Loop through results for each template
FOR i=0L, ntemplates-1 DO BEGIN
;Check for negative coefficient and add sigma until positive to account for it 
IF a[i] LT 0. THEN BEGIN
dummy=a[i]
;Count how many sigma are added in order to make value positive
sigma_count=0.
WHILE dummy LT 0. DO BEGIN
   ;Add one sigma
   dummy=dummy+a_err[i]
   ;Count +1
   sigma_count++
ENDWHILE
;Set second row to number of sigma added so know to set error to zero, indicating that the
;value is an upper bound
beta_results[i,1]=sigma_count
ENDIF ELSE dummy=a[i]

;Convert to correct units (at this point results in uK/K or uK/R)
IF (strmatch(mapsinfo[i,4],'*mK*') EQ 1) OR (strmatch(mapsinfo[i,4],'*K*') EQ 1) THEN dummy=dummy/1.e6 $
ELSE IF strmatch(mapsinfo[i,4],'*R*') EQ 1 THEN BEGIN
   dummy = dummy/(13.60569*1.602e-19)
   dummy = dummy*1.38e-23
   dummy = dummy / 1.e6
ENDIF    

;Get frequencies
datafreq=double(datainfo[datacount,3])
templatefreq=double(mapsinfo[i,3])
;Check if need to swap them if template > data
betaconvert=1.
IF templatefreq GT datafreq THEN betaconvert =-1.
;Calculate beta/error
beta_results[i,0]=betaconvert*alog(dummy)/alog(datafreq/templatefreq)
beta_error[i]=betaconvert*a_err[i]/(alog(datafreq/templatefreq)*a[i])
IF beta_results[i,1] NE 0. THEN beta_error[i]=0.

ENDFOR
ENDIF  

;-------------------------------------
;Histograms for montecarlo simulations
;-------------------------------------
IF keyword_set(montecarlo) THEN BEGIN 
IF ~keyword_set(silent) THEN BEGIN
print,''
print, '--->Calculating and plotting  histograms...'
print,''
print,'Results:'
ENDIF
;Template coefficient histograms
FOR i=0L,ntemplates-1 DO BEGIN
set_plot, 'ps'
loadct,39
;Set File Name
file = graphDir+savename+'_'+string(i+1,format='(3I0)')
DEVICE, file=file, /color, /landscape

plothist, a_results[*, i], xbin, hist, /autobin, xtitle = 'Coefficient', ytitle = 'Frequency',  title='Monte Carlo, n= ' + string(montecarlo, format='(d8.0)') + ', Expected value = ' + string(coeffs[i]), psym=1, charsize=1

;Gaussian Fit
yfit1 = gaussfit(xbin, hist, gausscoeff, nterms=3)

xarray = arrgen((gausscoeff[1]-10*gausscoeff[2]), (gausscoeff[1]+10*gausscoeff[2]), 0.01*gausscoeff[2])
yarray = arrgen(0, montecarlo, 0.1)
zsigma1 = (xarray-gausscoeff[1])/gausscoeff[2]
ymodel = gausscoeff[0] * exp(-(zsigma1^2)/2)
centreline = yarray
centreline[*] = gausscoeff[1]

oplot, xarray, ymodel, color=240
oplot, centreline, yarray, color=50, linestyle=2
notesxcoord= gausscoeff[1] + gausscoeff[2]
fit_str = 'amplitude = ' + string(gausscoeff[0], format='(e11.4)')
xyouts, notesxcoord, gausscoeff[0]*1.05, fit_str, /data, align=0, color=240, charsize=1.2
fit_str = 'centre = ' + string(gausscoeff[1], format='(e11.4)')
xyouts, notesxcoord, gausscoeff[0]*0.95, fit_str, /data, align=0, color=240, charsize=1.2
fit_str = 'width = ' + string(gausscoeff[2], format='(e11.4)')
xyouts, notesxcoord, gausscoeff[0]*0.85, fit_str, /data, align=0, color=240, charsize=1.2

print,'Fitted Coefficient '+mapsinfo[i,2]+': ',gausscoeff[1],'  +/-',gausscoeff[2] 

DEVICE, /close
ENDFOR

;Chi Squared Histogram
set_plot, 'ps'
DEVICE, file= graphDir + saveName + '_redchisq', /color, /landscape

plothist, red_chisq_results, xbin4, hist4, /autobin, xtitle = 'Coefficient', ytitle = 'Frequency',  title='Monte Carlo, n= ' + string(montecarlo, format='(d8.0)') + ', Reduced chi-squareds', psym=1, charsize=1

;Gaussian Fit
yfit4 = gaussfit(xbin4, hist4, gausscoeff4, nterms=3)


xarray = arrgen((gausscoeff4[1]-10*gausscoeff4[2]), (gausscoeff4[1]+10*gausscoeff4[2]), 0.01*gausscoeff4[2])
yarray = arrgen(0, montecarlo, 0.1)
zsigma4 = (xarray-gausscoeff4[1])/gausscoeff4[2]
ymodel = gausscoeff4[0] * exp(-(zsigma4^2)/2)
centreline = yarray
centreline[*] = gausscoeff4[1]

oplot, xarray, ymodel, color=240
oplot, centreline, yarray, color=50, linestyle=2
notesxcoord= gausscoeff4[1] + gausscoeff4[2]

fit_str = 'amplitude = ' + string(gausscoeff4[0], format='(d10.4)')
xyouts, notesxcoord, gausscoeff4[0]*1.05, fit_str, /data, align=0, color=240, charsize=1.4
fit_str = 'centre = ' + string(gausscoeff4[1], format='(d10.4)')
xyouts, notesxcoord, gausscoeff4[0]*0.95, fit_str, /data, align=0, color=240, charsize=1.4
fit_str = 'width = ' + string(gausscoeff4[2], format='(d10.4)')
xyouts, notesxcoord, gausscoeff4[0]*0.85, fit_str, /data, align=0, color=240, charsize=1.4

DEVICE, /close

print,''
print,'Reduced Chi-Squared:',gausscoeff4[1]
print,''

FOR i=0L,ntemplates-1 DO print,'Coefficient '+mapsinfo[i,2]+': ',coeffs[i],'  +/-',error_results[0,i]

ENDIF ELSE BEGIN

;---------------------------------
;Text files for results from data
;---------------------------------

;If no histograms just output results to text file
;Only print the top line once
IF datacount EQ 0 THEN BEGIN
printf,unit1,'#Region:',region[0],format='(A,I)'
printf,unit1,'#Data:              ',format='(A,$)'
FOR i=0L,ntemplates-1 DO BEGIN
printf,unit1,mapsinfo[i,2],format='(A,$)'
printf,unit1,'         Error         ',format='(A,$)'
ENDFOR
printf,unit1,'Reduced Chi-Squared'
ENDIF

;Output results
printf,unit1,datainfo[datacount,2],format='(A,$)'
printf,unit1,' ',format='(A,$)'
printf,unit1,datainfo[datacount,3],format='(A,$)'
printf,unit1,'   ',format='(A,$)'
FOR i=0L,ntemplates-1 DO BEGIN
printf,unit1,a[i],format='(e11.4,$)'
printf,unit1,'    ',format='(A,$)'
printf,unit1,a_err[i],format='(e11.4,$)'
printf,unit1,'    ',format='(A,$)'
ENDFOR
printf,unit1,red_chisq,format='(e11.4)'

IF datacount EQ 0 THEN BEGIN
printf,unit2,'#Region:',region[0],format='(A,I)'
printf,unit2,'#Data:              ',format='(A,$)'
FOR i=0L,ntemplates-1 DO BEGIN
printf,unit2,mapsinfo[i,2],format='(A,$)'
printf,unit2,'         Error         ',format='(A,$)'
ENDFOR
printf,unit2,''
ENDIF

;Output Spectral results
printf,unit2,datainfo[datacount,2],format='(A,$)'
printf,unit2,' ',format='(A,$)'
printf,unit2,datainfo[datacount,3],format='(A,$)'
printf,unit2,'   ',format='(A,$)'
FOR i=0L,ntemplates-1 DO BEGIN
printf,unit2,beta_results[i,0],format='(e11.4,$)'
printf,unit2,'    ',format='(A,$)'
printf,unit2,beta_error[i],format='(e11.4,$)'
printf,unit2,'    ',format='(A,$)'
ENDFOR 
printf,unit2,''
ENDELSE

;If in silent mode output percentage completion
IF keyword_set(silent) THEN BEGIN
percentage = 100.*((double(ndata)*double(regioncount))+(double(datacount)+1.))/(double(nregions)*double(ndata))
print,string(13b),'Template Fitting: ',percentage,'%',format='(A,A,F5.1,A,$)'
ENDIF

;END OF DATA LOOP
ENDFOR  
;Formating
IF keyword_set(data)  THEN BEGIN
printf,unit1,''
printf,unit2,''
ENDIF

;-----------------
;Graphs if wanted
;-----------------
IF keyword_set(dographs) THEN BEGIN
;Setup file
set_plot,'ps'

file_mkdir,graphDir+savename
file = graphDir+savename+'/region_'+string(region[0],format='(3I0)')
DEVICE,file=file,/color,/landscape

;Set up multiplot depending on how many templates there are
CASE  ntemplates OF
   1:!P.MULTI=0
   2:!P.MULTI=[0,1,2]
   3:!P.MULTI=[0,2,2]
   4:!P.MULTI=[0,2,2]
ENDCASE
;Get x values for setting 
xrange=minmax(datainfo[*,3])
;Plot all of the templates
FOR i=0L,ntemplates-1 DO ploterror,graph_freqs,graph_results[*,i],graph_errors[*,i],xtitle='Frequency (GHz)',ytitle='Template Coefficient',title=mapsinfo[i,2],psym=1,xrange=[xrange[0]-10,xrange[1]+10],color=250,errcolor=250,charsize=1.

;Close device
DEVICE,/close
ENDIF 
;END OF REGION LOOP
ENDFOR 

;Close files if needed
IF keyword_set(data) THEN BEGIN 
free_lun,unit1 
free_lun,unit2
ENDIF
print,''
print, 'DONE!'
END





   

