; GvH:
; fitout[*,7,jj] = delta[indx,0]
;
;+
; NAME:
;   DEMOD.PRO
;
; PURPOSE:
;   This procedure is contains a general demodulation algorithm for
;   spectropolarimetric meauserments performed using the technique of
;   spectral modulation.
;
; CATEGORY:
;   ALGORITHM
;
; CALLING SEQUENCE:
;   DEMOD, Input, Lambdainp, Spectra 
;
; INPUTS:
;   Input:      structure describing the multiple-order retarder
;   Lambdainp:  array containing the wavelengths corresponding to the data
;   Spectra:    array[nlamb, n, nspect] containing the spectra to be demodulated
;               nlamb = number of array indices corresponding to Lambdainp
;               n = 1 or 2, depending of the use of the single or double beam method
;               nspect = number of spectra to be demodulated
;   
; KEYWORDS:
;   CONF:       array[3] containing parameters that set the configuration of the SPEX optics
;   LRANGE:     array[2] specifying the minimum and maximum range of the wavelength spectrum
;   AOLP:       scalar specifying the AoLP of the reference spectrum
;   NLAMB:      scalar specifying the number of points that will be used to approximate, by
;               interpolation, the offset, DoLP, retardance and AoLP in the retardance-fit             
;   WINSTEP:    scalar specifying the number of fits per spectral window in case of spectral
;               window fitting; defaults = 0 (all spectral windows)
;   DOPLOT:     boolean specifying if intermediate output will be plotted (1: yes, 0: no)
;
; OUTPUTS:
;   Fitout:     structure containing all relevant input parameters of SETES
;
; OUTPUT KEYWORDS:
;   DELTALIT:   array containing the spectral retardance calculated from the input structure
;               using literature parameters
;
; DESCRIPTION:
;   This routine performs the demodulation of spectrally modulated spectra.
;   As input must be specified the composition of the multiple-order retarder in order, the
;   wavelength array to which the data in sampled, and the spectrally modulated spectra.
;   The latter can be specified using either the single or double beam method. The first
;   spectrum must be a calibration spectrum with known AoLP and preferably high DoLP, in order
;   to be able to accurately fit the spectral retardance. All other spectra are fitted using
;   the spectral fit method and using the fitted retardance as a fixed parameter.
;   
; EXAMPLE:
;
; INTERNAL ROUTINES
;   OPTICSDATA:       function that returns a structure containing literature data of
;                     optical materials
;   DEMOD_REFRINDX:   function that returns the refractive index of a given material on
;                     a given wavelength array
;   DEMOD_GET_DELTA:  function that returns the spectral retardance of a given set of crystals
;   DEMOD_SPC_WIN:    function that returns the length and edge positions of the spectral windows
;                     for a given spectral retardance
;   DEMOD_SPC_NORM:   function that returns normalized modulated spectra
;   DEMOD_FUNC1:      fit-function for performing full spectrum fit 
;   DEMOD_DELTA:      function that performs a full spectrum fit, including the spectral retardance
;   DEMOD_SPCWINFIT:  function that performs a spectral window fit
;   
; MODIFICATION HISTORY:
;   Written by: J.H.H. Rietjens, 01/12/2011
;   Dec., 2011
;-
; Copyright (C) 2011 SRON Netherlands Institute for Space Research
; Jeroen H. H. Rietjens
; j.h.h.rietjens@sron.nl
; 
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

;------------------------------------------------------------------------------

FUNCTION OPTICSDATA, lambda
nlambda = N_ELEMENTS(lambda)

optics = REPLICATE({name:'', N:DCOMPLEXARR(nlambda,2) $  ; Refractive index
                           , dNdT:DBLARR(8,2) $          ; Thermal optical coefficient
                           , T:DCOMPLEXARR(nlambda) $    ; Transmission
                           , tec:DBLARR(2)},10)           ; Thermal expansion coefficient
                           
optics.name = [ 'al2o3', 'bab2o4', 'caco3', 'mgf2', 'sio2', 'tio2', 'yvo4', '', '', 'artif' ]

;------------------------------------------------------------------------------
; Sellmeier equation: N[lambda] = SQRT(P1 + P2*lambda^2 / (lambda^2 - P3) $
;                                         + P4*lambda^2 / (lambda^2 - P5) $
;                                         + ...
;                                         + P(N-1)*lambda^2 / (lambda^2 - PN) )
; N[*,0] : ordinairy refractive index
; N[*,1] : extraordinairy refractive index
;------------------------------------------------------------------------------

; Al2O3, Sapphire, 0.22 - 5.0 micron
; Malitson, I. H. and Dodge, M. J., Refractive index and birefringence of synthetic sapphire,
; J.Opt. Soc. Am. 62, 1405 (1972).
optics[0].N[0:6,0] = [1., 1.43134936, 0.0726631^2, 0.65054713, 0.1193242^2, 5.3414021, 18.028251^2]
optics[0].N[0:6,1] = [1., 1.5039759, 0.0740288^2, 0.55069141, 0.1216529^2,6.59273791, 20.072248^2]
optics[0].dNdT[*,0] = [293., 1.755, -45.2665E-06, 83.5457E-06, 8.27, 5.85,  7.21, -2.4]
optics[0].dNdT[*,1] = [293., 1.748, -39.8961E-06, 81.9579E-06, 8.00, 5.42,  6.47, -2.2]
optics[0].tec = 1E-6*[7.21, 6.47]
;optics[0].dndtor = 13.6*1e-6 ; @ 589 nm, RT
;optics[0].dndtex = 14.7*1e-6 ; @ 589 nm, RT

;; BaB2O4, BBO, 0.22 - 1.06 micron
;optics[1].N[0:5,0] = SQRT(2.7405 + 0.0184*lambda^2 / (lambda^2 - 0.0179) - 0.0155 * lambda^2)
;optics[1].N[0:5,1] = SQRT(2.3730 + 0.0128*lambda^2 / (lambda^2 - 0.0156) - 0.0044 * lambda^2)

; beta-BaB2O4, BBO, 0.22 - 1.06 micron - dn/dT values of beta-BBO
; Eimerl, D., Davis, L., and Velsko, S., Optical, mechanical, and thermal properties of barium borate,
; J. Appl. Phys. 62, 1968 (1987).
; Ghosh, G, Temperature dispersion of refractive indices in βBaB2O4 and LiB3O5 crystals for nonlinear optical devices
; J. Appl. Phys. 78, 6752 (1995)
optics[1].N[0:4,0] = [1., 1.73651, 0.0120649, 0.0758505, 264.897]
optics[1].N[0:4,1] = [1., 1.36847, 0.0100179,   1.29495, 187.560]
optics[1].dNdT[*,0] = [293., 1.610, -19.3007E-6, -34.9683E-6, 19.00, 6.43,  4.00,  1.40]
optics[1].dNdT[*,1] = [293., 1.520, -141.421E-6, 110.8630E-6, 17.00, 6.43, 36.00, -5.40]
optics[1].tec = 1E-6*[4.00, 36.00]

; CaCO3, Calcite, 0.2 - 2.2 micron
; Gray, D. E., (Ed.), American Institute of Physics Handbook, 3rd ed. (McGraw-Hill, New York, 1972).
optics[2].N[0:8,0] = [1., 0.8559, 0.0588^2, 0.83913, 0.141^2, 0.0009, 0.197^2, 0.6845, 7.005^2]
optics[2].N[0:6,1] = [1., 1.0856, 0.07897^2, 0.0988, 0.142^2, 0.317, 1.468^2]
optics[2].dNdT[*,0] = [293., 1.613, -121.689E-06, 122.494E-06, 10.80, 10.00, 25.00, -7.60]
optics[2].dNdT[*,1] = [293., 1.472,  12.7011E-06, 20.4803E-06,  9.05,  6.83, -3.70, -1.20]
optics[2].tec = 1E-6*[25.00, -3.70]

; MgF2, 0.4 - 3.1 micron
; Dodge, M. J., Refractive properties of magnesium fluoride,
; Appl. Opt. 23, 1980 (1984).
optics[3].N[0:6,0] = [1., 0.48755108, 0.04338408^2, 0.39875031, 0.09461442^2, 2.3120353, 23.793604^2]
optics[3].N[0:6,1] = [1., 0.41344023, 0.03684262^2, 0.50497499, 0.09076162^2, 2.4904862, 23.771995^2]
optics[3].dNdT[*,0] = [293., 1.290, -37.2043E-06, 39.3186E-06, 13.10, 8.00,  9.3, -4.70]
optics[3].dNdT[*,1] = [293., 1.290, -56.7859E-06, 57.3986E-06, 15.50, 8.00, 14.2, -6.90]
optics[3].tec = 1E-6*[9.3, 14.2]
; lattice/ionic contribution: (o: L=15.9824E-6, l_ip = 40.0) (e: L=8.9996E-6, l_ip = 40.0)
;optics[3].dndtor = 1.12*1e-6 ; @ 633 nm, RT
;optics[3].dndtex = 0.58*1e-6 ; @ 633 nm, RT

; SiO2, Quartz, 0.18 - 0.71 micron
; Radhakrishnan, T., Further studies on the temperature variation of ; the refractive index of crystals,
; Proc. Indian Acad. Sci., A33, 22 (1951).
;optics[4].N[0:10,0] = [1., 0.663044, 0.060^2, 0.517852, 0.106^2, 0.175912, 0.119^2, 0.565380, 8.844^2, 1.675299, 20.742^2]
;optics[4].N[0:10,1] = [1., 0.665721, 0.060^2, 0.503511, 0.106^2, 0.214792, 0.119^2, 0.539173, 8.792^2, 1.807613, 197.70^2]

; SiO2, Quartz, 0.18 - 0.71 micron
; Gosh 1999 (Halle)
optics[4].N[0:4,0] = [1.28604141, 1.07044083, 1.00585997*0.01, 1.10202242, 100.]
optics[4].N[0:4,1] = [1.28851804, 1.09509924, 1.02101864*0.01, 1.15662475, 100.]
optics[4].dNdT[*,0] = [293., 1.515, -61.184E-06, 43.9990E-06, 10.30, 8.90,  6.88, -3.02]
optics[4].dNdT[*,1] = [293., 1.520, -70.1182E-06, 49.2875E-06, 10.30, 8.90, 12.38, -3.32]
optics[4].tec = 1E-6*[6.88, 12.38]
; fit expansion coef. 15.70 (meas. 6.88) and 17.80 (meas. 12.38)

; TiO2, Rutile, 0.43 - 1.5 micron
; DeVore,J. R., Refractive index of rutile and sphalerite,
; J. Opt. Soc. Am. 41, 416 (1951).
optics[5].N[0:4,0] = [1., 4.913, 0.0, 0.2441, 0.0803]
optics[5].N[0:4,1] = [1., 6.097, 0.0, 0.3322, 0.0843]
optics[5].dNdT[*,0] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
optics[5].dNdT[*,1] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
optics[5].tec = 1E-6*[8.98, 6.87]

; YVO4, 0.5 - 1.06 micron, HoOM Ref 112/113 - dn/dT values of TiO2!
; Maunder, E. A. and DeShazer, L. G., Use of yttrium orthovanadate for ; polarizers,
; J. Opt. Soc. Am. 61, 684A (1971).
; Lomheim, T. S. and DeShazer, L. G., Optical absorption intensities
; of trivalent neodymium in the uniaxial crystal yttrium orthovanadate,
; J. Appl. Phys. 49, 5517 (1978).
optics[6].N[0:2,0] = [1., 2.7665, 0.026884]
optics[6].N[0:2,1] = [1., 3.5930, 0.032103]
optics[6].dNdT[*,0] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
optics[6].dNdT[*,0] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
optics[6].tec = 1E-6*[8.98, 6.87]

; Artificial material with spectrally constant refractive index of 1.5 and birefringence of 0.01
optics[9].N[0,0] = [2.25]
optics[9].N[0,1] = [2.2801]
optics[9].dNdT[*,0] = [293., 1.5, 0.0E-06, 0.0E-06, 4.10, 3.50,  8.98, -0.46]
optics[9].dNdT[*,0] = [293., 1.51, 0.0E-06, 0.0E-06, 4.10, 3.50,  6.87, -0.26]
optics[9].tec = 1E-6*[0.0, 0.0]

RETURN, optics

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_REFRINDX, lambda, params1, params2, temp

nlambda = N_ELEMENTS(lambda)
lamb = lambda/1000.

Nabs = DCOMPLEXARR(nlambda,2)
dNabsdT = DBLARR(nlambda,2)

; Calculate refractive index as a function of wavelength via Sellmeier equation
Nabs[*,0] += params1[0,0]
Nabs[*,1] += params1[0,1]



for jj = 0, (N_ELEMENTS(WHERE(params1[*,0] ne 0))-1)/2 - 1 do begin
   Nabs[*,0] = Nabs[*,0] + (params1[2*jj+1,0]*lamb^2. / (lamb^2. - params1[2*jj+2,0]))
   Nabs[*,1] = Nabs[*,1] + (params1[2*jj+1,1]*lamb^2. / (lamb^2. - params1[2*jj+2,1]))
endfor

Nabs = SQRT(Nabs)

; Calculate temperature correction w.r.t. 293 K
dT = temp - params2[0,0]




if (params2[7,0] eq 0) then begin
  ; Using Schott's Sellmeier-type equation for glasses
  dNabsdT[*,0] = (Nabs[*,0]^2. - 1.) / (2.*Nabs[*,0]) * (params2[1,0]*dT + params2[2,0]*dT^2. + params2[3,0]*dT^3. $
             + (params2[4,0]*dT + params2[5,0]*dT^2.) / (lamb^2. - (params2[6,0])^2.) )

  dNabsdT[*,1] = (Nabs[*,1]^2. - 1.) / (2.*Nabs[*,1]) * (params2[1,1]*dT + params2[2,1]*dT  ^2. + params2[3,1]*dT^3. $
             + (params2[4,1]*dT + params2[5,1]*dT^2.) / (lamb^2. - (params2[6,1])^2.) )

endif else begin
  ; Using dNdT equation of Gosh
  Econv = 1.E9 * 6.626E-34 * 2.9979E8 / 1.6E-19
  lambdaig = Econv / params2[4,*]
  Ro = lambda^2. / (lambda^2. - lambdaig[0]^2.)
  Re = lambda^2. / (lambda^2. - lambdaig[1]^2.)
  dNabsdT[*,0] = (params2[2,0] * Ro + params2[3,0] * Ro^2) / (2. * Nabs[*,0])
  dNabsdT[*,1] = (params2[2,1] * Re + params2[3,1] * Re^2) / (2. * Nabs[*,1])
endelse

N = Nabs + dNabsdT*dT

RETURN, N

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_GET_DELTA, input, lambda, optics, N = n

; Complement input structure
tags = TAG_NAMES(input)
comp = (WHERE(tags eq 'COMP') eq -1) ? 'al2o3' : input.comp
ncomp = N_ELEMENTS(comp)
thick = (WHERE(tags eq 'THICK') eq -1) ? FLTARR(ncomp)+1 : input.thick
par = (WHERE(tags eq 'PAR') eq -1) ? FLTARR(ncomp > 2,2)+0 : input.par
rot = (WHERE(tags eq 'ROT') eq -1) ? FLTARR(ncomp)+0 : input.rot
temp = (WHERE(tags eq 'TEMP') eq -1) ? FLTARR(ncomp)+293 : input.temp


if (N_ELEMENTS(N) eq 0) then Nset = 0 else Nset = 1



nlambda = N_ELEMENTS(lambda)
delta = FLTARR(nlambda)

; Loop over optical components
for jj = 0, ncomp-1 do begin
  print,comp[jj]
  ; Calculate refractive index
  params1 = Nset ? 0 : optics[WHERE(optics.name eq comp[jj])].N
  params2 = Nset ? 0 : optics[WHERE(optics.name eq comp[jj])].dNdT
  N1 = Nset ? N : DEMOD_REFRINDX(lambda, params1, params2, temp[jj])
  ;print,'0',N1[0,0]
  ;print,'1',N1[0,1]
  ; Calculate birefringence
  biref1 = REFORM(N1[*,1] - N1[*,0])
  ;print,biref1[0:2]
  ; Calculate thermal expansion
  a = optics[WHERE(optics.name eq comp[jj])].tec[0]
  
  d1 = 1E6*thick[jj]*(1.+a*(temp[jj]-293.))
  ;print,d1
  ; Calculate retardance including field-of-view behavior and thermal expansion
  phi = par[jj,0]
  print,par[jj,1] ,rot[jj]
  theta = (biref1[0] gt 0) ? par[jj,1] + rot[jj]: par[jj,1] + !pi/2. + rot[jj]
  
  c1 = SIN(phi)*SIN(theta)
  c2 = SIN(phi)*COS(theta)
  p0 = N1[*,0] * SQRT(1 - (c1/N1[*,0])^2 - (c2/N1[*,0])^2)
  p1 = N1[*,1] * SQRT(1 - (c1/N1[*,0])^2 - (c2/N1[*,1])^2)
  delta1 = COS(!pi*rot[jj]/90.) * (biref1 * d1 / lambda) * ( p1 - p0 ) / biref1
  
  ; Add retardance of each component        
  delta = delta + delta1
  ;print,(delta*lambda)[0:5]

endfor

RETURN, (1-0.0*FINDGEN(nlambda)/(nlambda-1))*DOUBLE(delta*lambda)

end 

;------------------------------------------------------------------------------

FUNCTION DEMOD_SPC_WIN, lambda, delta

nlambda = N_ELEMENTS(lambda)
spcwin = {min:FLTARR(nlambda), max:FLTARR(nlambda), length:FLTARR(nlambda)}

; Calculate spectral resolution of wavelength array
res = ROUND(100./(lambda[(FINDGEN(nlambda)+1) < (nlambda - 1)]-lambda))/100.
res[nlambda-1] = res[nlambda-2]

; Calculate spectral window length
;length = res*lambda^2/(ABS(delta)*(1 + (lambda/(2*ABS(delta)))^2))
; Calculate position of left and right boundary of each spectral window
;spcwin.min = ROUND(FINDGEN(nlambda) - length/2) > 0.
;spcwin.max = ROUND(FINDGEN(nlambda) + length/2) < (nlambda-1)
;spcwin.length = spcwin.max - spcwin.min + 1

; GvH boerenlullen methode
k = delta/lambda
for i=0, nlambda-1 do begin
   spcwin.min[i] = where(abs(k-(k[i]+0.5)) eq min(abs(k-(k[i]+0.5))))
   spcwin.max[i] = where(abs(k-(k[i]-0.5)) eq min(abs(k-(k[i]-0.5))))
endfor
spcwin.length = spcwin.max-spcwin.min+1

RETURN, spcwin

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_FUNC0, p, XDATA = xdata, YDATA = ydata, FITWAV = fitwav

npar = N_ELEMENTS(p)
trat = INTERPOL(p[0*npar:1*npar-1], fitwav, xdata, /SPLINE)
Iest = [trat*REFORM(ydata[*,0])] + [REFORM(ydata[*,1])]

chis = (Iest - SMOOTH(Iest, 30))

RETURN, chis

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_SPC_NORM, spcwin, lambda, spectra, DOPLOT = doplot, TRAT = trat

nlambda = N_ELEMENTS(lambda)
lmin = MIN(WHERE(lambda ge (lambda[0]+10) ) )
lmax = MIN(WHERE(lambda ge (lambda[nlambda-1]-20) ) )
nspectra = N_ELEMENTS(spectra)/nlambda
spcnorm = FLTARR(nlambda,2)

; Calculate second spectrum in case of single beam data
if (nspectra eq 1) then begin
       Itmp = FLTARR(nlambda)
       Iest = FLTARR(nlambda)
;       for jj=0L, nlambda-1 do Itmp[jj] = 2.*MEAN(spectra[spcwin.min[jj]:spcwin.max[jj]])
       for jj=0L, nlambda-1 do Itmp[jj] = $
          max(spectra[spcwin.min[jj]:spcwin.max[jj]]) + $
          min(spectra[spcwin.min[jj]:spcwin.max[jj]])

       Iest = SMOOTH(REFORM(Itmp), MAX(spcwin.length)/2, EDGE_TRUNCATE = 1) > 0.001
       ;Iest = Itmp
       spectra = [ [spectra] , [Iest - REFORM(spectra[*,0]) ] ]
endif



; Calculate intensity spectrum
Iest = spectra[*,0] + spectra[*,1]
print,'before',MEAN(Iest)
; Calculate normalized spectra
spcnorm = spectra / [ [Iest] , [Iest] ]

; Calculate transmission ratio
normtmp = FLTARR(nlambda)

;;;; GvH for jj=0L, nlambda-1 do normtmp[jj] = MEAN(spcnorm[spcwin.min[(jj-1)>0]:spcwin.max[(jj+1)<nlambda-1],0])
;for jj=0L, nlambda-1 do normtmp[jj] = MEAN(spcnorm[spcwin.min[(jj+1)<(nlambda-1)]:spcwin.max[jj<(nlambda-1)],0])
;normfitcoef = POLY_FIT(lambda[lmin:lmax], normtmp[lmin:lmax], 3, YFIT = normfit)
;normfit = normfitcoef[0] + normfitcoef[1]*lambda + normfitcoef[2]*lambda^2 + normfitcoef[3]*lambda^3
;normfitcoef = POLY_FIT(lambda[lmin:lmax], normtmp[lmin:lmax], 5, YFIT = normfit) ; GvH higher order
;normfit = poly(lambda, normfitcoef) ; GvH higher order
;trat = REFORM(1./normfit - 1.)


;++++++++++++++++++++++++++++++++++++++++++
; GvH ITERATE TRANSMISSION RATIO CORRECTION
trat = make_array(n_elements(lambda),value=1d)


pmulti = !p.multi
!P.MULTI = [0, 1, 2]
if doplot ne 0 then window, 10, xsize=600, ysize=900
;cgplot, lambda, trat, /xs,yrange=[0.9,1.1]
for t=0, 20 do begin
   print,t
   ; Correct spectra
   spectras = [ [trat*REFORM(spectra[*,0])] , [REFORM(spectra[*,1])] ]
   ; Recalculate intensity spectrum
   Iest = spectras[*,0] + spectras[*,1]
   ; Recalculate normalized spectra
   spcnorm = spectras / [ [Iest] , [Iest] ]
   ; Calculate transmission ratio
; GvH   for jj=0L, nlambda-1 do normtmp[jj] = MEAN((spcnorm)[spcwin.min[(jj+1)>0]:sp;cwin.max[jj<(nlambda-1)],0])


   ; as soon as t is not too extreme anymore, it's better to use the modulation mean
   ; especially at higher DoLP, the (max-min)/2 method breaks down (see Mathematica nb).
   ; in the first few iterations we use (max-min)/2, then the mean.
   if (t le 999) then begin
      ; average is not the center of the modulation, because modulation is not perfect
      ; cosine; it is 'spiky', so the center of mass is below the center of modulation!
      ; Now let's calculate the center using (max+min)/2:
      ; warning: use with caution: it will wash out contrast in O2A band!!!
      for jj=0L, nlambda-1 do normtmp[jj] = ($
         max(spcnorm[(spcwin.min[jj]+1)>0:spcwin.max[jj]<(nlambda-1),0]) + $
                       min(spcnorm[(spcwin.min[jj]+1)>0:spcwin.max[jj]<(nlambda-1),0]) )/2
   endif else begin
      for jj=0L, nlambda-1 do normtmp[jj] = MEAN(spcnorm[(spcwin.min[jj]+1)>0:$
                                                         spcwin.max[jj]<(nlambda-1),0])
   endelse
   
   ; spectral smooth
      for jj=0L, nlambda-1 do normtmp[jj] = mean( $
      normtmp[ (spcwin.min[jj]+1)>0 : spcwin.max[jj]<(nlambda-1),0] )
;   trat *= REFORM(1./normtmp - 1.)
   
   ; or fit
   normfitcoef = POLY_FIT(lambda[lmin:lmax], normtmp[lmin:lmax], 5, YFIT = normfit) ; GvH hig;her order
   
   normfit = poly(lambda, normfitcoef)                                              ; GvH hig;her order
   if t lt 999 then trat *= REFORM(1./normfit - 1.) else trat *= REFORM(1./normtmp - 1.)

   for jj=0L, nlambda-1 do trat[jj] = mean( $
      trat[ (spcwin.min[jj]+1)>0 : spcwin.max[jj]<(nlambda-1),0] )
   print,mean(Iest),mean(spectras[*,0])
if doplot ne 0 then begin
   ;cgplot, lambda, trat, /overplot, thick=t
;   cgplot, lambda, spcnorm[*,0], /xs, yrange=[0,1], title='t='+strtrim(t,2)
   cgplot, lambda, spcnorm[*,0], /xs, yrange=[0,1], /ys, $
           title='t='+strtrim(t-1,2)
   cgplot, lambda, spcnorm[*,1], color='red', /overplot
;   cgplot, lambda, make_array(n_elements(lambda),value=0.5), /overplot
   cgplot, lambda, normfit, /overplot, color='green'
   cgplot, lambda, normtmp, /overplot, color='yellow'

   cgplot, lambda, spectras[*,0], /xs, /ys
   cgplot, lambda, spectras[*,1], /overplot, color='red'
   cgplot, lambda, Iest/2d, /overplot, color='purple'
   wait, 1.0
endif
endfor
!p.multi = pmulti
;++++++++++++++++++++++++++++++++++++++++++


if (N_ELEMENTS(trat) eq 0) then begin
print,'HERE'
; Definition of fit-parameters
nlamb = 10
fitpar1 = REPLICATE({fixed:0, limited:[0,0], limits:[0.D,0.D]},nlamb)
fitpar1[0*nlamb:1*nlamb-1].limited = [1, 1]
fitpar1[0*nlamb:1*nlamb-1].limits  = [0.5, 1.5]

; Initialize fit
nlambd = nlamb + 2
fitwav = lambda[[0,FINDGEN(nlamb-2)*0.98*nlambda/(nlamb-3)+0.01*nlambda, nlambda-1]]
fittrat = trat[[0,FINDGEN(nlamb-2)*0.98*nlambda/(nlamb-3)+0.01*nlambda, nlambda-1]]
functargs1 = {XDATA:lambda, YDATA:spectra, FITWAV:fitwav}

; Full trat curvefit
fitparout = MPFIT('DEMOD_FUNC0', fittrat, functargs=functargs1, parinfo=fitpar1 $
              , BESTNORM = bestnorm, NITER = niter, STATUS = status, ERRMSG = errmsg, /QUIET)
if (status gt 0) then PRINT, niter, bestnorm, FORMAT='("Fit successful with ",I3," iterations and chi-sqr = ",F8.5)' $
  else PRINT, errmsg

trat = INTERPOL(fitparout[0*nlamb:1*nlamb-1], fitwav, lambda, /SPLINE)
endif

;trat /= 0.95 ; little experiment o2a paper... GvH
;trat /= trat ; little experiment o2a paper... GvH
; Correct spectra
spectras = [ [trat*REFORM(spectra[*,0])] , [REFORM(spectra[*,1])] ]

; Recalculate intensity spectrum
Iest = spectras[*,0] + spectras[*,1]

; Recalculate normalized spectra
spcnorm = spectras / [ [Iest] , [Iest] ]

; Plot if DOPLOT = 1
if doplot ne 0 then begin
  WSET, 0
  !P.MULTI = [0,3,3]
  PLOT, lambda, spcnorm[*,0] + spcnorm[*,1], /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Normalized spectra', YRANGE = [0,1]
  OPLOT, lambda, spcnorm[*,0], COLOR = 100 ;'0000ff'x
  OPLOT, lambda, spcnorm[*,1], COLOR = 200 ;'00ff00'x
  PLOT, lambda, normtmp, /XSTYLE, /YSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Window average', yrange=[0,1];, yrange=[0.3, 0.7]
  OPLOT, lambda, normfit, COLOR = 100 ; '0000ff'x
  PLOT, lambda, trat, /XSTYLE, /YSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Transmission ratio', yrange=[0,1];, yrange=[0.7, 1.3]
  OPLOT, lambda, spcnorm[*,0] + spcnorm[*,1], COLOR = 100 ; 'ff0000'x
  ;OPLOT, lambda, trat2, COLOR = '0000ff'x

endif

RETURN, spcnorm

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_FUNC1, p, XDATA = xdata, YDATA = ydata, FITWAV = fitwav, DELTA = delta, CONF = conf, LRANGE = lrange

npar = N_ELEMENTS(p)/4.

Off   = INTERPOL(p[0*npar:1*npar-1], fitwav, xdata, /SPLINE)
DoLP  = INTERPOL(p[1*npar:2*npar-1], fitwav, xdata, /SPLINE)
delta = INTERPOL(p[2*npar:3*npar-1], fitwav, xdata, /SPLINE) * delta
AoLP  = FLTARR(N_ELEMENTS(xdata)) + p[3*npar]

model = Off + (-1)^conf[2]*DoLP/2. * COS(2.*!pi*delta / xdata - (-1)^(conf[0]+conf[1])*2.*AoLP)

ydiff = ydata[lrange[0]:lrange[1]] - model[lrange[0]:lrange[1]]

RETURN, ydiff

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_DELTA, nlamb, delta, lambda, spectra, aolp, CONF = conf, FITRES = fitres, DOPLOT = doplot

if (N_ELEMENTS(conf) eq 0) then conf = [0, 0, 0]

nlambda = N_ELEMENTS(lambda)
lmin = MIN(WHERE(lambda ge (lambda[0]+10) ) )
;lmax = MIN(WHERE(lambda ge (lambda[nlambda-1]-20) ) )
lmax = MIN(WHERE(lambda ge (lambda[nlambda-1]-0) ) ) ; GvH
delta = REFORM(delta[*,0])

PRINT, 'Performing retardance-fit...'

; Definition of fit-parameters
fitpar1 = REPLICATE({fixed:0, limited:[0,0], limits:[0.D,0.D]},4*nlamb)
fitpar1[3*nlamb:4*nlamb-1].fixed = 1 ; calibration AoLP fixed (as it should be)
;fitpar1[3*nlamb:4*nlamb-1].fixed = 0 ; calibration AoLP not fixed, but spectrally constant
fitpar1[0*nlamb:1*nlamb-1].limited = [1, 1]
fitpar1[0*nlamb:1*nlamb-1].limits  = [0.45, 0.55]
fitpar1[1*nlamb:2*nlamb-1].limited = [1, 1]
fitpar1[1*nlamb:2*nlamb-1].limits  = [0., 1.]
fitpar1[2*nlamb:3*nlamb-1].limited = [1, 1]
;fitpar1[2*nlamb:3*nlamb-1].limits  = [0.98, 1.02] ; GvH
fitpar1[2*nlamb:3*nlamb-1].limits  = [0.995, 1.005]

; Initialize fit
nlambd = nlamb + 2
fitwav = lambda[[0,FINDGEN(nlamb-2)*0.98*nlambda/(nlamb-3)+0.01*nlambda, nlambda-1]]
fitstartdelta = 1
functargs1 = {XDATA:lambda, YDATA:spectra, FITWAV:fitwav, DELTA:delta, CONF:conf, LRANGE:[lmin,lmax]}
fitstart = [ FLTARR(nlamb) + 0.5, FLTARR(nlamb) + 0.8, FLTARR(nlamb) + fitstartdelta, FLTARR(nlamb) + aolp]

; Full spectrum curvefit
fitparout = MPFIT('DEMOD_FUNC1', fitstart, functargs=functargs1, parinfo=fitpar1 $
              , BESTNORM = bestnorm, NITER = niter, STATUS = status, ERRMSG = errmsg, /QUIET)
if (status gt 0) then PRINT, niter, bestnorm, FORMAT='("Fit successful with ",I3," iterations and chi-sqr = ",F8.5)' $
  else PRINT, errmsg

; GvH: discard unrealistic deviations at 580--670 nm and >740 nm
; GvH2: actually, this yields unrealistically fluctuating AoLPs, so let's use Jeroen's code.
; ind = where((fitwav lt 580) or (fitwav gt 670 and fitwav lt 740))

; Calculate fit-spectrum
fitout = FLTARR(nlambda, 4)
fitout[*,0] = INTERPOL(fitparout[0*nlamb:1*nlamb-1], fitwav, lambda, /SPLINE)
fitout[*,1] = INTERPOL(fitparout[1*nlamb:2*nlamb-1], fitwav, lambda, /SPLINE)
fitout[*,2] = INTERPOL(fitparout[2*nlamb:3*nlamb-1], fitwav, lambda, /SPLINE)
;fitout[*,2] = INTERPOL((fitparout[2*nlamb:3*nlamb-1])[ind], fitwav[ind], lambda) ; GvH
fitout[*,3] = FLTARR(nlambda) + fitparout[3*nlamb]
fitres = fitout[*,0] + (-1)^conf[2]*fitout[*,1]/2. * COS(2.*!pi*delta*fitout[*,2]/lambda $
         - (-1)^(conf[0]+conf[1])*2.*fitout[*,3])

; Plot if DOPLOT != 0
if doplot ne 0 then begin
  WINDOW, 2, xsize=600, ysize=675;, RETAIN=2 ; GvH
  !P.MULTI = [0,1,3]
  cgPLOT, lambda, INTERPOL(fitparout[2*nlamb:3*nlamb-1], fitwav, lambda, /SPLINE) $
      , XTITLE = 'Wavelength (nm)', YTITLE = 'Fit result delta', /XSTYLE, /YSTYLE;, CHARSIZE = !P.CHARSIZE/2.
  cgPLOT, lambda, fitout[*,2], /overplot, color='green'
  cgPLOT, lambda, spectra, /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Normalized signal', /YSTYLE, thick=3 ;, CHARSIZE = !P.CHARSIZE/2.
  cgPLOT, lambda, fitres, COLOR = 'blue', /overplot
  cgPLOT, lambda, spectra-fitres, /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Fit residuals', /YSTYLE, thick=3;, CHARSIZE = !P.CHARSIZE/2.
endif
!p.multi = 0


RETURN, delta*REFORM(fitout[*,2])

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_FUNC2, p, XDATA = xdata, YDATA = ydata, DELTA = delta, CONF = conf

nxdata = N_ELEMENTS(xdata)
xdatan = FINDGEN(nxdata)/(nxdata-1) - 0.5
DoLP = p[1] + p[2]*xdatan; + p[3]*xdatan^2 ; GvH set linear term off
model = p[0] + (-1)^conf[2] * DoLP/2. * COS(2.*!pi*delta / xdata - (-1)^(conf[0]+conf[1])*2.*p[3])

ydiff = ydata - model

RETURN, ydiff

end

;------------------------------------------------------------------------------

FUNCTION DEMOD_SPCWINFIT, delta, spcwin, lambda, spectra, CONF = conf, DOPLOT = doplot, AOLP = aolp, QUIET = quiet

if (N_ELEMENTS(conf) eq 0) then conf = [0, 0, 0]
if (N_ELEMENTS(aolp) eq 0) then aolp = 0.
nlambda = N_ELEMENTS(lambda)
lmin = MIN(WHERE(lambda ge (lambda[0]+10) ) )
lmax = MIN(WHERE(lambda ge (lambda[nlambda-1]-20) ) )
nindx = N_ELEMENTS(spcwin.length)

; Definition of fit-parameters
fitpar2 = REPLICATE({fixed:0, limited:[0,0], limits:[0.D,0.D], relstep:0.},4)
fitpar2[0].limited = [1, 1] ; GvH
fitpar2[0].limits = [0.4, 0.6] ; GvH
;fitpar2[0].fixed = 1 ; GvH
fitpar2[1].limited = [1, 1]
fitpar2[1].limits = [0.0001, 1.]
fitpar2[2].limited = [1, 1]
fitpar2[2].limits = [-1., 1.]
fitpar2[3].limited = [0, 0]
fitpar2[3].limits = [0., 2*!pi]
fitpar2[3].relstep = 0.
fitstart = [0.5, 0.8, 0.00001, aolp+!pi]
fitout = FLTARR(nindx,5)

if doplot ne 0 then begin
  !P.MULTI = 0
  ;WINDOW, 1, RETAIN = 2, xsize=600, ysize=450
  WINDOW, 1, xsize=600, ysize=450 ; GvH
endif

; Loop over all selected spectral windows
if doplot eq 2 then WSET, 1
for jj=0, nindx-1 do begin
  if (spcwin.length[jj] gt 5) then begin
    ; Initialize fit
     spcwinels = indgen(spcwin.max[jj]-spcwin.min[jj]+1)+spcwin.min[jj] ; exclude O2A data points for better continuum fit, 
     o2aels = where( (lambda ge 757) and (lambda le 764) )
     noo2a = where((spcwinels le o2aels[0]) or (spcwinels ge o2aels[n_elements(o2aels)-1])) ; and hence better line efficiency correction GvH
    functargs2 = {XDATA:lambda[spcwinels[noo2a]], YDATA:spectra[spcwinels[noo2a],0] $
                , DELTA:delta[spcwinels[noo2a]], CONF:conf}
;    functargs2 = {XDATA:lambda[spcwin.min[jj]:spcwin.max[jj]], YDATA:spectra[spcwin.min[jj]:spcwin.max[jj],0] $
;                , DELTA:delta[spcwin.min[jj]:spcwin.max[jj]], CONF:conf}
    ; Full spectrum curvefit
    fitout[jj,0:3] = MPFIT('DEMOD_FUNC2', fitstart, functargs=functargs2, parinfo=fitpar2 $
              , BESTNORM = bestnorm, NITER = niter, STATUS = status, ERRMSG = errmsg, /QUIET)
    if ~(quiet) then $
      if (status gt 0) then PRINT, niter, bestnorm, FORMAT='("Fit successful with ",I3," iterations and chi-sqr = ",F8.5)' $
        else PRINT, errmsg
    fitout[jj,4] = bestnorm
    ; Set start values for next fit
    fitstart = [REFORM(fitout[jj,0:3])]
    ; Plot if DOPLOT = 1
    if doplot eq 2 then begin
;      if (jj MOD 20) eq 0 then begin
   ;   WSET, 1
;      !P.MULTI = [6,3,3]
;      POLYFILL, [0,0,0.33,0.33], [0.33,0.66,0.66,0.33], COLOR = '000000'x, /NORMAL ; GvH super slow
;      !P.MULTI = [6,3,3]
      CGPLOT, lambda[spcwin.min[jj]:spcwin.max[jj]], spectra[spcwin.min[jj]:spcwin.max[jj],0], YRANGE = [0,1] $
          , XTITLE = 'Wavelength (nm)', YTITLE = 'Normalized signal', psym=2, /xs
      xdatan = FINDGEN(spcwin.length[jj])/(spcwin.length[jj]-1) - 0.5
      CGPLOT, lambda[spcwin.min[jj]:spcwin.max[jj]], fitout[jj,0] + (-1)^conf[2] * (fitout[jj,1] + fitout[jj,2]*xdatan)/2. $
             * COS(2*!pi*delta[spcwin.min[jj]:spcwin.max[jj]]/lambda[spcwin.min[jj]:spcwin.max[jj]] - (-1)^(conf[0]+conf[1]) $
             * 2.*fitout[jj,3]), COLOR = 'red', thick=3, /overplot
      wait, 0.01
;      endif
    endif
  endif
endfor

if doplot ne 0 then if n_elements(where(finite(fitout[*,4]) eq 1)) gt 1 then $
   PLOT, lambda[(spcwin.max+spcwin.min)/2.], fitout[*,4], /XSTYLE

RETURN, fitout[*,0:4]

end

;------------------------------------------------------------------------------

PRO DEMOD, input, lambdainp, spectra, fitout $
         , LRANGE = lrange, NODELTAFIT = nodeltafit, CONF = conf, AOLP = aolp, NLAMB = nlamb, WINSTEP = winstep $
         , DOPLOT = doplot, DOSAVE = dosave, DELTALIT = deltalit, QUIET = quiet

; Keyword check     
if (N_ELEMENTS(input) eq 0) then input = {comp:['al2o3','mgf2'], thick:[1.2,2.88]};{comp:['al2o3','mgf2'], thick:[1.1,2.46]}
if (N_ELEMENTS(conf) eq 0) then conf = [1, 0, 0]
if (N_ELEMENTS(lrange) eq 0) then lrange = [390, 820] else lrange = [lrange[0]-10,lrange[1]+20]
if (N_ELEMENTS(doplot) eq 0) then doplot = 0
if (N_ELEMENTS(dosave) eq 0) then dosave = 0
if (N_ELEMENTS(aolp) eq 0) then aolp = 0
if (N_ELEMENTS(winstep) eq 0) then winstep = 0
if (N_ELEMENTS(nodeltafit) eq 0) then nodeltafit = 0
if (N_ELEMENTS(quiet) eq 0) then quiet = 0

; Select data in wavelength range
nlambda = N_ELEMENTS(lambdainp)

;print,'lrange',lrange
;print,'lambdainp',lambdainp[0],lambdainp[nlambda-1]

if (lrange[0] lt lambdainp[0]) then begin
  PRINT, '   Warning: adjusting lambda_min to avoid edge effect: ' + STRING(lambdainp[0]+10, FORMAT='(I3)') + ' nm, (was ' $
         + STRING(lambdainp[0], FORMAT='(I3)') + ' nm)'
  lrange[0] = lambdainp[0]+10
endif
if (lrange[1] gt lambdainp[nlambda-1]) then begin
  PRINT, '   Warning: adjusting lambda_max to avoid edge effect ' + STRING(lambdainp[nlambda-1]-20, FORMAT='(I3)') + ' nm, (was ' $
         + STRING(lambdainp[nlambda-1], FORMAT='(I3)') + ' nm)'
  lrange[1] = lambdainp[nlambda-1]-20
endif
lmin = MIN(WHERE(lambdainp ge (lrange[0]-10)))
lmax = MIN(WHERE(lambdainp ge (lrange[1]+20)))
lminout = MIN(WHERE(lambdainp ge (lrange[0])))
lmaxout = MIN(WHERE(lambdainp ge (lrange[1])))
lambda = lambdainp[lmin:lmax]
nlambda = N_ELEMENTS(lambda)
spectra = spectra[lmin:lmax,*,*]


szspectra = SIZE(spectra, /DIMENSIONS)

if (N_ELEMENTS(szspectra) eq 3) then nspect = szspectra[2] else nspect = 1



; Setup plot window if DOPLOT = 1

;; if doplot ne 0 then begin
;;   sysvar = {p:!P, x:!X, y:!Y, z:!Z, w:!D.WINDOW, dev:!D.NAME}
;;   SET_PLOT, 'x', /COPY
;; ;  WINDOW, 0, RETAIN = 2, XSIZE = 900, YSIZE = 700
;;   WINDOW, 0, XSIZE = 900, YSIZE = 700 ; GvH
;;   !P.MULTI = [0,3,3]
;;   !P.CHARSIZE = 2.5
;;   device, decomposed=0
;;   cgloadct, 4, ncolors=254
;; endif

; Load optical literature parameters
optics = OPTICSDATA(lambda)

; Calculate spectral retardance based on input composition of multiple-order retarder
deltalit = DEMOD_GET_DELTA(input, lambda, optics)

if ~(quiet) then PRINT, MEAN(deltalit), FORMAT = '("Mean retardance = ",I7)'
; Calculate spectral window lengths and positions
spcwin = DEMOD_SPC_WIN(lambda, deltalit)


;plot,spectra[*,1,0]+spectra[*,0,0]
;oplot,spectra[*,1,1]+spectra[*,0,1]
;oplot,spectra[*,1,2]+spectra[*,0,2]


if (nodeltafit) then delta = deltalit else begin
  ; Calculate normalized spectrum of calibration spectrum
  spcnorm = DEMOD_SPC_NORM(spcwin, lambda, REFORM(spectra[*,*,0]), DOPLOT = 1)
  ; Calculate number of point to be used in the interpolation in the full spectrum fit
  if (N_ELEMENTS(nlamb) eq 0) then $ 
    nlamb = (ABS(CEIL(0.6*MEAN(deltalit)*(1./lambda[0]-1./lambda[N_ELEMENTS(lambda)-1]))) < 40) > 4
  ; Perform full spectrum fit in order to obtain best fit of the spectral retardance
  delta = DEMOD_DELTA(nlamb, deltalit, lambda, spcnorm[*,0], aolp, CONF = conf, FITRES = fitres, DOPLOT = doplot)
  ; Recalculate spectral window lengths and positions
  spcwin = DEMOD_SPC_WIN(lambda, delta)
endelse
stop
END




; Calculate positions at which to perform spectral window fitting
indx = INTARR(nlambda)
if (winstep gt 0) then begin
    jj = nlambda-spcwin.length[nlambda-1]/(2*winstep)
    while (jj gt 0) do begin
      indx[jj] = jj
      jj = jj - spcwin.length[jj]/(2*winstep)
    endwhile
  indx = indx[WHERE(indx gt 0)]
  indx = indx[UNIQ(indx, SORT(indx))]
  spcwinfit = {min:spcwin.min[indx], max:spcwin.max[indx], length:spcwin.length[indx]}
endif else begin
  spcwinfit = spcwin
  indx = INDGEN(N_ELEMENTS(lambda))
endelse
nindx = N_ELEMENTS(indx)

fitout = (nspect eq 1) ? FLTARR(nindx, 8) : FLTARR(nindx, 8, szspectra[2])
for jj=0, nspect-1 do begin
  spcnorm = DEMOD_SPC_NORM(spcwin, lambda, REFORM(spectra[*,*,jj]), DOPLOT = doplot, TRAT = trat)
  PRINT, jj+1, nspect, FORMAT='("Performing spectral window fitting ",I4,"/",I4)'
  fitout[*,0:4,jj] = DEMOD_SPCWINFIT(delta, spcwinfit, lambda, spcnorm $
                                 , CONF = conf, DOPLOT = doplot, AOLP = aolp, QUIET = quiet)
  fitout[*,5,jj] = lambda[indx]
  fitout[*,6,jj] = trat[indx,0]
  fitout[*,7,jj] = delta[indx,0]
; Plot if DOPLOT = 1

;; if doplot ne 0 then begin
;;   WSET, 0
;;   !P.MULTI = [3,3,3]
;;   PLOT, lambda[indx], fitout[*,0,jj], /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Offset / DoLP', YRANGE = [0,1]
;;   OPLOT, lambda[indx], fitout[*,1,jj], COLOR = 'ff6666'x
;; ;  for jj=0, nspect-1 do begin
;; ;     oplot, lambda[indx], fitout[*,0,jj]
;; ;     oplot, lambda[indx], fitout[*,1,jj], color = 'ff6666'x
;; ;  endfor
;;   PLOT, lambda[indx], fitout[*,2,jj], /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'Linear term DoLP', /YSTYLE
;; ;  for jj=0, nspect-1 do begin
;; ;     oplot, lambda[indx], fitout[*,2,jj]
;; ;  endfor
;;   PLOT, lambda[indx], fitout[*,3,jj], /XSTYLE, XTITLE = 'Wavelength (nm)', YTITLE = 'AoLP'
;; ;  for jj=0, nspect-1 do begin
;; ;     oplot, lambda[indx], fitout[*,3,jj]
;; ;  endfor
;;   WAIT, 0.001
;; endif
  
endfor
;fitout = fitout[lminout:lmaxout,*,*]

;; if doplot ne 0 then begin
;;   ; Restore graphics system variables
;;   !P = sysvar.p
;;   !X = sysvar.x
;;   !Y = sysvar.y
;;   !Z = sysvar.z
;;   SET_PLOT, sysvar.dev
;;   WSET, sysvar.w
;; endif

; Save output if DOSAVE = 1
; [lambda, indx, spectrum]
; indx: 0-offset, 1-DoLP, 2-Linear term DoLP, 3-AoLP, 4-Chi^2, 5:lambda
if dosave then SAVE, fitout, FILENAME = 'DEMOD_' + STRCOMPRESS(STRING(SYSTIME(/JULIAN), FORMAT='(F11.3)'),/REMOVE_ALL) + '.sav'

return

end  
