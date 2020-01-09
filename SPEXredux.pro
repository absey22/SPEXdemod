PRO SPEXredux, dir

;dir='/home/seymour/Documents/SPEX/demodulation_pipeline/SPEXredux/ScriptRun_09072013_105035'

;default data items as per GvH:
;'/home/seymour/Documents/SPEX/demodulation_pipeline/SPEXredux/ScriptRun_09072013_105035'

if strmid(dir,0,/reverse_offset) ne path_sep() then dir += path_sep()


;print, '   '
;print, 'specified dir:'
;print, dir




;function defining method changed in current version of IDL?
;SPEXredux crashes on compiling with these internal functions:

;;================
;;================ FUNCTIONS
;;================
;(moved to separate script poly2d.pro)
;FUNCTION poly2d, x, y, coefs
;
;z = make_array(n_elements(x), n_elements(y), /double)
;
;for xx=0, n_elements(x)-1 do begin
;   for yy=0, n_elements(y)-1 do begin
;      for i=0, (size(coefs))[2]-1 do begin
;         for j=0, (size(coefs))[1]-1 do begin
;            z[xx,yy] += coefs[j,i] *x[xx]^i *y[yy]^j
;         endfor
;      endfor
;   endfor
;endfor
;
;RETURN, z
;end


;(moved to separate script astm.pro)
;FUNCTION ASTM, lambdamin, lambdamax
;
;if (size(findfile('solarspec.sav')))[0] ne 0 then $
;   restore, 'solarspec.sav' else print, '*** error: solarspec.sav file not found ***'
;
;minel = (sort(abs(solar[*,0]-lambdamin)))[0]
;maxel = (sort(abs(solar[*,0]-lambdamax)))[0]
;
;irrad = mean(solar[minel:maxel,1])
;
;RETURN, irrad
;end








;;================ IMPORT META DATA
meta = file_search(dir+'Meta*.*')
temp = 0.
pantilt = make_array(2,9999,/float)
temperature = make_array(2,9999,/float)
texp = make_array(9999,/float)
openr, 1, meta

on_ioerror, abortreadf ; read file untill error occurs
   for i=0, 10 do readf, 1, temp ; header
   for sl=0, 999 do begin
      readf, 1, temp            ; ! Scriptline
      for pt=0, 1 do begin
         readf, 1, temp
         pantilt[pt,sl] = temp
      endfor
      readf, 1, temp
      texp[sl] = temp
      readf, 1, temp ; at scriptline
      for c=0, 1 do begin
         readf, 1, temp
         temperature[c,sl] = temp
      endfor
      readf, 1, temp ; ! RainSensor
   endfor
abortreadf: close, 1
pantilt = pantilt[*,0:sl-1]
temperature = temperature[*,0:sl-1]
texp = texp[0:sl-1]



;;================ IMPORT SPECTRA
; find and sort file names
specfiles = make_array(2,n_elements(texp),/string)
sline = make_array(n_elements(texp))
for c=0, 1 do begin
   specfiles[c,*] = file_search(dir+'Spectrometer_110516'+strtrim(c+1,2)+'U2_*pix.txt')
   print,c
   print, 'filename:',specfiles
   for sl=0, n_elements(texp)-1 do begin
      temp = strsplit(specfiles[c,sl],'_',/extract)
      sline[sl] = temp[n_elements(temp)-8]
   endfor
   specfiles[c,*] = specfiles[c,sort(sline)]
endfor
; import data
spec = make_array(3648,2,n_elements(texp))
temp = make_array(3648)
for sl=0, n_elements(texp)-1 do begin
   for c=0, 1 do begin
      openr, 1, specfiles[c,sl]
      readf, 1, temp
      spec[*,c,sl] = temp
      close, 1
   endfor
endfor



;;================ IMPORT DATE + TIMES
times = make_array(n_elements(texp),/float)
date = make_array(3,/float) ; [day,month,year]
for sl=0, n_elements(texp)-1 do begin
   temp = strsplit(specfiles[0,sl],'_',/extract)
   times[sl] = strmid(temp[n_elements(temp)-2],9,2) +$
               strmid(temp[n_elements(temp)-2],11,2) /60d +$
               strmid(temp[n_elements(temp)-2],13,2) /60d^2d
endfor
date[0] = strmid(temp[n_elements(temp)-2],0,2)
date[1] = strmid(temp[n_elements(temp)-2],2,2)
date[2] = strmid(temp[n_elements(temp)-2],4,4)




;;================ IMPORT BLACK PIXELS
; find and sort file names
blackfiles = make_array(2,n_elements(texp),/string)
sline = make_array(n_elements(texp))
for c=0, 1 do begin
   blackfiles[c,*] = file_search(dir+'Spectrometer_110516'+strtrim(c+1,2)+'U2_*pix.dark13.txt')
   for sl=0, n_elements(texp)-1 do begin
      temp = strsplit(blackfiles[c,sl],'_',/extract)
      sline[sl] = temp[n_elements(temp)-8]
   endfor
   blackfiles[c,*] = blackfiles[c,sort(sline)]
endfor
; import data
black = make_array(13,2,n_elements(texp))
temp = make_array(13)
for sl=0, n_elements(texp)-1 do begin
   for c=0, 1 do begin
      openr, 1, blackfiles[c,sl]
      readf, 1, temp
      black[*,c,sl] = temp
      close, 1
   endfor
endfor





;;================ DARK SUBTRACTION
if (size(findfile('darkmap.sav')))[0] ne 0 then $
   restore, 'darkmap.sav' else print, '*** error: darkmap.sav file not found ***'
darkspec = make_array(3648,2,n_elements(texp))
darkblack = make_array(13,2,n_elements(texp))
for sl=0, n_elements(texp)-1 do begin
   for c=0, 1 do begin
      for p=0, 13-1 do darkblack[p,c,sl] = poly2d(texp[sl],temperature[c,sl],darkmodblack[*,*,p,c])
      for p=0, 3648-1 do begin
         darkspec[p,c,sl] = poly2d(texp[sl],temperature[c,sl],darkmodspec[*,*,p,c])
         spec[p,c,sl] = (spec[p,c,sl]-darkspec[p,c,sl]) -mean(black[*,c,sl]-darkblack[*,c,sl])
      endfor
   endfor
endfor




;;================ WAVELENGTH CALIBRATION
coefs = [ [355.688, 0.167436, -2.93242e-06, -2.22549e-10], $
          [360.071, 0.165454, -3.35036e-06, -1.88750e-10] ]
wavs = make_array(2, 3648)
for c=0, 1 do wavs[c,*] = poly(findgen(3648),coefs[*,c])
plot,wavs
for sl=0, n_elements(texp)-1 do $
   spec[*,0,sl] = interpol(spec[*,0,sl],wavs[0,*],wavs[1,*])

stop
END


;;================ "FLAT FIELDING" / DIFFERENTIAL TRANSMISSION CORRECTION
if (size(findfile('transmission.sav')))[0] ne 0 then $
   restore, 'transmission.sav' else print, '*** error: transmission.sav file not found ***'
for sl=0, n_elements(texp)-1 do spec[*,1,sl] /= T2
; T2=spec[*,1,33]/spec[*,0,33]
; dir='Cabauw2013/08072013/ScriptRun_08072013 151130'


;;================ EFFICIENCY CORRECTION INIT
; dir='Cabauw2013/05072013/POL_1a'
if (size(findfile('efficiency.sav')))[0] ne 0 then $
   restore, 'efficiency.sav' else print, '*** error: efficiency.sav file not found ***'
fitoutpol = fitout
; correct sheet polarizer performance above 672 nm
wl660 = (sort(abs(fitoutpol[*,5,0]-660)))[0]
wl672 = (sort(abs(fitoutpol[*,5,0]-672)))[0]
for sl=1, (size(fitoutpol))[3]-1 do $
   fitoutpol[wl672:*,1,sl] = mean(fitoutpol[wl660:wl672,1,sl])
; correct sheet polarizer aolp (determined) above 672 nm
wl400 = (sort(abs(fitoutpol[*,5,0]-400)))[0]
for sl=1, (size(fitoutpol))[3]-1 do $
   fitoutpol[wl672:*,3,sl] = mean(fitoutpol[wl400:wl672,3,sl])
; duplicate aolp right above 0 to right above 180, and right below 180 to right below 0
; to allow for aolp interpolation close to 0 and 180
aolpwrap = transpose( $
           [ transpose([(fitoutpol[*,3,1:18] mod !dpi +!dpi) mod !dpi]), $
             reform([(fitoutpol[*,3,7] mod !dpi +!dpi) mod !dpi -!dpi],1,1,3613), $
             reform([(fitoutpol[*,3,8] mod !dpi +!dpi) mod !dpi +!dpi],1,1,3613) ] /!dtor)
fitoutwrap = transpose( $
           [ transpose(fitoutpol[*,*,1:18]), $
             reform(transpose(fitoutpol[*,*,7]),1,8,3613), $
             reform(transpose(fitoutpol[*,*,8]),1,8,3613) ] )


;;================ DEMODULATE POLARIZATION
; dir='Cabauw2013/08072013/ScriptRun_08072013 151130'
; polspec=spec[*,*,90]
if (size(findfile('polspec.sav')))[0] ne 0 then $
   restore, 'polspec.sav' else print, '*** error: polspec.sav file not found ***'
inp = {comp:['mgf2','sio2'], thick:[3.82,-1.63]}
wavsDM = reform(wavs[1,*])
specDM = transpose([ reform(transpose(polspec),1,2,3648), transpose(spec) ])
demod, inp, wavsDM, specDM, fitout, $
       doplot=1, aolp=90*!dtor, lrange=[370., 865.];, /quiet

;;================ EFFICIENCY CORRECTION APPLY
fitoutcal = fitout
eff = make_array( (size(fitout))[1], n_elements(texp)+1 )
for sl=0, n_elements(texp) do begin
   for i=0, (size(fitout))[1]-1 do begin
      aolp = (fitout[i,3,sl]/!dtor mod 180 +180) mod 180
      aolpels = (sort(abs(aolpwrap[i,0,*]-aolp)))[0:1] ; linear interpolation
        eff[i,sl] =  ( (aolp-aolpwrap[i,0,aolpels[0]]) *fitoutwrap[i,1,aolpels[1]] - $
                       (aolp-aolpwrap[i,0,aolpels[1]]) *fitoutwrap[i,1,aolpels[0]] ) / $
                     (aolpwrap[i,0,aolpels[1]]-aolpwrap[i,0,aolpels[0]])
   endfor
   fitoutcal[*,1,sl] /= eff[*,sl]
endfor


;;================ SMOOTH POLARIZATION
k = fitoutcal[*,7,0]/fitoutcal[*,5,0]
for i=0, (size(fitoutcal))[1]-1 do begin
   winmin = min(where(abs(k-(k[i]+0.5)) eq min(abs(k-(k[i]+0.5)))))
   winmax = max(where(abs(k-(k[i]-0.5)) eq min(abs(k-(k[i]-0.5)))))
   for sl=0, n_elements(texp) do begin
      fitoutcal[i,1,sl] = mean(fitoutcal[winmin:winmax,1,sl])
      fitoutcal[i,3,sl] = mean(fitoutcal[winmin:winmax,3,sl])
   endfor
endfor


;;================ POINTING
; zero-point calibration
Dtilt = 4.00
Dpan = 203.87
; viewing angles
vzavazi = make_array(2,n_elements(texp))
for sl=0, n_elements(texp)-1 do begin
   vzavazi[0,*] = 90-(pantilt[1,*]-Dtilt)
   vzavazi[1,*] = (pantilt[0,*]-Dpan +360) mod 360
endfor
; solar angles
szasazi = make_array(2,n_elements(texp))
for sl=0, n_elements(texp)-1 do begin
   temp = strsplit(specfiles[0,sl],'_',/extract)
   julian = julday( date[1], date[0], date[2], $
                    strmid(temp[n_elements(temp)-2],9,2) -2d, $ ; we are GMT+2 in summer
                    strmid(temp[n_elements(temp)-2],11,2), $
                    strmid(temp[n_elements(temp)-2],13,2) )
   sunpos, julian, ra, dec
   eq2hor, ra, dec, julian, alt, azi, lat=51.968298, lon=4.927692
   szasazi[*,sl] = [90.-alt, azi]
endfor
; scattering angles
razi = (vzavazi[1,*] - szasazi[1,*] +360) mod 360
scat = make_array(n_elements(texp))
for sl=0, n_elements(texp)-1 do begin
   scat[sl] = acos( [sin(szasazi[0,sl]*!dtor), 0, cos(szasazi[0,sl]*!dtor)] ## $
         transpose([sin(vzavazi[0,sl]*!dtor)*cos(razi[sl]*!dtor), $
         sin(vzavazi[0,sl]*!dtor)*cos(razi[sl]*!dtor), cos(vzavazi[0,sl]*!dtor)]) ) /!dtor
   if (abs(razi[sl]-180) gt 90) and (vzavazi[0,sl] gt szasazi[0,sl]) then scat[sl] *= -1
endfor


;;================ RETRIEVAL INPUT @ AERONET WAVELENGTHS INIT
wl = [ [441.0, 10, 0.0428466], $
       [675.0, 10, 0.0228766], $
       [870.2, 10, 0.1101490] ]    ; AERONET wl, filter width, gain
scat[n_elements(texp)/2:*] = abs(scat[n_elements(texp)/2:*]) ; fix signs (ad hoc)
scat6 = where(scat ge 6)
dat = make_array(4, n_elements(scat6), (size(wl))[2])
for k=0, (size(wl))[2]-1 do begin
   win = [ (sort(abs(fitoutcal[*,5,0]-(wl[0,k]-wl[1,k]/2))))[0], $
           (sort(abs(fitoutcal[*,5,0]-(wl[0,k]+wl[1,k]/2))))[0] ]
   dat[0,*,k] = vzavazi[0,scat6]
   dat[1,*,k] = razi[scat6]
   for j=0, n_elements(scat6)-1 do begin
      dat[2,j,k] = $
      mean( specdm[win[0]:win[1],0,scat6[j]+1] + $
            specdm[win[0]:win[1],1,scat6[j]+1] /$
            fitoutcal[win[0]:win[1],6,scat6[j]+1] ,/double) / texp[scat6[j]] *wl[2,k] $
      *0.01 / ASTM(wl[0,k]-wl[1,k]/2,wl[0,k]+wl[1,k]/2) ; normalized to ASTM standard solar spectrum
      dat[3,j,k] = mean( fitoutcal[win[0]:win[1],1,scat6[j]+1], /double )
   endfor
endfor


;;================ EXPORT
temp = strsplit(specfiles[0,0],'_',/extract)
tstamp = strmid(temp[n_elements(temp)-2],0,8) + '_' + $
         strmid(temp[n_elements(temp)-2],9,6)
; IDL .sav: complete output
save, dir, date, times, texp, scat, pantilt, Dpan, Dtilt, $
      fitoutcal, specdm, eff, fitout, filename=dir+'SPEX_output_'+tstamp+'.sav'
; .dat: output
openw, 1, dir+'SPEX_output_'+tstamp+'.dat'
printf, 1, '! directory', dir, $
           '! date [day, month, year]', date, $
           '! times [hour (local time)]', times, $
           '! exposure times [ms]', texp, $
           '! scattering angles [deg]', scat, $
           '! PTU coords [pan, tilt]', pantilt, $
           '! PTU offset [Dpan, Dtilt]', [Dpan, Dtilt], $
           '! fitoutcal [offset, DoLP (calibrated and smoothed), slope(DoLP), AoLP [rad] (smoothed), chi^2(best fit), wavelength [nm], T1/T0, MOR [mm]]', fitoutcal, $
           '! spectra [counts]', specdm, $
           '! polarimetric efficiency', eff, $
           '! fitout [offset, DoLP (raw), slope(DoLP), AoLP [rad] (raw), chi^2(best fit), wavelength [nm], T1/T0, MOR [mm]]', fitout
close, 1
; .dat: retrieval input @ AERONET wavelengths
openw, 1, dir+'AERO_input_'+tstamp+'.dat'
printf, 1, 'Number of viewing zenith angles'
printf, 1, n_elements(scat6)
printf, 1, 'Number of wavelengths'
printf, 1, (size(wl))[2]
printf, 1, 'Wavelengths'
printf, 1, transpose(wl[0,*])
printf, 1, 'Solar zenith angle'
printf, 1, mean(szasazi[0,scat6],/double)
printf, 1, ['          VZA','          AZ','          I','          P']
for k=0, (size(wl))[2]-1 do printf, 1, dat[*,*,k]
close, 1


stop
END



