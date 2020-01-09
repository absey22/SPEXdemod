
;function that SPEXredux.pro calls to make a 2d polynomial?


function poly2d, x, y, coefs

;print, 'x',x
;print,'y',y
;print, 'coefs',coefs


z = make_array(n_elements(x), n_elements(y), /double)
  
for xx=0, n_elements(x)-1 do begin
   for yy=0, n_elements(y)-1 do begin
      for i=0, (size(coefs))[2]-1 do begin
         for j=0, (size(coefs))[1]-1 do begin
            z[xx,yy] += coefs[j,i] *x[xx]^i *y[yy]^j
         endfor
      endfor
   endfor
endfor


return, z
end
