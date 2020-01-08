
;
FUNCTION ASTM, lambdamin, lambdamax

if (size(findfile('solarspec.sav')))[0] ne 0 then $
   restore, 'solarspec.sav' else print, '*** error: solarspec.sav file not found ***'

minel = (sort(abs(solar[*,0]-lambdamin)))[0]
maxel = (sort(abs(solar[*,0]-lambdamax)))[0]

irrad = mean(solar[minel:maxel,1])

RETURN, irrad
end
