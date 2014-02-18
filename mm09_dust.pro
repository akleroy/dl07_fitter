function mm09_dust $
   , ir24 = ir24 $
   , ir70 = ir70 $
   , ir160 = ir160 $  
   , color = color $
   , units = units $
   , avg_u = avg_u

  @constants.bat
  nu24 = c/(24.d-4)
  nu70 = c/(70.d-4)
  nu160 = c/(160.d-4)

  units = 'Msun/pc2'
  
  if n_elements(color) eq 0 then $
     color = ir70/ir160

; ESTIMATE THE AVERAGE RADIATION FIELD
  r70 = (nu70*ir70) / (nu160*ir160)
  avg_u = 10.^(0.468 + 1.801*alog10(r70))

; ESTIMATE F_PDR (GAMMA)
 
; ... tbd

; ESTIMATE THE DUST SURFACE DENSITY
  sigma_dust = $
     4.*!pi/1.616d-13*(1.00)^2 * $
     (nu70/nu160*color)^(-1.801) * $
     (1.559*(nu24*ir24) + $
      0.7686*(nu70*ir70) + $
      1.347*(nu160*ir160)) * $
     (1d6*1d-23) / $
     (1000.d3)^2
  
  return, sigma_dust

end
   
