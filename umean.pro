function umean, umin=umin, umax=umax, gamma=gamma
  
  if n_elements(umax) eq 0 then $
     umax=1e6

; DRAINE+ 07 EQ 17
  umean = (1.0 - gamma)*umin + $
           gamma*(alog(umax/umin)/(umin^(-1.)+umax^(-1.)))

  return, umean

end
