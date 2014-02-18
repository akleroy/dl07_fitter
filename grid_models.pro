pro grid_models, gal=gal, nospec=nospec

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD THE MODELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; DEFAULT MILKY WAY MODEL
  if n_elements(gal) eq 0 then $
     gal = "MW"

; RESTORE IDL FILE
  if gal eq "MW" then begin
     restore, "draine_models_MW.idl"
     outfile = 'draine_grid_MW.idl'
  endif else if gal eq "LMC" then begin
     restore, "draine_models_LMC.idl"
     outfile = 'draine_grid_LMC.idl'
  endif else if gal eq "SMC" then begin
     restore, "draine_models_SMC.idl"
     outfile = 'draine_grid_SMC.idl'
  endif else begin
     message, "Invalid galaxy string."
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SEPERATE THE SINGLE-U AND DISTRIBUTION OF U MODELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GRAB THE SINGLE-U DISTRIBUTIONS
  single_u = draine_models[where(draine_models.umin eq draine_models.umax)]

; GRAB THE POWER-LAW DISTRIBUTIONS...
  distrib_u = draine_models[where(draine_models.umin ne draine_models.umax)]
; ... AND DON'T WORRY ABOUT UMAX
  distrib_u = distrib_u[where(distrib_u.umax eq 1.00000e+06)]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONSTRUCT THE AXES - PAH FRACTION, RADIATION FIELD
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; IDENTIFY THE DIFFERENT VALUES OF QPAH IN THE ORIGINAL MODELS
  orig_qpah_vec = single_u.qpah
  ind = sort(orig_qpah_vec)
  orig_qpah_vec = orig_qpah_vec[ind]
  orig_qpah_vec = orig_qpah_vec[uniq(orig_qpah_vec)]

; MAKE A TARGET VECTOR OF PAH FRACTIONS THAT WE WILL INTERPOLATE TO
  min_qpah = min(orig_qpah_vec)
  max_qpah = max(orig_qpah_vec)
  n_qpah=41
  qpah_vec = findgen(n_qpah)*(max_qpah - min_qpah)/(n_qpah-1.0) + min_qpah

; IDENTIFY THE DIFFERENT VALUES OF MINIMUM U-FIELD IN THE MODELS
  umin_vec = single_u.umin
  ind = sort(umin_vec)
  umin_vec = umin_vec[ind]
  umin_vec = umin_vec[uniq(umin_vec)]
  umin_vec = umin_vec[where(umin_vec le max(distrib_u.umin))]
  n_umin = n_elements(umin_vec)

; MAKE A TARGET VECTOR USED TO MIX POWER LAW AND SINGLE U-FIELD MODELS
  n_gamma = 51
  gamma_step = 2.
  gamma_vec = ((findgen(n_gamma))*gamma_step)/100.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; INITIALIZE AN EMPTY CUBE OF MODELS
  dummy = new_draine_model(nospec=nospec)
  model_cube = replicate(dummy,n_umin,n_qpah,n_gamma)

; LOOP OVER U-FIELD AND QPAH
  for i = 0, n_umin-1 do begin    
     for j = 0, n_qpah-1 do begin        

;       AT EACH COMBINATION, INTERPOLATE TO THE DESIRED VALUE OF QPAH IN BOTH
;       ... THE SINGLE U-FIELD MODEL...
        up_ind = where(single_u.umin eq umin_vec[i] and $
                       single_u.qpah ge qpah_vec[j], ct)
        if ct eq 0 then stop
        high_qpah = min(single_u[up_ind].qpah,minind)
        single_qpah_above = single_u[up_ind[minind]]

        low_ind = where(single_u.umin eq umin_vec[i] and $
                        single_u.qpah le qpah_vec[j], ct)
        if ct eq 0 then stop
        low_qpah = max(single_u[low_ind].qpah,maxind)
        single_qpah_below = single_u[low_ind[maxind]]
        
;        print, high_qpah, low_qpah, qpah_vec[j]

        if (high_qpah eq low_qpah) then $
           weight_above = 0.5 $
        else $
           weight_above = $
           1.0 - (high_qpah - qpah_vec[j])/(high_qpah - low_qpah)

        umin_model = $
           merge_draine_models(single_qpah_above $
                               , single_qpah_below $
                               ,weight_1=weight_above $
                               , /quiet, nospec=nospec)

;       ... AND THE POWER LAW MODEL...
        up_ind = where(distrib_u.umin eq umin_vec[i] and $
                       distrib_u.qpah ge qpah_vec[j], ct)
        if ct eq 0 then begin
           message, 'Warning! Floating point comparisons causing problems.', /info
           dummy = max(distrib_u.qpah(where(distrib_u.umin eq umin_vec[i])),maxind)
           up_ind = [maxind]
        endif
        high_qpah = min(distrib_u[up_ind].qpah,minind)
        distrib_qpah_above = distrib_u[up_ind[minind]]

        low_ind = where(distrib_u.umin eq umin_vec[i] and $
                        distrib_u.qpah le qpah_vec[j], ct)
        if ct eq 0 then stop
        low_qpah = max(distrib_u[low_ind].qpah,maxind)
        distrib_qpah_below = distrib_u[low_ind[maxind]]

        if (high_qpah eq low_qpah) then $
           weight_above = 0.5 $
        else $
           weight_above = $
           1.0 - (high_qpah - qpah_vec[j])/(high_qpah - low_qpah)
        
        distrib_model = $
           merge_draine_models(distrib_qpah_above $
                               , distrib_qpah_below $
                               , weight_1=weight_above, /quiet)

        for k = 0, n_gamma-1 do begin
           dummy = $
              merge_draine_models(umin_model $
                                  ,distrib_model $
                                  ,weight_1=(1.0-gamma_vec[k]) $
                                  ,nospec=nospec)
           model_cube[i,j,k] = dummy
        endfor
     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT PROJECTIONS OF 8 MICRON EMISSION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 33

  irac8 = where(model_cube[0].bands eq 'irac_7.9')
  cube8 = (model_cube.band_vals[irac8])
  sz = size(cube8)
  cube8 = reform(cube8,sz[2],sz[3],sz[4])

  !p.multi=[0,2,2]
  disp, total(alog10(cube8),3)
  disp, total(alog10(cube8),2)
  disp, total(alog10(cube8),1)
  !p.multi=0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RENAME AXES AND SAVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  stop

  model_gal = gal
  save, umin_vec, gamma_vec, qpah_vec, model_cube, model_gal, file=outfile

end
