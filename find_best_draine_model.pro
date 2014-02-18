function find_best_draine_model $
;    X, Y, ERROR
   , bands=data_bands $
   , data=data_vec_in $
   , error=error_vec_in $
;   THE MODEL GRID (WHEN RERUNNING, USE THESE PARAMETERS TO SAVE TIME)
   , mw_umin = mw_umin $
   , mw_gamma = mw_gamma $
   , mw_qpah = mw_qpah $
   , mw_cube = mw_cube $
   , load_models = load_models $
;  THROW THIS FLAG IF YOU HAVE SUBMM DATA TO USE A WIDER RANGE OF U-FIELDS
   , submm = submm $
;   SUPPRESS INFORMATIONAL MESSAGES
   , quiet = quiet

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS AND DEFINITIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  nan = !values.d_nan

; DEFAULT TO THE FULL SET OF MODELS
  if keyword_set(load_models) or $
     n_elements(mw_gamma) eq 0 or $
     n_elements(mw_umin) eq 0 or $
     n_elements(mw_qpah) eq 0 or $
     n_elements(mw_cube) eq 0 then begin
     if keyword_set(quiet) eq 0 then $
        message, 'Reloading MW data from disk.', /info
     restore, '$MY_IDL/draine/draine_mw_grid.idl'
; IF WE LACK SUBMM DATA, EXCISE MINIMUM RADIATION FIELDS BELOW 0.7
     if keyword_set(submm) eq 0 then begin
        mw_cube = mw_cube[6:*,*,*]
        mw_umin = mw_umin[6:*]
     endif     
  endif

; GET SIZE OF LOOP/CUBE
  n_umin = n_elements(mw_umin)
  n_qpah = n_elements(mw_qpah)
  n_gamma = n_elements(mw_gamma)

; INIALIZE OUTPUT ARRAYS
  gamma_ra = fltarr(n_umin, n_qpah, n_gamma)*!values.f_nan
  scale_ra = fltarr(n_umin, n_qpah, n_gamma)*!values.f_nan
  chisq_ra = fltarr(n_umin, n_qpah, n_gamma)*!values.f_nan
  red_chisq_ra = fltarr(n_umin, n_qpah, n_gamma)*!values.f_nan

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DATA AND ERROR CHECKING
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; CATCH NON-FINITE DATA AND CRY FOUL
  if total(finite(data_vec_in) eq 0) gt 0 then begin
     message, 'Cannot fit non-finite data. Returning empty output.',/info

     results = {$
               bands: data_bands $
               ,data: double(data_vec_in) $
               ,unc: nan*data_vec_in $
               ,fit: double(nan*data_vec_in) $
               ,model_lam: mw_cube[0].lam $
               ,model_int: mw_cube[0].int $
               ,chisq: nan, red_chisq: nan $
               ,qpah: nan $
               ,min_qpah: nan, max_qpah: nan, exp_qpah:nan $
               ,umin: nan $
               ,min_umin: nan,max_umin: nan, exp_umin:nan $
               ,gamma: nan $
               , min_gamma:nan, max_gamma: nan, exp_gamma:nan $
               ,mdust_per_mh: nan $
               ,min_mdust_per_mh: nan,max_mdust_per_mh: nan $
               ,exp_mdust_per_mh: nan $
               ,dust_sd: nan $
               ,min_dust_sd: nan,max_dust_sd: nan $
               ,exp_dust_sd: nan $
               }
     
     return, results
  endif
  
; COPY... JUST IN CASE
  data_vec = data_vec_in

; CHECK THAT WE HAVE SOME DATA
  ndata = n_elements(data_vec)
  if ndata eq 0 then begin 
     message, 'Require data to fit! Returning.', /info
     return, -1
  endif

; INITIALIZE THE ARRAY TO CROSS-REFERENCE THE MODEL
  band_inds = lonarr(ndata)
  for i = 0, ndata-1 do begin
     band_inds[i] = where(mw_cube[0].bands eq data_bands[i],ct)
     if ct ne 1 then begin
        message, 'Cross referencing problem. Stopping.', /info
        stop
     endif    
  endfor

; MAKE A DEFAULT ERROR VECTOR
  if n_elements(error_vec_in) ne n_elements(data_vec) then begin
     message, 'Bad or missing error vector. Using 10% value as weight.', /info
     error_vec = data_vec*0.1
  endif else begin
;    ADD 10% IN QUADRATURE TO THE UNCERTAINTY
;    ... CLOSE, BUT NOT EXACTLY, THE DRAINE APPROACH (+10% OF MODEL)
     if keyword_set(quiet) eq 0 then $
        message, 'Adding 10% to uncertainties in quadrature.', /info
     error_vec = error_vec_in
  endelse

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; FIND THE BEST-FIT SCALING AND ASSOCIATED CHI SQUARED FOR EACH MODEL
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; AKL - IN THEORY, WE COULD INTERPOLATE BETWEEN ADJACENT MODELS TO IMPROVE THE
;       RESOLUTION. IN PRACTICE, WE'RE ALREADY GETTING KILLED ON PROCESSOR
;       CYCLES. WITHOUT A SMART SEARCH ALGORITHM OR A FASTER FITTER, LET'S
;       STICK TO THE MODELS WE HAVE.

  ndata = n_elements(data_vec)

  for i = 0, n_umin-1 do begin
     for j = 0, n_qpah-1 do begin
        for k = 0, n_gamma-1 do begin
           model_vec = (mw_cube[i,j,k].band_vals)[band_inds]

           b = dblarr(1,ndata)
           b[0,*] = model_vec / error_vec
           y = data_vec / error_vec
           psi = invert(transpose(b) ## b)
           result = (psi ## transpose(b) ## transpose(y))[0]
           chisq = total((data_vec - result*model_vec)^2/error_vec^2)
           
           scale_ra[i,j,k] = result
           chisq_ra[i,j,k] = chisq
           red_chisq_ra[i,j,k] = chisq / (n_elements(model_vec) - 1)
           gamma_ra[i,j,k] = mw_gamma[k]
        endfor
     endfor
  endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; WORK OUT KEY RESULTS FROM THE CHI-SQUARED CUBE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; THE LOWEST CHISQUARED VALUE (AND ITS LOCATION)
  chisq = min(chisq_ra, minind)
  red_chisq = chisq / n_elements(data_vec - 3.)

; THE SURFACE ENCOMPASSING CHI-SQUARED = 1
  chisq1_ind = where(chisq_ra le (chisq+1.), chisq1_ct)

; THE BEST-FIT VALUES
  best_scale = scale_ra[minind]
  fit_vals = (mw_cube[minind].band_vals)[band_inds]*best_scale

; DUST SURFACE DENSITIES
  mh = 1.6735340e-24            ; HYDROGEN MASS
  dust_sd =  $
     1e-23 * scale_ra * 1e6 $    ; FLUX PER STER IN ERG/S/CM^-2/HZ
     * 4. * !pi $                ; FLUX -> LUMINOSITY
     / (1e-23*4.*!pi) $          ; POWER PER H
     * mh * (mw_cube.mdust_per_mh) ; DUST/H IN THE MODEL

; NOW WORK OUT EXPECTATION VALUES APPLYING EQUAL WEIGHTING ACROSS THE WHOLE
; CHI-SQUARED = 1 SPACE
  exp_qpah = $
     10.^(total(alog10(mw_cube[chisq1_ind].qpah))/(1.0d*chisq1_ct))
  exp_umin = $
     10.^(total(alog10(mw_cube[chisq1_ind].umin))/(1.0d*chisq1_ct))
  exp_gamma = $
     10.^(total(alog10(gamma_ra[chisq1_ind]))/(1.0d*chisq1_ct))
  exp_mdust_per_mh = $
     10.^(total(alog10(mw_cube[chisq1_ind].mdust_per_mh))/(1.0d*chisq1_ct))
  exp_dust_sd = $
     10.^(total(alog10(dust_sd[chisq1_ind]))/(1.0d*chisq1_ct))

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; COMPILE A FINAL "RESULTS" STRUCTURE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  results = {$
            bands: data_bands $
            ,data: double(data_vec) $
            ,unc: double(error_vec) $
            ,fit: double(fit_vals) $
            ,model_lam: mw_cube[minind].lam $
            ,model_int: mw_cube[minind].int*best_scale $
            ,chisq: double(chisq) $
            ,red_chisq: double(red_chisq) $              
            ,qpah: double(mw_cube[minind].qpah) $
            ,min_qpah: double(min(mw_cube[chisq1_ind].qpah)) $
            ,max_qpah: double(max(mw_cube[chisq1_ind].qpah)) $
            ,exp_qpah: double(exp_qpah) $
            ,umin: double(mw_cube[minind].umin) $
            ,min_umin: double(min(mw_cube[chisq1_ind].umin)) $
            ,max_umin: double(max(mw_cube[chisq1_ind].umin)) $
            ,exp_umin: double(exp_umin) $
            ,gamma: double(gamma_ra[minind]) $
            ,min_gamma: double(min(gamma_ra[chisq1_ind])) $
            ,max_gamma: double(max(gamma_ra[chisq1_ind])) $
            ,exp_gamma: double(exp_gamma) $
            ,mdust_per_mh: double(mw_cube[minind].mdust_per_mh) $
            ,min_mdust_per_mh: double(min(mw_cube[chisq1_ind].mdust_per_mh)) $
            ,max_mdust_per_mh: double(max(mw_cube[chisq1_ind].mdust_per_mh)) $
            ,exp_mdust_per_mh: double(exp_mdust_per_mh) $
            ,dust_sd: double(dust_sd[minind]) $
            ,min_dust_sd: double(min(dust_sd[chisq1_ind])) $
            ,max_dust_sd: double(max(dust_sd[chisq1_ind])) $
            ,exp_dust_sd: double(exp_dust_sd) $
            }

  return, results

end
