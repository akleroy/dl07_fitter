function find_best_draine_model $
   , bands=data_bands $
   , data=data_vec_in $
   , error=error_vec_in $
   , umin_vec = umin_vec $
   , gamma_vec = gamma_vec $
   , qpah_vec = qpah_vec $
   , model_cube = model_cube $
   , load_models = load_models $
   , gal = gal $
   , gof_cube = gof_cube $
   , add_error = add_error $
   , submm = submm $
   , quiet = quiet

;+
; NAME:
;
; find_best_draine_model
;
; PURPOSE:
;
; Identify the best-fit Draine & Li (2007) model given a set of
; measurements and errors in broad bands from space telescopes.
;
; CATEGORY:
;
; Science!
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
; Revised and documented feb 14 - aleroy
;
;-

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS AND DEFINITIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; CONVENIENCE
  nan = !values.d_nan

; COPY INPUT VECTOR
  data_vec = data_vec_in  

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; LOAD MODEL DATA
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; DEFAULT MILKY WAY MODEL
  if n_elements(gal) eq 0 then $
     gal = "MW"

; IDENTIFY IDL FILE TO LOAD
  if gal eq "MW" then begin
     model_file = 'draine_grid_MW.idl'
  endif else if gal eq "LMC" then begin
     model_file = 'draine_grid_LMC.idl'
  endif else if gal eq "SMC" then begin
     model_file = 'draine_grid_SMC.idl'
  endif else begin
     message, "Invalid galaxy string."
  endelse

; RELOAD THE MODELS IF NECESSARY OR REQUESTED
  if keyword_set(load_models) or $
     n_elements(gamma_vec) eq 0 or $
     n_elements(umin_vec) eq 0 or $
     n_elements(qpah_vec) eq 0 or $
     n_elements(model_cube) eq 0 then begin
     
     if keyword_set(quiet) eq 0 then $
        message, 'Reloading MW data from disk.', /info

     if n_elements(model_dir) eq 0 then begin
        model_dir = '$MY_IDL/dl07_fitter/'
     endif
     restore, model_dir+model_file

; IF WE LACK SUBMM DATA, EXCISE MODELS WITH MINIMUM RADIATION FIELDS
; BELOW 0.7 FROM THE GRID
     if keyword_set(submm) eq 0 then begin
        umin_lt_0p7 = where(umin_vec lt 0.7, ct_low)
        if ct_low gt 0 then begin
           min_umin_ind = (max(umin_lt_0p7)+1)
           model_cube = model_cube[min_umin_ind:*,*,*]
           umin_vec = umin_vec[min_umin_ind:*]
        endif     
     endif
  endif
  
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; INITIALIZE OUTPUT
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; GET SIZE OF LOOP/CUBE
  n_umin = n_elements(umin_vec)
  n_qpah = n_elements(qpah_vec)
  n_gamma = n_elements(gamma_vec)

; INIALIZE OUTPUT ARRAYS

; ... GAMMA PARAMETER
  gamma_cube = dblarr(n_umin, n_qpah, n_gamma)*!values.f_nan

; ... SCALING VALUE TO GET FROM MODEL TO DATA (RELATED TO MASS)
  scale_cube = dblarr(n_umin, n_qpah, n_gamma)*!values.f_nan

; ... GOODNESS OF FIT PARAMETER (CAN BE TUNED)
  gof_cube = dblarr(n_umin, n_qpah, n_gamma)*!values.f_nan
    
; INITIALIZE OUTPUT

; ... STRUCTURE HOLDS ONE PARAMETER
  fparm = {$
          best: nan $
          , min: nan $
          , max: nan $
          , mean: nan $
          }

; ... STRUCTURE HOLDS WHOLE FIT
  results = {$
            bands: data_bands $
            , data: double(data_vec) $
            , unc: nan*data_vec $
            , fit: double(nan*data_vec) $
            , gal: gal $
            , gof: nan $
            , tol_gof: nan $
            , qpah: fparm $
            , umin: fparm $
            , gamma: fparm $
            , mdust_per_mh: fparm $
            , dust_sd: fparm $
            }
  
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ERROR CHECK
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; CHECK THAT WE HAVE SOME DATA
  ndata = n_elements(data_vec)
  if ndata eq 0 then begin 
     message, 'Require data to fit! Returning.', /info
     return, -1
  endif

; CATCH NON-FINITE DATA AND CRY FOUL
  if total(finite(data_vec) eq 0) gt 0 then begin
     message, 'Cannot fit non-finite data. Returning empty output.',/info     
     return, results
  endif
  
; FIGURE OUT WHAT BAND INDICES IN THE MODEL CORRESPOND TO THE INPUT
  band_inds = lonarr(ndata)
  for i = 0, ndata-1 do begin
     band_inds[i] = where(model_cube[0].bands eq data_bands[i],ct)
     if ct ne 1 then begin
        message, 'Band identification problem. Stopping.', /info
        stop
     endif    
  endfor

; MAKE A DEFAULT ERROR VECTOR IF NONE ARE SUPPLIED
  if n_elements(error_vec_in) ne n_elements(data_vec) then begin
     message, 'Bad or missing error vector. Using 10% value as weight.', /info
     error_vec = data_vec*0.1
  endif else begin
     error_vec = error_vec_in  
  endelse

; ADD 10% IN QUADRATURE TO THE UNCERTAINTY ... CLOSE, BUT NOT
; EXACTLY, THE DRAINE APPROACH (+10% OF MODEL)
  if keyword_set(add_error) then begin
     if keyword_set(quiet) eq 0 then $
        message, 'Adding 10% to uncertainties in quadrature.', /info
     error_vec = sqrt(error_vec^2 + (0.1*data_vec)^2)
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; FIND THE BEST-FIT SCALING GOODNESS OF FIT FOR EACH MODEL
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; N.B. We could interpolate between adjacent models to improve the
; resolution, but the fitting is already slow and uncertain. The need
; for this is not clear.

; DEFINE DEGREES OF FREEDOM
  dof = 0

; LOOP OVER ALL MODELS
  for i = 0, n_umin-1 do begin
     for j = 0, n_qpah-1 do begin
        for k = 0, n_gamma-1 do begin

;          ... EXTRACT AN INDIVIDUAL MODEL
           model_vec = (model_cube[i,j,k].band_vals)[band_inds]

;          ... DO A LINEAR LEAST SQUARES TO FIND THE BEST SCALING OF
;          THIS MODEL TO THE DATA
           b = dblarr(1,ndata)
           b[0,*] = model_vec / error_vec
           y = data_vec / error_vec
           psi = invert(transpose(b) ## b)
           result = (psi ## transpose(b) ## transpose(y))[0]

;          ... CALCULATE THE GOODNESS OF FIT
           chisq = total((data_vec - result*model_vec)^2/error_vec^2)
           
;          ... SAVE THE RESULTS
           scale_cube[i,j,k] = result
           gamma_cube[i,j,k] = gamma_vec[k]
           gof_cube[i,j,k] = chisq / (n_elements(model_vec) - 1 - dof)

        endfor
     endfor
  endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ANALYZE THE GOODNESS OF FIT CUBE TO FIND THE BST 
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; THE LOWEST GOODNESS OF FIT VALUE (AND ITS LOCATION)
  best_gof = min(gof_cube, best_ind)

; DEFINE TOLERANCE FOR ERROR DEFINITION
  if n_elements(tol) eq 0 then $
     tol = 0.1*best_gof

; THE SURFACE OF INTEREST FOR DEFINING ERRORS
  within_tol = where(gof_cube le (gof_cube+tol), tol_ct)

; THE BEST-FIT VALUES
  best_scale = scale_cube[best_ind]
  fit_vals = (model_cube[best_ind].band_vals)[band_inds]*best_scale

; DUST SURFACE DENSITIES
  mh = 1.6735340e-24            ; HYDROGEN MASS

  dust_sd_cube =  $
     1e-23 * scale_cube * 1e6 $  ; FLUX PER STER IN ERG/S/CM^-2/HZ
     * 4. * !pi $                ; FLUX -> LUMINOSITY
     / (1e-23*4.*!pi) $          ; POWER PER H
     * mh $                      ; MASS OF HYDROGEN
     * (model_cube.mdust_per_mh) ; DUST/H IN THE MODEL
  
; NOW WORK OUT MEAN (LOG) VALUES OVER THE "WITHIN TOLERANCE" SPACE
  exp_qpah = $
     10.^(total(alog10(model_cube[within_tol].qpah))/(1.0d*tol_ct))
  exp_umin = $
     10.^(total(alog10(model_cube[within_tol].umin))/(1.0d*tol_ct))
  exp_gamma = $
     10.^(total(alog10(gamma_cube[within_tol]))/(1.0d*tol_ct))
  exp_mdust_per_mh = $
     10.^(total(alog10(model_cube[within_tol].mdust_per_mh))/(1.0d*tol_ct))
  exp_dust_sd = $
     10.^(total(alog10(dust_sd_cube[within_tol]))/(1.0d*tol_ct))
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; PLACE RESULTS INTO A STRUCTURE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; THE DATA
  results.bands = data_bands
  results.data = double(data_vec)
  results.unc = double(error_vec)

; THE FIT
  results.gof = double(best_gof)
  results.tol_gof = double(tol)
  results.fit = double(fit_vals)  

; FIT PARAMETERS:

; ... QPAH
  results.qpah.best = double(model_cube[best_ind].qpah)
  results.qpah.min = double(min(model_cube[within_tol].qpah))
  results.qpah.max = double(max(model_cube[within_tol].qpah))
  results.qpah.mean = double(exp_qpah)

; ... UMIN
  results.umin.best = double(model_cube[best_ind].umin)
  results.umin.min = double(min(model_cube[within_tol].umin))
  results.umin.max = double(max(model_cube[within_tol].umin))
  results.umin.mean = double(exp_umin)

; ... GAMMA
  results.gamma.best = double(gamma_cube[best_ind])
  results.gamma.min = double(min(gamma_cube[within_tol]))
  results.gamma.max = double(max(gamma_cube[within_tol]))
  results.gamma.mean = double(exp_gamma)

; ... MDUST PER MH
  results.mdust_per_mh.best = double(model_cube[best_ind].mdust_per_mh)
  results.mdust_per_mh.min = double(min(model_cube[within_tol].mdust_per_mh))
  results.mdust_per_mh.max = double(max(model_cube[within_tol].mdust_per_mh))
  results.mdust_per_mh.mean = double(exp_mdust_per_mh)

; ... DUST SURFACE DENSITY
  results.dust_sd.best = double(dust_sd_cube[best_ind])
  results.dust_sd.min = double(min(dust_sd_cube[within_tol]))
  results.dust_sd.max = double(max(dust_sd_cube[within_tol]))
  results.dust_sd.mean = double(exp_dust_sd)

  return, results

end
