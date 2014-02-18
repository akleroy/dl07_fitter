function merge_draine_models $
   , model_1 $
   , model_2 $
   , weight_1 = weight_1 $
   , quiet=quiet $
   , nospec=nospec

;+
;
; Merge two "Draine model" structures using linear interpolation. Give
; a weight of "weight_1" to model_1 and a weight of "1-weight_1" to
; the second one. Return a new, hybridized model.
;
; This program is mostly used by other programs to build grids but
; could have more general applications.
;
;-

; SET THE DEFAULT WEIGHT
if n_elements(weight_1) eq 0 then begin
    message, 'Defaulting to equal weights.', /info
    weight_1 = 0.5
endif

; ERROR TRAP
if weight_1 lt 0.0 or weight_1 gt 1.0 then begin
    message, 'WARNING! Weight outside 0-1, probably an error.', /info
endif

; CALCULATE THE WEIGHT FOR THE SECOND MODEL
weight_2 = 1.0 - weight_1

; INITIALIZE OUTPUT
new_model = new_draine_model(nospec=nospec)

; HEADER MATERIAL
new_model.filename = 'merge'
new_model.directory = ''
new_model.model_num = -1

; GRAIN MODEL
if (model_1.grain_model ne model_2.grain_model) then begin
    if keyword_set(quiet) eq 0 then begin
        message, 'Warning: attempting to merge different grain models.', /info
    endif
    new_model.grain_model = 'hybrid'
endif else begin
    new_model.grain_model = model_1.grain_model
endelse

; MODEL PARAMETERS
new_model.qpah = $
  weight_1*(model_1.qpah) + weight_2*(model_2.qpah)

new_model.umin = $
  weight_1*(model_1.umin) + weight_2*(model_2.umin)

new_model.umax = $
  weight_1*(model_1.umax) + weight_2*(model_2.umax)

new_model.beta = $
  weight_1*(model_1.beta) + weight_2*(model_2.beta)

new_model.power_per_hatom = $
  weight_1*(model_1.power_per_hatom) + weight_2*(model_2.power_per_hatom)

new_model.mdust_per_mh = $
  weight_1*(model_1.mdust_per_mh) + weight_2*(model_2.mdust_per_mh)

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; LINEARLY AVERAGE BANDS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

new_model.band_vals = $
  weight_1*(model_1.band_vals) + weight_2*(model_2.band_vals)

; IF WE ARE KEEPING THE SPECTRUM ALSO AVERAGE THAT
if keyword_set(nospec) eq 0 then begin
   new_model.int = $
      weight_1*(model_1.int) + weight_2*(model_2.int)
   
   new_model.lam = $
      model_1.lam   
endif

return, new_model

end

