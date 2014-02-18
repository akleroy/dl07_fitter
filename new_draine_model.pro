function new_draine_model $
   , nospec=nospec

;+
; NAME:
;
; new_draine_model
;
; PURPOSE:
;
; Create an empty IDL structure to hold a Draine & Li (2007) model.
;
; CATEGORY:
;
; Part of my DL07 fitting package.
;
; CALLING SEQUENCE:
;
; struct = new_draine_model()
;
; MODIFICATION HISTORY:
;
; written - summer 2009
; documented - winter 2010 aleroy@nrao.edu
; added /nospec - feb 14
; 
;-

nan = !values.f_nan

bands = $
  [ 'dirbe_1.3' $
    , 'dirbe_2.2' $
    , 'dirbe_3.5' $
    , 'dirbe_4.9' $
    , 'dirbe_12' $
    , 'dirbe_21' $
    , 'dirbe_56' $
    , 'dirbe_98' $
    , 'dirbe_148' $
    , 'dirbe_248' $
    , 'irac_3.6' $
    , 'irac_4.5' $
    , 'irac_5.7' $
    , 'irac_7.9' $
    , 'mips_24' $
    , 'mips_70' $
    , 'mips_160' $
    , 'pacs_70' $
    , 'pacs_100' $
    , 'pacs_160' $
    , 'spire_250' $
    , 'spire_350' $
    , 'spire_500' $
    , 'iras_12' $
    , 'iras_25' $
    , 'iras_60' $
    , 'iras_100' $
  ]

nan_ra = dblarr(n_elements(bands))*!values.f_nan

if keyword_set(nospec) then begin
   empty_model = { $
                 filename: '', $
                 directory: '', $
                 model_num: -1, $
                 grain_model: '', $
                 qpah: nan, $
                 umin: nan, $
                 umax: nan, $
                 beta: nan, $
                 avg_u: nan, $
                 mdust_per_mh: nan, $
                 power_per_hatom: nan, $
                 bands: bands, $
                 band_vals: nan_ra $
                 }
endif else begin
   empty_model = { $
                 filename: '', $
                 directory: '', $
                 model_num: -1, $
                 grain_model: '', $
                 qpah: nan, $
                 umin: nan, $
                 umax: nan, $
                 beta: nan, $
                 avg_u: nan, $
                 mdust_per_mh: nan, $
                 power_per_hatom: nan, $
                 bands: bands, $
                 band_vals: nan_ra, $
                 lam: fltarr(1001)*nan, $
                 int: fltarr(1001)*nan $
                 }
endelse

return, empty_model

end
