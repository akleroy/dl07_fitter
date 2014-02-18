pro read_all_models $
   , basedir=basedir $
   , gal=gal

  if n_elements(basedir) eq 0 then $
     basedir = getenv('MY_IDL')+'/draine/irem4/'

  if n_elements(gal) eq 0 then $
     gal = "MW"

  if gal eq "MW" then begin
     outfile = "draine_models_MW.idl"
     file_list = $
        file_search(basedir+'U*','*MW*.txt', count=ct)
  endif else if gal eq "LMC" then begin
     outfile = "draine_models_LMC.idl"
     file_list = $
        file_search(basedir+'U*','*LMC*.txt', count=ct)
  endif else if gal eq "SMC" then begin
     outfile = "draine_models_SMC.idl"
     file_list = $
        file_search(basedir+'U*','*smc*.txt', count=ct)     
  endif else begin
     message, "Invalid galaxy string."
  endelse
  

  message, "Reading models for "+gal, /info

  for i = 0, ct-1 do begin
     lines = numlines(file_list[i])
     if lines lt 10 then begin
        message, file_list[i]+" appears invalid. Skipping.", /info
        continue
     endif
     model = read_one_model(file_list[i], /minmax_only, fail=fail)
     if fail then $
        continue
     if n_elements(draine_models) eq 0 then $
        draine_models = model $
     else $
        draine_models = [draine_models, model]
  endfor

  save, draine_models, file=outfile

end
