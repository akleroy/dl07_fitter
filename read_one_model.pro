function read_one_model $
   , filename $
   , quiet=quiet $
   , minmax_only=minmax_only $
   , fail=fail

;+
; NAME:
;
; read_one_model
;
; PURPOSE:
;
; Read a single Draine & Li 2007 model from a text file and output a
; structure.
;
; CATEGORY:
;
; Part of my DL07 fitting package.
;
; CALLING SEQUENCE:
;
; struct = read_one_model(filename, quiet=quiet)
;
; MODIFICATION HISTORY:
;
; written - summer 2009
; documented - winter 2010 aleroy@nrao.edu
;
;-


; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
; WORK OUT DIRECTORIES
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$

  short_filename = filename
  directory=''
  slash = strpos(short_filename,'/')
  while (slash ne -1) do begin
      directory = directory + strmid(short_filename,0,slash+1)
      short_filename=$
        strmid(short_filename $
               , slash+1, strlen(short_filename)-(slash+1))
      slash = strpos(short_filename,'/')
  endwhile

  if keyword_set(quiet) eq 0 then begin
      message $
        , 'Reading file='+short_filename+' from directory='+directory $
        ,/info
  endif

; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
; INITIALIZE AN EMPTY STRUCTURE
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$

  model = new_draine_model()
  model.filename = short_filename
  model.directory = directory
  
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
; READ THE MODEL PROPERTIES
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$

  openr, 1, filename
  halt = 0
  while not(eof(1)) and (halt eq 0) do begin      
      line=''
      readf,1,line
      if strmid(line,0,6) eq ' band ' then halt = 1

;     LOOK UP THE NAME OF THE MODEL, WHICH MAPS TO THE PAH FRACTION      
      modelpos = strpos(line,'= grain model')
      if (modelpos ne -1) then begin
          modelstr = strmid(line,0,modelpos-1)  
          reads,modelstr,modelnum
          case modelnum of
              1.0: begin
                  model.model_num = 1
                  model.grain_model = 'MW3.1_00'
                  model.qpah = 0.47e-2
                  model.mdust_per_mh = 0.01
              end
              2.0: begin
                  model.model_num = 2
                  model.grain_model = 'MW3.1_10'
                  model.qpah = 1.12e-2
                  model.mdust_per_mh = 0.01
              end
              3.0: begin
                  model.model_num = 3
                  model.grain_model = 'MW3.1_20'
                  model.qpah = 1.77e-2
                  model.mdust_per_mh = 0.0101
              end
              4.0: begin
                  model.model_num = 4
                  model.grain_model = 'MW3.1_30'
                  model.qpah = 2.50e-2
                  model.mdust_per_mh = 0.0102
              end
              5.0: begin
                  model.model_num = 5
                  model.grain_model = 'MW3.1_40'
                  model.qpah = 3.19e-2
                  model.mdust_per_mh = 0.0102
              end
              6.0: begin
                  model.model_num = 6
                  model.grain_model = 'MW3.1_50'
                  model.qpah = 3.90e-2
                  model.mdust_per_mh = 0.0103
              end
              7.0: begin
                  model.model_num = 7
                  model.grain_model = 'MW3.1_60'
                  model.qpah = 4.58e-2
                  model.mdust_per_mh = 0.0104
              end
              8.0: begin
                  model.model_num = 8
                  model.grain_model = 'LMC2_00'
                  model.qpah = 0.75e-2
                  model.mdust_per_mh = 0.00343
              end
              9.0: begin
                  model.model_num = 9
                  model.grain_model = 'LMC2_05'
                  model.qpah = 1.49e-2
                  model.mdust_per_mh = 0.00344
              end
              10.0: begin
                  model.model_num = 10
                  model.grain_model = 'LMC2_10'
                  model.qpah = 2.37e-2
                  model.mdust_per_mh = 0.00359
              end
              11.0: begin
                  model.model_num = 11
                  model.grain_model = 'smc'
                  model.qpah = 0.10e-2
                  model.mdust_per_mh = 0.00206
              end
          endcase         
      endif

;     LOOK UP THE PARAMETERS OF THE RADIATION FIELD
      upos = strpos(line,'= Umin')
      if (upos ne -1) then begin
          ustr = strmid(line,0,upos-1)  
          reads, ustr, umin, umax, beta
          model.umin = umin
          model.umax = umax
          model.beta = beta
      endif

;     LOOK UP THE AVERAGE RADIATION FIELD
      uavgpos = strpos(line,'= <U>')
      if (uavgpos ne -1) then begin
          uavgstr = strmid(line,0,uavgpos-1)  
          reads, uavgstr, avg_u
          model.avg_u = avg_u
      endif

;     GET THE POWER PER HATOM
      powerpos = strpos(line,'= power/H (erg s-1)')
      if (powerpos ne -1) then begin
          powerstr = strmid(line,0,powerpos-1)  
          reads, powerstr, power_per_hatom
          model.power_per_hatom = power_per_hatom
      endif
  endwhile

  close,1

; IF ONLY THE MIN-MIN AND MIN-MAX MODELS ARE DESIRED (PER THE DL07
; RECOMMENDATIONS), CHECK IF THIS IS ONE OF THOSE
  if keyword_set(minmax_only) then begin
     if model.umax ne 1e6 and $
        model.umin ne model.umax then begin
        if keyword_set(quiet) eq 0 then $
           print, "Not a min-min or min-max model. Returning."
        fail = 1
        return, -1
     endif
  endif

; PRINT THE PROPERTIES OF THE MODEL IF WE'RE NOT IN HUSH-HUSH MODE
  if keyword_set(quiet) eq 0 then begin
      print, 'Model = ', model.grain_model
      print, 'Umin, Umax, Beta = ' $
        , model.umin, ', ', model.umax, ', ', model.beta
      print, 'Avg U = ', model.avg_u
      print, 'power per Hatom = ' $
        , model.power_per_hatom, ' ' $
        , 'mdust per mH = ', model.mdust_per_mh
  endif
    

; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
; READ THE MODEL CONVOLVED WITH BANDPASSES
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$

  readcol, filename $
    , skipline=11 $
    , band_microns $
    , power_erg_s_hatom $
    , intens_jy_ster_hatom $
    , telescope_name $
    , format='F,F,F,A' $
    , numline=29, /silent

; DIRBE/COBE
  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 1.270E+00, ct)
  if ct eq 0 then message, 'Missing DIRBE band 1' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_1.3')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 2.220E+00, ct)
  if ct eq 0 then message, 'Missing DIRBE band 2' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_2.2')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 3.530E+00, ct)
  if ct eq 0 then message, 'Missing DIRBE band 3' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_3.5')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 4.880E+00, ct)
  if ct eq 0 then message, 'Missing DIRBE band 4' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_4.9')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 1.229E+01, ct)
  if ct eq 0 then message, 'Missing DIRBE band 5' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_12')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 2.079E+01, ct)
  if ct eq 0 then message, 'Missing DIRBE band 6' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_21')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 5.599E+01, ct)
  if ct eq 0 then message, 'Missing DIRBE band 7' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_56')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 9.770E+01, ct)
  if ct eq 0 then message, 'Missing DIRBE band 8' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_98')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 1.479E+02, ct)
  if ct eq 0 then message, 'Missing DIRBE band 9' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_148')] = $
    intens_jy_ster_hatom[ind]

  ind = $
    where(telescope_name eq 'COBE-DIRBE' and band_microns eq 2.479E+02, ct)
  if ct eq 0 then message, 'Missing DIRBE band 10' $
  else $
    model.band_vals[where(model.bands eq 'dirbe_248')] = $
    intens_jy_ster_hatom[ind]

; IRAS
  ind = where(telescope_name eq 'IRAS' and band_microns eq 1.200E+01, ct)
  if ct eq 0 then message, 'Missing IRAS band 1' $
  else $
    model.band_vals[where(model.bands eq 'iras_12')] = $
    intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'IRAS' and band_microns eq 2.500E+01, ct)
  if ct eq 0 then message, 'Missing IRAS band 2' $
  else $
    model.band_vals[where(model.bands eq 'iras_25')] = $
    intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'IRAS' and band_microns eq 6.000E+01, ct)
  if ct eq 0 then message, 'Missing IRAS band 3' $
  else $
    model.band_vals[where(model.bands eq 'iras_60')] = $
    intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'IRAS' and band_microns eq 1.000E+02, ct)
  if ct eq 0 then message, 'Missing IRAS band 4' $
  else $
    model.band_vals[where(model.bands eq 'iras_100')] = $
    intens_jy_ster_hatom[ind]

; IRAC
  ind = where(telescope_name eq 'SST-IRAC' and band_microns eq 3.550E+00, ct)
  if ct eq 0 then message, 'Missing IRAC band 1' $
  else $
    model.band_vals[where(model.bands eq 'irac_3.6')] = $
    intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'SST-IRAC' and band_microns eq 4.493E+00, ct)
  if ct eq 0 then message, 'Missing IRAC band 2' $
  else $
    model.band_vals[where(model.bands eq 'irac_4.5')] = $
    intens_jy_ster_hatom[ind]
    
  ind = where(telescope_name eq 'SST-IRAC' and band_microns eq 5.731E+00, ct)
  if ct eq 0 then message, 'Missing IRAC band 3' $
  else $
    model.band_vals[where(model.bands eq 'irac_5.7')] = $
    intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'SST-IRAC' and band_microns eq 7.872E+00, ct)
  if ct eq 0 then message, 'Missing IRAC band 4' $
  else $
    model.band_vals[where(model.bands eq 'irac_7.9')] = $
    intens_jy_ster_hatom[ind]

; MIPS
  ind = where(telescope_name eq 'SST-MIPS' and band_microns eq 2.368E+01, ct)
  if ct eq 0 then $
     message, 'Missing MIPS band 1' $
  else $
    model.band_vals[where(model.bands eq 'mips_24')] = $
    intens_jy_ster_hatom[ind]    

  ind = where(telescope_name eq 'SST-MIPS' and band_microns eq 7.142E+01, ct)
  if ct eq 0 then $
     message, 'Missing MIPS band 2' $
  else $
    model.band_vals[where(model.bands eq 'mips_70')] = $
    intens_jy_ster_hatom[ind]    

  ind = where(telescope_name eq 'SST-MIPS' and band_microns eq 1.559E+02, ct)
  if ct eq 0 then $
     message, 'Missing MIPS band 3' $
  else $
    model.band_vals[where(model.bands eq 'mips_160')] = $
    intens_jy_ster_hatom[ind]

; HERSCHEL PACS
  ind = where(telescope_name eq 'Herschel-PACS' and band_microns eq 7.500E+01, ct)
  if ct eq 0 then begin
     message, 'Missing PACS band 1 : '+filename, /info
  endif else begin
     model.band_vals[where(model.bands eq 'pacs_70')] = $
        intens_jy_ster_hatom[ind]
  endelse

  ind = where(telescope_name eq 'Herschel-PACS' and band_microns eq 1.100E+02, ct)
  if ct eq 0 then $
     message, 'Missing PACS band 2 : '+filename, /info $
  else $
     model.band_vals[where(model.bands eq 'pacs_100')] = $
     intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'Herschel-PACS' and band_microns eq 1.700E+02, ct)
  if ct eq 0 then $
     message, 'Missing PACS band 3 : '+filename, /info $
  else $
     model.band_vals[where(model.bands eq 'pacs_160')] = $
     intens_jy_ster_hatom[ind]

; HERSCHEL SPIRE
  ind = where(telescope_name eq 'Herschel-SPIRE' and band_microns eq 2.500E+02, ct)
  if ct eq 0 then $
     message, 'Missing SPIRE band 1 : '+filename, /info $
  else $
     model.band_vals[where(model.bands eq 'spire_250')] = $
     intens_jy_ster_hatom[ind]
  
  ind = where(telescope_name eq 'Herschel-SPIRE' and band_microns eq 3.600E+02, ct)
  if ct eq 0 then $
     message, 'Missing SPIRE band 2 : '+filename, /info $
  else $
     model.band_vals[where(model.bands eq 'spire_350')] = $
     intens_jy_ster_hatom[ind]

  ind = where(telescope_name eq 'Herschel-SPIRE' and band_microns eq 5.200E+02, ct)
  if ct eq 0 then $
     message, 'Missing SPIRE band 3 : '+filename, /info $
  else $
     model.band_vals[where(model.bands eq 'spire_500')] = $
     intens_jy_ster_hatom[ind]

; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
; READ THE MODEL ITSELF
; %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$

  openr, 1, filename
  foundit = 0
  counter = 0
  while not(eof(1)) and (foundit eq 0) do begin      
      line=''
      readf,1,line
      if strmid(line,0,6) eq 'lambda' then foundit = 1 $
      else counter+=1
  endwhile
  close,1

  if foundit then begin
      readcol, filename $
        , skipline=counter $
        , lam_microns $
        , power_erg_s_hatom $
        , intens_jy_ster_hatom $
        , format='F,F,F' $
        , numline=1010, /silent
      
      model.lam = lam_microns
      model.int = intens_jy_ster_hatom
  endif else begin
      message, 'Did not find start of spectrum.'
  endelse

  fail = 0
  return, model

end
