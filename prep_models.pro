pro prep_models

; READ THE MODELS INTO IDL SAVE FILES
  read_all_models, gal="MW"
  read_all_models, gal="LMC"
  read_all_models, gal="SMC"

; CREATE MODEL GRIDS
  grid_models, gal="MW", /nospec
  grid_models, gal="LMC", /nospec
  grid_models, gal="SMC", /nospec

end
