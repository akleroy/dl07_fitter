This directory contains a package for turning the Draine & Li (2007)
text file models into a grid of IDL structures. The area routines to
read the text files provided by B. Draine's website, to interpolate
between models and to define a grid of interpolated model values
spanning a range of radiation field formulations and PAH
fraction. There are also routines to compare measurements to a grid of
models to identify the best fit model.

The models include the IRAS, IRAC (Spitzer), MIPS (Spitzer), PACS
(Herschel), SPIRE (Herschel), and DIRBE (COBE) bandpasses.

EXAMPLE USE:

(1) Download
(2) Run "prep_models"
(3) Use "find_best_draine_model" to fit your data

Note that you may (probably do, in fact) want to modify prep_models to
turn off the /nospec option. Having this on only computes the model
grid for the various bands. Turning it off allows the grid to be
calculated for the full spectrum (and then convolved with an arbitrary
bandpass later).

PROGRAMS FOR USERS:

- prep models : read in the models for the Milky Way, LMC, and SMC;
  save them to disk; and then build them into model grids.

- mm09_dust : implement's Juan-Carlos Munoz-Mateos's fits to the
  Draine model grid. These give simple good enough estimates without
  requiring fitting using the Spitzer MIPS bands. See Munoz-Mateos+
  (2009).

SUPPORT PROGRAMS:

- grid_models: grid a set of models into a cube appropriate for
  fitting. The cube has axes of PAH fraction, PDR fraction (fPDR), and
  minimum radiation field.

- merge_draine_models : mostly a support program, but potentially of
  general use. Takes two Draine model structures and linearly averages
  them according to the provided weights. Used under the hood for
  building a model grid but could in general be used to mix any pair
  of models.

- read_one_model : reads on DL07 text file into a structure created by
  "new_draine_model".

- new_draine_model : makes an empty structure to contain a DL07 model.

TO DO:

- the fitter needs the most work

- a monte carlo tool to kick back realistic error estimates
