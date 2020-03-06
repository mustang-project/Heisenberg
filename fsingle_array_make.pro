;----------------------------------------------------------------------------------------
pro fsingle_array_make, array_dir, fcl, fgmc
;----------------------------------------------------------------------------------------
; when fcl and fgmc are single-valued (not calculated using fourier filtering) this
; routine writes out arrays of size(n_elements(lambdarr)) for use in derivephys.pro
;--(input)-------------------------------------------------------------------------------
; *** array_dir = directory to restore lambda array from
; *** fcl       = non-diffuse fraction in the gas image
; *** fgmc      = non-diffuse fraction in the star image
;--(output files)------------------------------------------------------------------------
; *** fclarr.sav  = array fcl of size(n_elements(lambdarr)) for use in derivephys.pro
; *** fgmcarr.sav = array fgmc of size(n_elements(lambdarr)) for use in derivephys.pro
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  restore, filename = array_dir + 'lambdaarr.sav' ; restores array [lambdarr], which is a grid of lambda values of size ntry spanning the parameter search space.

  star_flux_frac_arr = fltarr(n_elements(lambdaarr)) + fcl
  gas_flux_frac_arr = fltarr(n_elements(lambdaarr)) + fgmc

  save, filename = array_dir + 'fclarr.sav', star_flux_frac_arr ; fcl is kl14 naming convention
  save, filename = array_dir + 'fgmcarr.sav', gas_flux_frac_arr ; fgmc is kl14 naming convention

end
