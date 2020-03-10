;----------------------------------------------------------------------------------------
pro power_fraction_calc, masked_path_gas, masked_path_star $ ; input image (use masked image, but not regridded/smoothened)
                       , distance, inclination, astr_tolerance $
                       , array_dir $ ; directory to restore lambda array from
                       , fourier_length_conv $  ; fourier length conversion factor
                       , lambda $
                       , gas_flux_frac, star_flux_frac $
                       , gas_flux_frac_arr , star_flux_frac_arr $
                       , save_arrays = save_arrays
;----------------------------------------------------------------------------------------

;--(dependencies)------------------------------------------------------------------------
; *** To run, flux_fration_calc requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; ***
;--(output)------------------------------------------------------------------------------
; ***
;--(keywords)-----------------------------------------------------------------------------
; ***
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range


  fourier_power_fraction, masked_path_star $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , array_dir $ ; directory to restore lambda array from
    , fourier_length_conv $  ; fourier length conversion factor
    , lambda $
    , star_flux_frac $
    , star_flux_frac_arr

  fourier_power_fraction, masked_path_gas $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , array_dir $ ; directory to restore lambda array from
    , fourier_length_conv $  ; fourier length conversion factor
    , lambda $
    , gas_flux_frac $
    , gas_flux_frac_arr


  if save_arrays eq 1 then begin
    save, filename = array_dir + 'fgmcarr.sav', gas_flux_frac_arr ; fgmc is Heisenberg naming convention
    save, filename = array_dir + 'fclarr.sav', star_flux_frac_arr ; fcl is Heisenberg naming convention
  endif




end
