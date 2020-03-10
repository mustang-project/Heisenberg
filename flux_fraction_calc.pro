;----------------------------------------------------------------------------------------
pro flux_fraction_calc, kernel, butterworth_order, pass  $ ; description of the filter   $
                      , masked_path_gas, masked_path_star $ ; input image (use masked image, but not regridded/smoothened)
                      , distance, inclination, astr_tolerance $
                      , array_dir $ ; directory to restore lambda array from
                      , fourier_length_conv $  ; fourier length conversion factor
                      , lambda $
                      , gas_flux_frac, star_flux_frac $
                      , gas_flux_frac_arr , star_flux_frac_arr $
                      , save_arrays = save_arrays $
                      , gas_noise_threshold = gas_noise_threshold $
                      , star_noise_threshold = star_noise_threshold $
                      , mask_file = mask_file
;----------------------------------------------------------------------------------------
; calculates the non-diffuse-fraction in a .fits image over the range of values of lambda
; NOTE errors and multiple cut-lengths not compatible -> innefficient multi-fft for multi-cuts
;--(dependencies)------------------------------------------------------------------------
; *** To run, flux_fration_calc requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; *** kernel              = the kernel with which to filter to filter
; ***                       * 'gaussian'
; ***                       * 'butterworth'
; ***                       * 'ideal'
; *** butterworth_order   = The order of the butterworth filter (only relevant if
; ***                       kernel = butterworth, but a dummy value must be provided
; ***                       regardless)
; *** pass                = 'low' or 'high', i.e. lowpass filtering or highpass filtering
; *** masked_path_gas     = file path to the gas image to be filtered. The image should be
; ***                       a linear multiple of the flux
; *** masked_path_star    = file path to the star image to be filtered. The image should
; ***                       be a linear multiple of the flux
; *** distance            = distance to the galaxy
; *** inclination         = inclination of the galaxy [radians]
; *** astr_tolerance      = allowable tolerance in decimal degrees between astrometric
; ***                       position of image pixels (e.g. 1.8" = 5e-4 decimal degrees)
; *** array_dir           = directory to restore lambda array from
; *** fourier_length_conv = fourier length conversion factor
; *** lambda              = best fitting lambda (region seperation length)
;--(output)------------------------------------------------------------------------------
; *** gas_flux_frac      = non-diffuse fraction in the gas image given the measured
; *** star_flux_frac     = non-diffuse fraction in the star image given the measured
; *** gas_flux_frac_arr  = array of gas flux fractions at values of lambda in lambdaarr
; *** star_flux_frac_arr = array of star flux fractions at values of lambda in lambdaarr
;--(keywords)-----------------------------------------------------------------------------
; ***  save_arrays = (1/0) [1] write gas_flux_frac_arr and star_flux_frac_arr to a file.
; ***                      [0] do nothing
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  ; ************************************************
  ; get flux fractions at best fitting lambda
  ; ************************************************
  cut_length = lambda ; in pc (fourier_diffuse_fraction handles conversion to pixels)

  fourier_diffuse_fraction, kernel, butterworth_order, pass  $ ; description of the filter   $
    , masked_path_gas $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , fourier_length_conv $  ; fourier length conversion factor
    , cut_length $ ; lambda (minimal call without error calculation)
    , gas_flux_frac $
    , noise_threshold = gas_noise_threshold $
    , mask_file = mask_file


  fourier_diffuse_fraction, kernel, butterworth_order, pass $ ; description of the filter   $
    , masked_path_star $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , fourier_length_conv $  ; fourier length conversion factor
    , cut_length $ ; lambda (minimal call without error calculation)
    , star_flux_frac $
    , noise_threshold = star_noise_threshold $
    , mask_file = mask_file

  ; ************************************************
  ; evaluate flux fractions over lambdaarr values
  ; ************************************************
  restore, filename= array_dir + 'lambdaarr.sav' ; restores array [lambdarr], which is a grid of lambda values of size ntry spanning the parameter search space.

  cut_lengths = lambdaarr ; multiple values in pc (fourier_diffuse_fraction handles conversion to pixels)

  fourier_diffuse_fraction, kernel, butterworth_order, pass $ ; description of the filter   $
    , masked_path_gas $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , fourier_length_conv $  ; fourier length conversion factor
    , cut_lengths $ ; lambda (minimal call without error calculation)
    , gas_flux_frac_arr $
    , noise_threshold = gas_noise_threshold $
    , mask_file = mask_file


  fourier_diffuse_fraction, kernel, butterworth_order, pass $ ; description of the filter   $
    , masked_path_star $ ; input image (use masked image, but not regridded/smoothened)
    , distance, inclination, astr_tolerance $
    , fourier_length_conv $  ; fourier length conversion factor
    , cut_lengths $ ; lambda (minimal call without error calculation)
    , star_flux_frac_arr $
    , noise_threshold = star_noise_threshold $
    , mask_file = mask_file

    if save_arrays eq 1 then begin
      save, filename = array_dir + 'fgmcarr.sav', gas_flux_frac_arr ; fgmc is Heisenberg naming convention
      save, filename = array_dir + 'fclarr.sav', star_flux_frac_arr ; fcl is Heisenberg naming convention
    endif

end
