;----------------------------------------------------------------------------------------
pro fourier_power_fraction, image_path $ ; input image (use masked image, but not regridded/smoothened)
                          , distance, inclination, astr_tolerance $
                          , array_dir $ ; directory to restore lambda array from
                          , fourier_length_conv $  ; fourier length conversion factor
                          , lambda $
                          , fit_power_frac $
                          , power_frac_arr
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


  ; ************************************************
  ; input image and get important parameters
  ; ************************************************
  image_arr = readfits(image_path, image_hdr)


  getrot, image_hdr, rotation, cdelt_var ; get cdelt value ; only single precision
  if n_elements(cdelt_var) eq 2 && abs(abs(cdelt_var[0]) - abs(cdelt_var[1])) le astr_tolerance then begin
    cdelt = mean(cdelt_var)
  endif else begin
    print, "image not square"
    stop
  endelse

  pix_to_pc = distance*tan(!dtor*cdelt)/sqrt(cos(inclination))


  ; ************************************************
  ; cut_length = lambda
  ; ************************************************
  fourier_power_fraction_vec, (lambda/pix_to_pc) * fourier_length_conv $ ; vector of cut_lengths
                            , image_arr, image_hdr $ ; input image
                            , partial_power_frac = partial_power_frac

  fit_power_frac = partial_power_frac[0]

  ; ************************************************
  ; evaluate power fractions over lambdaarr values
  ; ************************************************

  restore, filename=array_dir + 'lambdaarr.sav' ; restores array [lambdarr], which is a grid of lambda values of size ntry spanning the parameter search space.

  fourier_power_fraction_vec, (lambdaarr/pix_to_pc) * fourier_length_conv $ ; vector of cut_lengths
                            , image_arr, image_hdr $ ; input image
                            , partial_power_frac = power_frac_arr





end
