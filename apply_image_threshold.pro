;----------------------------------------------------------------------------------------
pro apply_image_threshold, raw_arr, noise_threshold, zero_negatives = zero_negatives
;----------------------------------------------------------------------------------------
;
;--(input)-------------------------------------------------------------------------------
; ***   =
;--(keywords)----------------------------------------------------------------------------
; ***     =
;--(output)------------------------------------------------------------------------------
; ***            =
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  if n_elements(noise_threshold) eq 1 then begin ; zero values below the noise threshold
    if noise_threshold lt 0.0e0 then f_error, 'iter_image_postprocess: a negative noise_threshold has been supplied for: ' + raw_image_file + '. a positive (or zero) threshold must be suplied'
    zero_list = where(raw_arr lt noise_threshold, zero_count)
    if zero_count ne 0 then raw_arr[zero_list] = 0.0
  endif

  if n_elements(zero_negatives) eq 1 && zero_negatives eq 1 then begin ; zero negative values
    zero_list = where(raw_arr lt 0.0e0, zero_count)
    if zero_count ne 0 then raw_arr[zero_list] = 0.0
  endif



end
