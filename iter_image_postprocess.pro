;----------------------------------------------------------------------------------------
pro iter_image_postprocess, raw_image_file, processed_image_file, noise_threshold, zero_negatives = zero_negatives
;----------------------------------------------------------------------------------------
;
;--(dependencies)------------------------------------------------------------------------
; ***
;--(input)-------------------------------------------------------------------------------
; ***
;----------------------------------------------------------------------------------------
  ; read in raw file:
  raw_arr = readfits(raw_image_file, raw_hdr)

  apply_image_threshold, raw_arr, noise_threshold, zero_negatives = zero_negatives ; apply the threshold

  ; write out processed file:
  writefits, processed_image_file, raw_arr, raw_hdr

end
