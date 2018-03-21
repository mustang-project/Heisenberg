;----------------------------------------------------------------------------------------
function fourier_zeta_correction, kernel, eta
;----------------------------------------------------------------------------------------
; apply correction for
;--(dependencies)------------------------------------------------------------------------
; *** To run, xxx requires:
; ***
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'

;--(keywords)-----------------------------------------------------------------------------
; *** filtered_image_path = filepath to output the filtered image to
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if kernel eq 'butterworth' || kernel eq 'b' then begin
    qzeta = 1.0 ; not yet implemented
  endif else if kernel eq 'gaussian' || kernel eq 'gauss' || kernel eq 'g' then begin
    slope = -1.3157485
    qzeta = lin_unity_intercept(eta,slope)
  endif else if kernel eq 'ideal' || kernel eq 'i' then begin
    qzeta = 1.0 ; not yet implemented
  endif



  ; zeta_slope =
  ; zeta_intercept =
  ; qzeta = (zeta * zeta_slope) + zeta_intercept

  return, qzeta
end
