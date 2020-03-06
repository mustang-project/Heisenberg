;----------------------------------------------------------------------------------------
function fourier_zeta_correction, kernel, eta
;----------------------------------------------------------------------------------------
; calculate the corrective factor to apply to the measured compact fraction due to flux
; loss as a result of overlapping peaks
; NOTE: this is the method presented in Hygate (2019) and is depreciated in favour of
;       fourier_overlap_correction.pro
;--(dependencies)------------------------------------------------------------------------
; *** none
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'
; ***                     * 'ideal'
; *** eta               = the evolutionary phase lifetime adjusted filling factor
; ***                     equations 33 and 34 (Hygate 2019)
;--(output)-----------------------------------------------------------------------------
; *** qeta              = the compact fraction corrective factor to account for
; ***                     overlapping peaks as caculated based on eta
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if kernel eq 'butterworth' || kernel eq 'b' then begin
    qzeta = 1.0 ; not yet implemented
  endif else if kernel eq 'gaussian' || kernel eq 'gauss' || kernel eq 'g' then begin
    intercept = 1.0439514d0
    slope = -1.1516489d0
    qzeta = (slope * eta) + intercept
    if qzeta gt 1.0d0 then qzeta = 1.0d0
  endif else if kernel eq 'ideal' || kernel eq 'i' then begin
    qzeta = 1.0 ; not yet implemented
  endif

  return, qzeta
end
