;----------------------------------------------------------------------------------------
function fourier_flux_loss_correction, kernel, cut_ratio
;----------------------------------------------------------------------------------------
; calculates the corrective factor to apply to the measured compact fraction due
; to flux loss as a result of the application of a filter to a single isolated gaussian
; region. The function implements the fit presented as equtation 31 in Hygate (2019)
;--(dependencies)------------------------------------------------------------------------
; *** unity_symmetric_sigmoidal.pro
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'
;--(output)-----------------------------------------------------------------------------
; *** qcon             = the compact fraction corrective factor to account for flux loss
; ***                    from single peaks
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if kernel eq 'butterworth' || kernel eq 'b' then begin
    qcon = 1.0 ; not yet implemented

    ; 1st order B: -0.0381051      1.30460      4.44501
    ; 2nd order B: 0.00191006      2.37707      3.43289

  endif else if kernel eq 'gaussian' || kernel eq 'gauss' || kernel eq 'g' then begin
    a = -0.0159072
    b = 1.68958
    c = 4.85505
    qcon = unity_symmetric_sigmoidal(cut_ratio,a,b,c)
  endif else if kernel eq 'ideal' || kernel eq 'i' then begin
    qcon = 1.0 ; not yet implemented
  endif

  return, qcon


end
