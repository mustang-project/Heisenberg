;----------------------------------------------------------------------------------------
function fourier_flux_loss_correction, kernel, cut_ratio
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
    qcon = 1.0 ; not yet implemented
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



; function gwfhm_multiple_flux_frac, gwfhm_multiple, path_prefix   ; , path_prefix
;   readcol, path_prefix + '/Idldir/F_paper_routines_new/F_paper_git/gauss_fraction_data.dat', skipline =1, mult, frac
;
;   msort = mult[ sort(mult)]
;   fsort = frac[ sort(mult)]
;
;
;   return, interpol(fsort, msort, gwfhm_multiple) ; return fraction of signal that remians
;
; end
