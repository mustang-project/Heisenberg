;----------------------------------------------------------------------------------------
function fourier_overlap_correction, kernel, peak_id_file_path, dist_stat_pass = dist_stat_pass
;----------------------------------------------------------------------------------------
; calculate the corrective factor to apply to the measured compact fraction due
; to flux loss as a result of overlapping peaks
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'
; ***                     * 'ideal'
; *** peak_id_file_path = clumpfind peak output file for the image to be analysed
;--(output)-----------------------------------------------------------------------------
; *** qover             = the compact fraction corrective factor to account for
; ***                     overlapping peaks
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if kernel eq 'butterworth' || kernel eq 'b' then begin
    qover = 1.0 ; not yet implemented
  endif else if kernel eq 'gaussian' || kernel eq 'gauss' || kernel eq 'g' then begin
    a = -0.090760350
    b = 3.0986142
    c = 1.6054474

    if n_elements(dist_stat_pass) eq 1 then dist_stat = dist_stat_pass else dist_stat = med_peak_relative_nearest_neighbour_dist_calc(peak_id_file_path)
    qover = unity_symmetric_sigmoidal(dist_stat,a,b,c)
    if qover gt 1.0d0 then qover = 1.0d0
  endif else if kernel eq 'ideal' || kernel eq 'i' then begin
    qover = 1.0 ; not yet implemented
  endif

  return, qover
end
