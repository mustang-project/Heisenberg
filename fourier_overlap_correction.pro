;----------------------------------------------------------------------------------------
function fourier_overlap_correction, kernel, peak_id_file_path, dist_stat_pass = dist_stat_pass, $
  qover_sigma = qover_sigma, dist_stat_sigma_pass = dist_stat_sigma_pass
;----------------------------------------------------------------------------------------
; calculate the corrective factor to apply to the measured compact fraction due
; to flux loss as a result of overlapping peaks
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'
; ***                     * 'ideal'
; *** peak_id_file_path = clumpfind peak output file for the image to be analysed
;--(keywords)----------------------------------------------------------------------------
; *** dist_stat_pass    = allows the value of dist_stat to be passed to routine
; ***                     (if not set dist_stat is calculated)
; *** qover_sigma       = the 1-sigma error on qover
;--(output)-----------------------------------------------------------------------------
; *** qover             = the compact fraction corrective factor to account for
; ***                     overlapping peaks
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if kernel eq 'butterworth' || kernel eq 'b' then begin
    qover = 1.0 ; not yet implemented
  endif else if kernel eq 'gaussian' || kernel eq 'gauss' || kernel eq 'g' then begin
    a = -0.35443467
    b = 3.1018254
    c = 1.40684605
    d = 0.98407271


    if n_elements(dist_stat_pass) eq 1 then dist_stat = dist_stat_pass else dist_stat = med_peak_relative_nearest_neighbour_dist_calc(peak_id_file_path, dist_med_sigma = dist_stat_sigma)

    qover = symmetric_sigmoidal(dist_stat,a,b,c,d)

    if qover gt 1.0d0 then qover = 1.0d0

    ; covariance matrix for the fitting parameters a,b,c and d:
    covariance_matrix = [[ 2.15940168e-01, 4.86565784e-02, 1.02120890e-01, -1.36567146e-04], $
     [ 4.86565784e-02, 1.30124103e-02, 2.41226648e-02, -3.82292032e-05], $
     [ 1.02120890e-01, 2.41226648e-02, 4.89232032e-02, -6.81607721e-05], $
     [-1.36567146e-04, -3.82292032e-05, -6.81607721e-05, 1.22697461e-07]]

    qover_sigma = fourier_overlap_sigma_calc(dist_stat, a, b, c, d, covariance_matrix, dist_stat_sigma)


  endif else if kernel eq 'ideal' || kernel eq 'i' then begin
    qover = 1.0 ; not yet implemented
  endif

  return, qover
end
