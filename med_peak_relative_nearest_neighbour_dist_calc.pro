;----------------------------------------------------------------------------------------
function med_peak_relative_nearest_neighbour_dist_calc, cords_filepath, $
  dist_val = dist_val, fwhm_val = fwhm_val, dist_med_sigma = dist_med_sigma
;----------------------------------------------------------------------------------------
; cacluate median peak nearest neighbour distance for an image
;--(input)-------------------------------------------------------------------------------
; *** cords_filepath  =  clumpfind peak output file for the image to be analysed
;--(keywords)----------------------------------------------------------------------------
; *** dist_val        = the median nearest neighbour distance not divided by the
; ***                   fwhm
; *** fwhm_val        = the mean fwhm as calculated from clumpfind
;--(output)------------------------------------------------------------------------------
; ***  dist_med       = the median nearest neighbour distance for peaks in the
; ***                   image being analysed
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  if n_params() lt 1 then begin ; check for parameters
   print, 'Syntax: med_peak_relative_nearest_neighbour_dist_calc, cords_filepath,'
   print, 'dist_val = , fwhm_val = , dist_med_sigma = '
   return, 'error'
  endif


  ; *** read the co-ordinates filepath
  readcol, cords_filepath, ncl, xx, yy, flux, xfwhm_vec, yfwhm_vec, /silent


  ; *** calculate the FWHM
  fwhm_vec = (xfwhm_vec +  yfwhm_vec)/2
  fwhm_val = mean([xfwhm_vec, yfwhm_vec]) ; mean fwhm of the peaks in the image
  fwhm_std = robust_sigma([xfwhm_vec, yfwhm_vec]) ; standard deviation of the fwhm of the peaks

  ; *** bootstrap fwhm to get the error
  f_vec = [xfwhm_vec, yfwhm_vec]
  nfs = n_elements(f_vec)
  n_bootstrap = 10000
  fboot_quant = fltarr(n_bootstrap)
  for bb =  0, n_bootstrap-1, 1 do begin
    randomindices = round(randomu(bootseed,nfs) * (nfs-1.0)) ; update this
    fboot_quant[bb] = mean(f_vec[randomindices])
  endfor


  bfwhm_val = median(fboot_quant)
  bfwhm_lq = cgpercentiles(fboot_quant, percentiles =[0.16])
  bfwhm_uq = cgpercentiles(fboot_quant, percentiles =[0.84])
  fwhm_sigma = mean(abs([bfwhm_lq - bfwhm_val, bfwhm_uq - bfwhm_val ]) )



  ; *** calculate distances between peaks
  npoints = n_elements(xx)
  dist_arr = fltarr(npoints,npoints)

  for ii = 0, npoints -1, 1 do begin ; loop over points
    xval = xx[ii]
    yval = yy[ii]

    for jj = 0, npoints -1, 1 do begin ; loop over all other points
      if ii eq jj then begin
        dist_arr[ii,jj] = !values.f_nan
        continue ; skip the point itself
      endif
      xcom = xx[jj]
      ycom = yy[jj]
      dist = sqrt((xval - xcom)^2 + (yval - ycom)^2)
      dist_arr[ii,jj] = dist

    endfor
  endfor

  ; *** get the median "nearest neighbour" distance for a peak
  min_dist_vec = min(dist_arr, dimension = 2, /nan)
  dist_val = median(min_dist_vec)

  ; *** bootstrap min_dist_vec to get the error on dist_val
  ndists = n_elements(min_dist_vec)
  n_bootstrap = 10000
  boot_quant = fltarr(n_bootstrap)
  for bb =  0, n_bootstrap-1, 1 do begin
    randomindices = round(randomu(bootseed,ndists) * (ndists-1)) ; update this
    boot_quant[bb] = median(min_dist_vec[randomindices])
  endfor

  ; dist_val = median(boot_quant)
  dist_lq = cgpercentiles(boot_quant, percentiles =[0.16])
  dist_uq = cgpercentiles(boot_quant, percentiles =[0.84])
  dist_sigma = mean(abs([dist_lq -dist_val, dist_uq -dist_val ]) )


  ; *** calculate final value and error
  dist_med = dist_val/fwhm_val ; normalise to the fwhm_val
  dist_med_sigma = dist_med * sqrt(  (dist_sigma/dist_val)^2 + (fwhm_sigma/fwhm_val)^2   )

  return, dist_med

end
