;----------------------------------------------------------------------------------------
function med_peak_relative_nearest_neighbour_dist_calc, cords_filepath, $
  dist_val = dist_val, fwhm_val = fwhm_val
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

  readcol, cords_filepath, ncl, xx, yy, flux, xfwhm_vec, yfwhm_vec, /silent
  fwhm_vec = (xfwhm_vec +  yfwhm_vec)/2
  fwhm_val = mean([xfwhm_vec, yfwhm_vec]) ; mean fwhm of the peaks in the image

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

  ; dist_med = median(min(dist_arr, dimension = 2, /nan)) ; get the median "nearest neighbour" distance for a peak
  dist_val = median(min(dist_arr, dimension = 2, /nan))
  dist_med = dist_val/fwhm_val ; normalise to the fwhm_val

  return, dist_med
  
end
