;----------------------------------------------------------------------------------------
pro fourier_diffuse_fraction, kernel, butterworth_order, pass $ ; description of the filter   $
                            , image_file_path $ ; input image
                            , distance, inclination, astr_tolerance $ ;
                            , fourier_length_conv $  ; fourier length conversion factor
                            , lambda $ ; lambda
                            , flux_frac, diffuse_frac $ ; output flux and diffuse fraction
                            , lambda_errmin = lambda_errmin, lambda_errmax = lambda_errmax $ ; lambda errors
                            , diffuse_frac_errmax = diffuse_frac_errmax, diffuse_frac_errmin = diffuse_frac_errmin $ ; output diffuse fraction errors
                            , flux_frac_errmax = flux_frac_errmax, flux_frac_errmin = flux_frac_errmin ; output flux fraction errors
;----------------------------------------------------------------------------------------
; calculates the diffuse fraciton in a .fits image
; NOTE errors and multiple cut-lengths not compatible -> innefficient multi-fft for multi-cuts
;--(dependencies)------------------------------------------------------------------------
; *** To run, fourier_diffuse_fraction requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; *** kernel              = the kernel with which to filter to filter
; ***                       * 'gaussian'
; ***                       * 'butterworth'
; ***                       * 'ideal'
; *** butterworth_order   = The order of the butterworth filter (only relevant if
; ***                       kernel = butterworth, but a dummy value must be provided
; ***                       regardless)
; *** pass                = 'low' or 'high', i.e. lowpass filtering or highpass filtering
; *** image_file_path     = file path to the image to be filtered. The image should be a
; ***                       linear multiple of the flux
; *** distance            = distance to the galaxy
; *** inclination         = inclination of the galaxy [radians]
; *** astr_tolerance      = allowable tolerance in decimal degrees between astrometric
; ***                       position of image pixels (e.g. 1.8" = 5e-4 decimal degrees)
; *** fourier_length_conv = fourier length conversion factor
; *** lambda              = best fitting lambda (region seperation length)
;--(output)------------------------------------------------------------------------------
; *** flux_frac    = non-diffuse fraction in the input image given the measured
; *** diffuse_frac = diffuse fraction in the input image given the measured
;--(keywords)-----------------------------------------------------------------------------
; ***  lambda_errmin      = upward error on lambda         [input] (neccessary for error calculations)
; ***  lambda_errmax      = downward error on lambda       [input] (neccessary for error calculations)
; *** diffuse_frac_errmax = upward error on diffuse_frac   [output]
; *** diffuse_frac_errmin = downward error on diffuse_frac [output]
; *** flux_frac_errmax    = upward error on flux_frac      [output]
; *** flux_frac_errmin    = downward error on flux_frac    [output]
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  ; if n_params() eq 10 then calc_errors = 1 else $
  ;   if n_params() eq 8 then calc_errors = 0 else begin
  ;     print, "fourier_diffuse_fraction: incorrect number of parameters specified "
  ;     stop
  ;   endelse

  if n_elements(lambda_errmin) eq 1 && n_elements(lambda_errmax) eq 1 then calc_errors = 1 else calc_errors = 0

  if inclination gt (2*!pi) then print, "fourier_diffuse_fraction: warning: inclination must be in radians. supplied inclination is greater than 2*Pi"


  ; ************************************************
  ; input image and get important parameters
  ; ************************************************
  image_arr = readfits(image_file_path, image_hdr)

  if size(image_arr, /type) eq 4 then begin
    zero_flux = 0.0e0
    one_flux = 1.0d0
    image_double = 0
  endif else if size(image_arr, /type) eq 5 then begin
    zero_flux = 0.0d0
    one_flux = 1.0d0
    image_double = 1
  endif


  getrot, image_hdr, rotation, cdelt_var ; get cdelt value ; only single precision
  if n_elements(cdelt_var) eq 2 && abs(abs(cdelt_var[0]) - abs(cdelt_var[1])) le astr_tolerance then begin
    cdelt = mean(abs(cdelt_var))
  endif else begin
    print, "image not square"
    stop
  endelse


  pix_to_pc = distance*tan(!dtor*cdelt)/sqrt(cos(inclination))

  pixel_cut_length_mid = (lambda/pix_to_pc) * fourier_length_conv ; convert from pc to pixels and then apply length conversion factor

  if calc_errors eq 1 then begin
    pixel_cut_length_min = ((lambda - lambda_errmin)/pix_to_pc) * fourier_length_conv ; convert from pc to pixels and then apply length conversion factor
    pixel_cut_length_max = ((lambda + lambda_errmax)/pix_to_pc) * fourier_length_conv ; convert from pc to pixels and then apply length conversion factor
  endif


  image_total = total(image_arr, /nan)

  ; ************************************************
  ; cut_length = lambda
  ; ************************************************

  if n_elements(pixel_cut_length_mid) eq 1 then begin

    fourier_filter_tool, kernel, butterworth_order, pass, pixel_cut_length_mid $ ; description of the filter
     , image_arr, image_hdr $ ; input image
     , filtered_image_arr = filtered_image_arr $ ; filtered image output variablse
     , filtered_image_path = filtered_image_path ; filtered image output paths

    lambda_mid_total = total(filtered_image_arr[where(filtered_image_arr ge zero_flux)], /nan) ; total of only non_negative flux
    lambda_mid_frac = lambda_mid_total/image_total ; % remaining flux
  endif else if n_elements(pixel_cut_length_mid) gt 1 then begin

    if image_double eq 0 then begin
      lambda_mid_frac = fltarr(n_elements(pixel_cut_length_mid))
    endif else if image_double eq 1 then begin
      lambda_mid_frac = dblarr(n_elements(pixel_cut_length_mid))
    endif

    for ii = 0, n_elements(pixel_cut_length_mid)-1 do begin ; for loop leads to re-fourier filtering of the image each time, more efficient version is possilbe
      pclm_i = pixel_cut_length_mid[ii]
      fourier_filter_tool, kernel, butterworth_order, pass, pclm_i $ ; description of the filter
        , image_arr, image_hdr $ ; input image
        , filtered_image_arr = filtered_image_arr $ ; filtered image output variablse
        , filtered_image_path = filtered_image_path ; filtered image output paths

      lambda_mid_total = total(filtered_image_arr[where(filtered_image_arr ge zero_flux)], /nan) ; total of only non_negative flux
      lambda_mid_frac[ii] = lambda_mid_total/image_total ; % remaining flux
    endfor

  endif else begin
    print, 'error'
    stop
  endelse


  ; ************************************************
  ; cut_length = lambda + lambda_errmax
  ; ************************************************
  if calc_errors eq 1 then begin
    fourier_filter_tool, kernel, butterworth_order, pass, pixel_cut_length_max $ ; description of the filter
                      , image_arr, image_hdr $ ; input image
                      , filtered_image_arr = filtered_image_arr, filtered_image_hdr = filtered_image_hdr $ ; filtered image output variablse
                      , filtered_image_path = filtered_image_path ; filtered image output paths

    lambda_max_total = total(filtered_image_arr[where(filtered_image_arr ge zero_flux)], /nan) ; total of only non_negative flux (have to remove negative flux else lambda_max_total~0 due to removal of dc component )
    lambda_max_frac = lambda_max_total/image_total
  endif


  ; ************************************************
  ; cut_length = lambda - lambda_errmin
  ; ************************************************
  if calc_errors eq 1 then begin
    fourier_filter_tool, kernel, butterworth_order, pass, pixel_cut_length_min $ ; description of the filter
                      , image_arr, image_hdr $ ; input image
                      , filtered_image_arr = filtered_image_arr, filtered_image_hdr = filtered_image_hdr $ ; filtered image output variablse
                      , filtered_image_path = filtered_image_path ; filtered image output paths

    lambda_min_total = total(filtered_image_arr[where(filtered_image_arr ge zero_flux)], /nan) ; total of only non_negative flux
    lambda_min_frac = lambda_min_total/image_total   ; minimum signal flux left -> max diffuse flux removed
  endif

  ; ************************************************
  ; calculate flux fraction + errors
  ; ************************************************
  flux_frac = lambda_mid_frac
  if calc_errors eq 1 then begin
    flux_frac_errmax = lambda_max_frac - lambda_mid_frac
    flux_frac_errmin = lambda_mid_frac - lambda_min_frac
  endif else begin
    flux_frac_errmax = !values.f_nan ; consider rewrite to avoid this
    flux_frac_errmin = !values.f_nan ; consider rewrite to avoid this
  endelse

  ; ************************************************
  ; calculate diffuse fraction + errors
  ; ************************************************
  diffuse_frac_mid = one_flux - lambda_mid_frac ; 1-lambda_mid_frac ; one_flux takes account of double/single precision
  diffuse_frac = diffuse_frac_mid


  if calc_errors eq 1 then begin
    diffuse_frac_max = one_flux - lambda_min_frac  ; minimum signal flux left -> max diffuse flux removed
    diffuse_frac_min = one_flux - lambda_max_frac  ; maximum signal flux left -> min diffuse flux removed

    diffuse_frac_errmax = diffuse_frac_max - diffuse_frac_mid
    diffuse_frac_errmin = diffuse_frac_mid - diffuse_frac_min
  endif else begin
    diffuse_frac_errmax = !values.f_nan ; consider rewrite to avoid this
    diffuse_frac_errmin = !values.f_nan ; consider rewrite to avoid this
  endelse


end
