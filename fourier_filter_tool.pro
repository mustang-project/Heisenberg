;----------------------------------------------------------------------------------------
pro fourier_filter_tool, kernel, butterworth_order, pass, cut_length $ ; description of the filter
                       , image_input, image_hdr $ ; input image
                       , filtered_image_arr = filtered_image_arr, filtered_image_hdr = filtered_image_hdr $ ; filtered image output variablse
                       , filtered_image_path = filtered_image_path $ ; filtered image output paths
                       , negative_image_arr = negative_image_arr, negative_image_hdr = negative_image_hdr $ ; filtered image output variablse
                       , negative_image_path = negative_image_path ; filtered image output paths
;----------------------------------------------------------------------------------------
; Tool for Fourier filtering a .fits image in Fourier space by applying a lowpass or
; highpass filter with a number of different kernel definitions. This has the effect of
; attenuating the small-scale structure (lowpass filter) or attenuating the large-scale
; structure (highpass filter)
;--(dependencies)------------------------------------------------------------------------
; *** To run, fourier_filter_tool requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; *** kernel            = the kernel with which to filter to filter
; ***                     * 'gaussian'
; ***                     * 'butterworth'
; ***                     * 'ideal'
; *** butterworth_order = The order of the butterworth filter (only relevant if
; ***                     kernel = butterworth, but a dummy value must be provided
; ***                     regardless)
; *** pass              = 'low' or 'high', i.e. lowpass filtering or highpass filtering
; *** cut_length        = [pixels] the "cosine" wavelength used to determine the
; ***                     crititicalfrequency (freq = 1/cut_length) of the filter
; *** image_input       = the file path of the image to be fourier filtered e.g
; ***                     './images/image.fits'
; *** image_hdr         = (optional) The header of the image. If this paramter is suppied
; ***                     then image_input is the image in an array and not a filepath
;--(keywords)-----------------------------------------------------------------------------
; *** filtered_image_path = filepath to output the filtered image to
; *** filtered_image_arr  = variable to output the filtered image array to
; *** filtered_image_hdr  = variable to output the filtered image .fits header to
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  ; ************************************************
  ; test input validity
  ; ************************************************
  if n_params() lt 5 then begin
     print,'Syntax: fourier_filter_tool, kernel, butterworth_order, pass, cut_length, image_input, [image_hdr'
     print,'                             filtered_image_arr = , filtered_image_hdr =      '
     print,'                             filtered_image_path = , filtered_image_path =    '
     return
  endif else if (n_params() eq 5) then image_arr = readfits(image_input, image_hdr) else begin ; read in image
    image_arr = image_input ; image and hdr supplied to command
  endelse

  if size(image_arr, /n_dimensions) ne 2 then begin
    print, ' error: Fourier_filter.pro only works on a two dimensional image'
    print, ' quitting...'
    stop
  endif

  ; ************************************************
  ; deal with precision of image
  ; ************************************************
  if size(image_arr, /type) eq 4 then begin
    zero_flux = 0.0e0
    image_double = 0
  endif else if size(image_arr, /type) eq 5 then begin
    zero_flux = 0.0d0
    image_double = 1
  endif

  ; ************************************************
  ; take fourier transform and deal with nans
  ; ************************************************
  nan_list = where(finite(image_arr, /nan) eq 1, nan_count)
  if (nan_count gt 1.0) then image_arr[nan_list] = zero_flux ; set nan_pixes to zero

  fft_arr = fft(image_arr)

  ; ************************************************
  ; Create frequency vectors
  ; ************************************************

  image_dimensions = size(image_arr, /dimensions)
  image_xsize      = image_dimensions[0]
  image_ysize      = image_dimensions[1]


  ; deal with integer cut lengths to prevent errors, convert to double to be safe
  ; and protect original input variable
  ; this assumes the desired cut lengths in pixels is an integer with fewer than
  ; ~15 significant digits (limit of double precision significance)
  if is_an_integer(cut_length) eq 1 then cut_length_temp = double(cut_length) else cut_length_temp = cut_length


  if size(cut_length_temp, /type) eq 4 then begin
    image_xvect = findgen((image_xsize - 1)/2) + 1
    image_yvect = findgen((image_ysize - 1)/2) + 1

    if ~keyword_set(image_xsamp_interval) then image_xsamp_interval = 1.0e0
    if ~keyword_set(image_ysamp_interval) then image_ysamp_interval = 1.0e0

    zero_freq = 0.0e0
    one_freq = 1.0e0
    freq_double = 0

  endif else if size(cut_length_temp, /type) eq 5 then begin
    image_xvect = dindgen((image_xsize - 1)/2) + 1
    image_yvect = dindgen((image_ysize - 1)/2) + 1

    if ~keyword_set(image_xsamp_interval) then image_xsamp_interval = 1.0d0
    if ~keyword_set(image_ysamp_interval) then image_ysamp_interval = 1.0d0

    zero_freq = 0.0d0
    one_freq = 1.0d0
    freq_double = 1

  endif

  ;http://www.harrisgeospatial.com/docs/fft.html
  if ((image_xsize mod 2) eq 0) then begin
    image_xfreq = [zero_freq, image_xvect, image_xsize/2, -image_xsize/2 + image_xvect]/(image_xsize *image_xsamp_interval)
  endif else begin
    image_xfreq = [zero_freq, image_xvect, -(image_xsize/2 + 1) + image_xvect]/(image_xsize *  image_xsamp_interval)
  endelse



  if ((image_ysize mod 2) eq 0) then begin
    image_yfreq = [zero_freq, image_yvect, image_ysize/2, -image_ysize/2 + image_yvect]/(image_ysize *image_ysamp_interval)
  endif else begin
    image_yfreq = [zero_freq, image_yvect, -(image_ysize/2 + 1) + image_yvect]/(image_ysize *  image_ysamp_interval)
  endelse


  ; ************************************************
  ; Create frequency distance array   ; hive off into function
  ; ************************************************
  if freq_double eq 0 then begin
    fft_freq_dist = fltarr(image_xsize,image_ysize)
  endif else if freq_double eq 1 then begin
    fft_freq_dist = dblarr(image_xsize,image_ysize)
  endif

  for xx = 0, image_xsize-1, 1 do begin  ; must be more efficient way of doing this...
    for yy = 0, image_ysize-1, 1 do begin
      fft_freq_dist[xx,yy] = sqrt(image_xfreq[xx]^2 + image_yfreq[yy]^2)     ; needs to be checked
    endfor
  endfor

  ; ************************************************
  ; calculate filter
  ; ************************************************
  cut_freq = one_freq/cut_length_temp  ; 1/cut_length one_freq takes account of double/single precision

  kchoice = strcompress(strlowcase(kernel),/remove_all)
  pchoice = strcompress(strlowcase(pass),/remove_all)

  if pchoice eq 'high' || pchoice eq 'h' then begin
    if kchoice eq 'butterworth' || kchoice eq 'b' then begin
      taper_vals = fourier_highpass_butterworth(fft_freq_dist, cut_freq, butterworth_order, double = freq_double)
    endif else if kchoice eq 'gaussian' || kchoice eq 'gauss' || kchoice eq 'g' then begin
      taper_vals = fourier_highpass_gaussian(fft_freq_dist, cut_freq, double = freq_double)
    endif else if kchoice eq 'ideal' || kchoice eq 'i' then begin
      taper_vals = fourier_highpass_ideal(fft_freq_dist, cut_freq, double = freq_double)
    endif else begin
      print, ' error: the specified kernel: ' + string(kernel) +  ' is not valid. Please set the input variable kernel equal to "butterworth", "gaussian" or "ideal"'
      print, ' quitting...'
      stop
    endelse

  endif else if pchoice eq 'low' || pchoice eq 'l' then begin
    if kchoice eq 'gaussian' || kchoice eq 'gauss' || kchoice eq 'g' then begin
      taper_vals = fourier_lowpass_gaussian(fft_freq_dist, cut_freq, double = freq_double)
    endif else if kchoice eq 'butterworth' || kchoice eq 'b' then begin
      taper_vals = fourier_lowpass_butterworth(fft_freq_dist, cut_freq, butterworth_order, double = freq_double)
    endif else if kchoice eq 'ideal' || kchoice eq 'i' then begin
      taper_vals = fourier_lowpass_ideal(fft_freq_dist, cut_freq, double = freq_double)
    endif else begin
      print, ' error: the specified kernel: ' + string(kernel) +  ' is not valid. Please set the input variable kernel equal to "butterworth", "gaussian" or "ideal"'
      print, ' quitting...'
      stop
    endelse

  endif else begin
    print, ' error: the specified pass: ' + string(pass) +  ' is not valid. Please set the variable pass equal to "low" or "high" for lowpass or highpass filtering, respectively'
    print, ' quitting...'
    stop
  endelse

  ; ************************************************
  ; apply filter and transform to real space
  ; ************************************************
  fft_arr_tapered = fft_arr * taper_vals ; taper the fft array (appy filter)
  filtered_image_arr = real_part(fft(fft_arr_tapered, /inverse)) ; transform back to real space and discard complex parts

  if (nan_count gt 1) then begin
    if image_double eq 0 then filtered_image_arr[nan_list] = !values.f_nan else if image_double eq 1 then filtered_image_arr[nan_list] = !values.d_nan
  endif

  ; ************************************************
  ; write out filtered image
  ; ************************************************
  filtered_image_hdr = image_hdr
  if n_elements(filtered_image_path) ne 0 then writefits, filtered_image_path, filtered_image_arr, filtered_image_hdr

  ; ************************************************
  ; write out negative image
  ; ************************************************
  if n_elements(negative_image_path) ne 0 || n_elements(negative_image_hdr) ne 0 || n_elements(negative_image_arr) ne 0 then begin
    negative_image_hdr = image_hdr
    negative_image_arr = image_arr - filtered_image_arr
  endif
  if n_elements(negative_image_path) ne 0 then writefits, negative_image_path, negative_image_arr, negative_image_hdr

end
