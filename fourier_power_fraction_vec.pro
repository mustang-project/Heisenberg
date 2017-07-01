;----------------------------------------------------------------------------------------
pro fourier_power_fraction_vec, cut_vec $ ; vector of cut_lengths
                              , image_input, image_hdr $ ; input image
                              , partial_power_frac = partial_power_frac $
                              , partial_power = partial_power, total_power = total_power
;----------------------------------------------------------------------------------------
;
;--(dependencies)------------------------------------------------------------------------
; *** To run, fourier_power_fraction requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; ***
;--(output)------------------------------------------------------------------------------
; ***
;--(keywords)-----------------------------------------------------------------------------
; ***
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  ; ************************************************
  ; test input validity
  ; ************************************************
  if n_params() lt 2 then begin
     print,'Syntax: fourier_power_fraction, cut_vec, image_input, [image_hdr'
     print,'                                partial_power_frac = , partial_power ='
     print,'                                total_power ='
     return
  endif else if (n_params() eq 2) then image_arr = readfits(image_input, image_hdr) else begin ; read in image
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
  powerspectrum_arr = abs(fft_arr)^2     ; calculate the powerspectrum for the original image
  total_power = total(powerspectrum_arr) ; total power in the original image



  ; ************************************************
  ; Create frequency vectors !CHECK!
  ; ************************************************

  image_dimensions = size(image_arr, /dimensions)
  image_xsize      = image_dimensions[0]
  image_ysize      = image_dimensions[1]

  if size(cut_vec, /type) eq 4 then begin
    image_xvect = findgen((image_xsize - 1)/2) + 1
    image_yvect = findgen((image_ysize - 1)/2) + 1

    if ~keyword_set(image_xsamp_interval) then image_xsamp_interval = 1.0e0
    if ~keyword_set(image_ysamp_interval) then image_ysamp_interval = 1.0e0

    zero_freq = 0.0e0
    one_freq = 1.0e0
    freq_double = 0

  endif else if size(cut_vec, /type) eq 5 then begin
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
  ; Create frequency distance array   ; hive off into function !CHECK!
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
  ; iterate over the cut vector:
  ; ************************************************
  if image_double eq 1 then partial_power = dblarr(n_elements(cut_vec)) else partial_power = fltarr(n_elements(cut_vec))
  ; filtered_powersepc = powerspectrum_arr
  for ii = 0, n_elements(cut_vec)-1, 1 do begin
    cut_freq = one_freq/cut_vec[ii]  ; 1/cut_length one_freq takes account of double/single precision
    ; ************************************************
    ; Use ideal filter to filter the powerspectrum
    ; ************************************************
    taper_vals = fourier_highpass_ideal(fft_freq_dist, cut_freq, double = freq_double)
    filtered_powersepc = powerspectrum_arr * taper_vals
    partial_power[ii] = total(filtered_powersepc)

  endfor

  partial_power_frac = partial_power/total_power

end
