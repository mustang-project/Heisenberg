;----------------------------------------------------------------------------------------
function fourier_lowpass_butterworth, freq_dist, cut_freq, butterworth_order, double = double
;----------------------------------------------------------------------------------------
; Function that returns a lowpass Butterworth filter for Fourier filtering.
;--(input)-------------------------------------------------------------------------------
; ***  freq_dist = An array of frequency distances
; ***  cut_freq  = the frequency at which to cut from the frequency distribution
;--(keywords)----------------------------------------------------------------------------
; ***  double    = (1/0) 1: taper_vals is returned as a double precision array else
; ***              0: (default) taper_vals is returned as a single precision (float) array
;--(output)------------------------------------------------------------------------------
; *** taper_vals = array of same dimensions as freq_dist, with which to multiply the
; ***              fourier spectrum of an image to apply the filter
;-----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  ;check that butterworth order is a positive interger greter than 1
  if butterworth_order lt 1 then begin
    if ~is_an_integer(butterworth_order) then begin
      if (round(butterworth_order) - butterworth_order) ne 0 then begin
        print, 'error: butterworth_order must be a positve integer greater than or equal to 1, not:', butterworth_order
        stop
      endif
    endif
  endif

  ; calculate the filter
  if keyword_set(double) then value_1 = 1.0d0 else value_1 = 1.0e0
  taper_vals = value_1/(value_1+((freq_dist / cut_freq)^(2*butterworth_order)))

  return, taper_vals

end
