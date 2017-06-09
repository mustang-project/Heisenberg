;----------------------------------------------------------------------------------------
function fourier_highpass_gaussian, freq_dist, cut_freq, double = double
;----------------------------------------------------------------------------------------
; Function that returns a highpass Gaussian filter for Fourier filtering.
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

  taper_vals = fourier_lowpass_gaussian(freq_dist, cut_freq, double = double)

  if keyword_set(double) then value_1 = 1.0d0 else value_1 = 1.0e0

  taper_vals = value_1 - taper_vals

  return, taper_vals

end
