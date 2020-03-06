;----------------------------------------------------------------------------------------
function fourier_lowpass_gaussian, freq_dist, cut_freq, double = double
;----------------------------------------------------------------------------------------
; Function that returns a lowpass Gaussian filter for Fourier filtering.
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

  if keyword_set(double) then taper_vals = exp( - ((freq_dist*freq_dist)  / (2.0d0 * (cut_freq*cut_freq))) ) $
  else taper_vals = exp( - ((freq_dist*freq_dist)  / (2.0e0 * (cut_freq*cut_freq))) )

  return, taper_vals

end
