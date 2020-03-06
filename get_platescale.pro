;----------------------------------------------------------------------------------------
function get_platescale, hdr_or_astr, astr_tolerance, $
  xplatescale = xplatescale, yplatescale = yplatescale
;----------------------------------------------------------------------------------------
; get the platescale (angular size of pixels) from a .fits header
;--(input)-------------------------------------------------------------------------------
; *** hdr_astr       = a fits header or astrometry structure (extracted from a hdr with
; ***                  extast)
; *** astr_tolerance = allowable tolerance in decimal degrees between astrometric
; ***                  position of image pixels (e.g. 1.8" = 5e-4 decimal degrees)
;--(keywords)----------------------------------------------------------------------------
; *** xplatescale    = return the plate scale in the "NAXIS1" dimension
; *** yplatescale    = return the plate scale in the "NAXIS2" dimension
;--(output)------------------------------------------------------------------------------
; *** cdelt          = the platescale in decimal degrees
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  if n_params() lt 2 then begin ; check for parameters
   print, 'Syntax: get_platescale(hdr_or_astr, astr_tolerance   '
   return, 'ERROR'
  endif


  getrot, hdr_or_astr, rotation, cdelt_var ; get cdelt value ; only single precision

  if n_elements(xplatescale) eq 1 && xplatescale eq 1 then begin
    cdelt = cdelt_var[0]
  endif else if n_elements(yplatescale) eq 1 && yplatescale eq 1 then begin
    cdelt = cdelt_var[1]
  endif else begin
    if n_elements(cdelt_var) eq 2 && abs(abs(cdelt_var[0]) - abs(cdelt_var[1])) le astr_tolerance then begin
      cdelt = mean(abs(cdelt_var))
    endif else begin ; check that the image is square
      print, "image not square"
      stop
    endelse
  endelse

  return, cdelt


end
