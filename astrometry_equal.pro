;----------------------------------------------------------------------------------------
function astrometry_equal, imagea, hdra, imageb, hdrb, $ ; input parameters
         astr=astr, strict_equinox=strict_equinox        ; keywords
;----------------------------------------------------------------------------------------
; *** astrometry_equal checks if two astronomical images with associate fits headers have
; *** the same astrometry. Returns 1 if true and 0 if false.
;----------------------------------------------------------------------------------------
; *** To run, astrometry_equal  requires the:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; *** imagea : first image to be compared (array)
; *** imageb : second image to be compared (array)
; *** hdra   : fits header of first image (string array)
; *** hdrb   : fits header of second image (string array)
;--(keywords)----------------------------------------------------------------------------
; *** astr          : (1/0) if set hdra and hdrb are taken to be astrometry structures
; ***               : as returned by the extast astrometry routine
; ***               : http://idlastro.gsfc.nasa.gov/ftp/pro/astrom/extast.pro
; *** strict_equinox: (1/0) if set force equinox to be specified in exactly the same way
; ***               : else allows fk5/IRC to be equivilent to J2000 and fk4 to be equal
; ***               ; to B1950 ***(currently not implemented)***
;----------------------------------------------------------------------------------------
; *** Notes: 1) update with get_equinox to take account of fk4 and fk5/ircs
; ***           equivilence withe numerical values.
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  ; 1) first a quick size check
  if array_equal(size(imagea),size(imageb)) ne 1 then return, 0   ; check if each image has same size


  ; 2) now full check:
  if keyword_set(astr) then begin  ; hdrx keywords contain the astrometry structure and not the header
    astra = hdra
    astrb = hdrb
  endif else begin ; extract astrometry parameters for images from header
    if get_equinox(hdra) ne get_equinox(hdrb) then return, 0 ; check if each image has same equinox (e.g. J2000)
    extast, hdra, astra
    extast, hdrb, astrb
  endelse



  diff_list = compare_struct(astra,astrb)
  if n_elements(diff_list) gt 0 then for ii=0,(n_elements(diff_list)-1),1 do begin
    difference = diff_list[ii]
    if difference.field ne '.DATEOBS' && difference.field ne '.MJDOBS' && difference.field ne '' then return, 0 ; allow observation date to be different, otherwise astrometry must be the same
  endfor

  return, 1 ; astrometry is equal

end
