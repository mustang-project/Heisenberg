;----------------------------------------------------------------------------------------
function astrometry_equal, imagea, hdra, imageb, hdrb, tolerance $  ; input parameters
                                      , astr=astr, strict_equinox=strict_equinox ; keywords
;----------------------------------------------------------------------------------------
; *** astrometry_equal checks if two 2D astronomical images with associated fits headers
; *** have the same astrometry. Returns 1 if true and 0 if false.
;----------------------------------------------------------------------------------------
; *** To run, astrometry_equal  requires the:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(input)-------------------------------------------------------------------------------
; *** imagea    : first image to be compared (array)
; *** imageb    : second image to be compared (array)
; *** hdra      : fits header of first image (string array)
; *** hdrb      : fits header of second image (string array)
; *** tolerance : allowable tolerance in decimal degrees between astrometric position of
; ***             image pixels (e.g. 1.8" = 5e-4 decimal degrees)
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
  if array_equal(size(imagea,/dimensions),size(imageb,/dimensions)) ne 1 then return, 0   ; check if each image has same size

  ; 2) now full check:
  if keyword_set(astr) then begin  ; hdrx keywords contain the astrometry structure and not the header
    astra = hdra
    astrb = hdrb
  endif else begin ; extract astrometry parameters for images from header
    if get_equinox(hdra) ne get_equinox(hdrb) then return, 0 ; check if each image has same equinox (e.g. J2000)
    extast, hdra, astra
    extast, hdrb, astrb
  endelse

  ; 2a) get platescale:
  getrot, astra, rota, cdelt_a
  getrot, astrb, rotb, cdelt_b

  ; if cdelt_a le 0.0d0 && cdelt_b le 0.0d0 && abs(cdelt_a-cdelt_b) gt tolerance then return, 0 ; check platescales of images are within tolerance
  if  abs(cdelt_a[0]-cdelt_b[0]) gt tolerance then return, 0 ; check platescales of images are within tolerance
  if  abs(cdelt_a[1]-cdelt_b[1]) gt tolerance then return, 0 ; check platescales of images are within tolerance

  ; 2b) compare co-ordinates of all pixels in images (brute-force approach)

  image_dims = size(imagea, /dimensions) ; have checked size in 1)

  image_xdim = image_dims[0]
  image_ydim = image_dims[1]

  image_xcents = (fltarr(image_ydim)+1.0) ## (findgen(image_xdim) ) ; pixel centres in idl convention
  image_ycents = (fltarr(image_xdim)+1.0) # (findgen(image_ydim) )  ; idl pixel (0,0) = .fits pixel (1,1) ; need to take account of this for idl astronomy routines that convert to idl pixels

  xy2ad,image_xcents, image_ycents, astra, ra_a, dec_a ; get astrometric co-ordinates in degrees
  xy2ad,image_xcents, image_ycents, astrb, ra_b, dec_b

  ra_list = where( (abs(ra_a - ra_b) gt tolerance) eq 1, ra_count)
  if ra_count gt 0 then return, 0
  dec_list = where( (abs(dec_a - dec_b) gt tolerance) eq 1, dec_count)
  if dec_count gt 0 then return, 0

  return, 1 ; astrometry is sufficiently equal given the specified tolerance

end
