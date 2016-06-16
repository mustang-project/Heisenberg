;----------------------------------------------------------------------------------------
pro mask_tool, image_input, image_hdr $ ; image variables
            , ds9_positive_path = ds9_positive_path, ds9_negative_path = ds9_negative_path $ ; ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
            , masked_image_path = masked_image_path, mask_path = mask_path $ ; path for the image with masked regions set to nans to be output to
            , image_masked = image_masked, mask_array=mask_array $ ;
            , header_output = header_output $ ; output the image header (useful for single parameter calling and with keyword image_masked)
            , convert = convert $ ; routine behaviour control keywords
            , conv_filepath = conv_filepath ;
;--(dependencies)------------------------------------------------------------------------
; *** To run, the mask_tool  requires:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;--(limitations)------------------------------------------------------------------------
; *** 1) if the ds9 region file has co-ordinates and distances that are not image pixels
; ***    then it is assumed that .fits header keyword 'cdelt' is in degrees per pixel
;--(Optional dependencies)---------------------------------------------------------------
; *** 1) if keyword "convert" is set ds9 is required: http://ds9.si.edu . Otherwise ds9
; ***    is not required.
;--(input)-------------------------------------------------------------------------------
; ***  image_input = the file path of the image to be masked e.g './images/image.fits'
; ***  image_hdr   = (optional) The header of the image. If this paramter is suppied then
; ***                image_input is the image in an array and not a filepath
;--(keywords)-----------------------------------------------------------------------------
; *** Note: one or both of ds9_positive_path and ds9_negative_path must be specified.
; *** ds9_positive_path = the file path for the ds9 region file specifying allowed regions
; ***                     of the image.
; *** ds9_negative_path = the file path for the ds9 region file specifying blocked regions
; ***                     of the image (regions to be masked).
; *** masked_image_path = filepath to output the masked image to
; *** mask_path         = filepath to output the mask as a fits file with 0-valued pixels
; ***                     being pixels that are masked and pixels with value 1 being
; ***                     allowed pixels
; *** image_masked      = variable to output masked image as an array to
; *** mask_array        = variable to output masked image as an array to. Has convention
; ***                     as descirbed in mask_path
; *** header_output     = Output the image header (useful for single parameter calling
; ***                     and with keyword image_masked)
; *** convert           = (1/0) if keyword is set convert ds9 region files to image
; ***                     co-ordinates (requries ds9). The routine will "convert" files
; ***                     already in the image format, but the resulting file should be
; ***                     exactly the same as the original. This also requires a single
; ***                     paramter call, as the filepath to the file is required to allow
; ***                     ds9 to
; *** conv_path         = filepath to an image used to convert the ds9 filepath. Used for
; ***                     two parameter calls where the filepath is not specified. Should
; ***                     be an image with the same astrometry as the image that is being
; ***                     masked
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  ; ***********************************
  ; 1) read in image file and/or convert ds9 files if neccessary
  ; ***********************************
  if (n_params() eq 1) then image = readfits(image_input, image_head) else begin ; read in image
    image = image_input ; image and hdr supplied to command
    image_head = image_hdr
  endelse

  if keyword_set(convert) then begin  ; if neccessary convert file to image co-ordinates using ds9
    if (n_params() eq 1) then begin
      conv_filepath = image_input ; for single parameter calls use specified image path
    endif else if n_elements(conv_path) ne 0 then begin
      conv_filepath = conv_path  ; if two parameter call need keyword conv_path to specify the image for ds9 to use to convert
    endif else begin
      print, "for mask_tool to convert a ds9 region file to image co-ordinates you must either specify a path to the file with a single parameter call or use the keyword conv_path"
      stop
    endelse

    if n_elements(ds9_positive_path) ne 0 then begin
      positive_outpath = ds9_positive_path + '_maskconverted'
      mask_ds9_file_convert, conv_filepath, ds9_positive_path, positive_outpath
      ds9_positive_path = positive_outpath
    endif
    if n_elements(ds9_negative_path) ne 0 then begin
      negative_outpath = ds9_negative_path + '_maskconverted'
      mask_ds9_file_convert, conv_filepath, ds9_negative_path, negative_outpath
      ds9_negative_path = negative_outpath
    endif
  endif

  ; ***********************************
  ; 2) interpret ds9 file
  ; ***********************************
  if (n_elements(ds9_positive_path) ne 0 and n_elements(ds9_negative_path) ne 0) then begin ; both positive and negative masks supplied
    mask_ds9_file_mask, image, image_head, ds9_positive_path, positive_mask ; allowed regions
    mask_ds9_file_mask, image, image_head, ds9_negative_path, negative_mask, /negative ; not allowed regions
    mask = positive_mask * negative_mask ; only regions allowed by both negative and positive mask allowed in final mask

  endif else if (n_elements(ds9_positive_path) ne 0) then begin ; allowed regions ; only positive mask supplied
    mask_ds9_file_mask, image, image_head, ds9_positive_path, mask ; create mask straight from subroutine

  endif else if (n_elements(ds9_negative_path) ne 0) then begin ; only negative mask supplied
    mask_ds9_file_mask, image, image_head, ds9_negative_path, mask, /negative ; not allowed regions ; create mask straight from subroutine

  endif else begin
    print, ' error: neither a positive or negative ds9 mask filepath has been supplied -- you must supply one or both'
    print, ' quitting...'
    stop
  endelse

  ; ***********************************
  ; 3) Output
  ; ***********************************
  mask_list = where(mask eq 0.0, mask_count)
  image_masked = image
  header_output = image_head
  if mask_count gt 0 then image_masked[mask_list] = !values.f_nan ; apply mask

  if n_elements(masked_path) ne 0 then writefits, mask_path, mask, image_head ; write a out the mask to .fits file
  if n_elements(masked_image_path) ne 0 then writefits, masked_image_path, image_masked, image_head ; write a masked (masked regions set to nans) image to .fits file

  mask_array = mask

end
