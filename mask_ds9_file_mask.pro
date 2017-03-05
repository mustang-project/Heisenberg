;----------------------------------------------------------------------------------------
pro mask_ds9_file_mask, image, image_head, ds9_filepath, ds9_mask, negative=negative
;----------------------------------------------------------------------------------------
; *** mask_read_ds9_region_file reads in a ds9 region file with image co-ordinates and
; *** creates a idl mask array of 1s and 0s. Default behaviour is that defined regions
; *** are accepted (marked with 1s) set keyword /negative to mask regions (mark with 0s)
;----------------------------------------------------------------------------------------
; *** To run, the mask_tool requires the:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;---input--------------------------------------------------------------------------------
; *** image        = the image (as a 2D array) to be masked
; *** image_head   = the fits header of the image
; *** ds9_filepath = the filepath of the ds9 region file
;---output-------------------------------------------------------------------------------
; *** ds9_mask     = an array of the same size as image that describes the mask.
; ***              = 0s are blocked pixels and 1s are allowed pixels
;---keywords-----------------------------------------------------------------------------
; *** negative     = (1/0) if set (/negative) then a mask is created that blocks the
; ***                regions defined in the ds9 file. Default behavior is the block
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  ;***************************************
  ; 0) aquire astrometry
  ;***************************************
  image_equinox = get_equinox(image_head)
  ;***************************************
  ; 1) define mask and masking values
  ;***************************************

  ds9_mask = image ; create a mask image the same size as the image to be masked

  image_dims = size(image, /dimensions)

  image_xdim = image_dims[0]
  image_ydim = image_dims[1]

  image_xcents = (fltarr(image_ydim)+1.0) ## (findgen(image_xdim) ) + 1.0 ; pixel centres in .fits (and ds9 convention)
  image_ycents = (fltarr(image_xdim)+1.0) # (findgen(image_ydim) ) + 1.0  ; idl pixel (0,0) = .fits pixel (1,1) ; need to take account of this for idl astronomy routines that convert to idl pixels

  if keyword_set(negative) then begin ; "negative mask" blocks defined regions
    ds9_mask[*] = 1.0
    mask_val    = 0.0 ; value of mask at regions defined in ds9 mask
  endif else begin ; "positve mask" only unblocks defined regions
    ds9_mask[*] = 0.0
    mask_val    = 1.0 ; value of mask at regions defined in ds9 mask
  endelse

  ;***************************************
  ; 2) read in region file and perform a couple of checks
  ;***************************************
  openr, ds9_lun, ds9_filepath, /get_lun ; open the ds9 region file and get an unused lun

  linevar = '' ; define empty string variable to hold lines from reading

  while (~eof(ds9_lun)) do begin ; start while loop. End when end of file reached
    readf, ds9_lun, linevar
    linevar_trim = strtrim(linevar,2) ; remove trailing and leading blanks
    ;***************************************
    ; x) Split line where there are semi-colons and loop over (; = linebreak in ds9 region files)
    ;***************************************
    scolon_split = strtrim(strsplit(linevar_trim, ';',/extract),2) ; split line for semi-colon definitions and remove leading and trailing blanks

    for ss = 0, (n_elements(scolon_split)-1),1 do begin
      scolon_element = scolon_split[ss]

      scolon_elength = strlen(scolon_element)
      ;***************************************
      ; a) deal with comments and composite regions (must do this before splitting post-definition comments)
      ;***************************************
      if (strmid(scolon_element,0,1) eq '#') then begin
        if (scolon_elength  ge 12 && strmatch(strmid(scolon_element, 0, 12), '# composite ', /fold_case)) then begin ; if both of these are true we have a composite region. Use && to short-circuit so 2nd test is only evaluated if string is long enough
          print, ' warning: ds9 region file contains a composite region; composite regions should be handled, however, behavour of interpreter with respect to composite regions is not fully tested'
          continue ; continue as composite definition is unnecissary (this needs to be confirmed)
        endif else continue ; line is a comment so ignor line and continue
      endif

      ;***************************************
      ; b) deal with post-definition comments and composite markers
      ;***************************************
      comment_split = strsplit(scolon_element, '#',/extract) ; split line post-definition comments (started with a #)
      composite_split = strsplit(comment_split[0], "\|\|",/regex, /extract) ; remove '||' that comes after lines in a composite definition using regular expressions.

      scolon_def = strtrim(composite_split[0],2) ; extract first element and remove leading and trailing blanks
      scolon_dlength =  strlen(scolon_def)

      ;***************************************
      ; c) deal with global parameter lines
      ;***************************************
      if (scolon_dlength ge 7 && strmatch(strmid(scolon_def,0,7), 'global ',/fold_case)) then continue ; these lines specify irrelevant global parameters such as font style and can be skipped

      ;***************************************
      ; d) deal with unsupported mosaic images and multiple world co-ordinate systems
      ;***************************************
      if (strmatch(scolon_def, 'tile ?',/fold_case)) then begin ; check for unsupported mosaic image definitions
        print, "A definition including the form 'tile n' has been detected. ds9 mosaic images functionality is not supported. The line is:"
        print, linevar
        print, "the routine will now stop at code point 2d1"
        stop
      endif else if (strmatch(scolon_def, 'wcs?',/fold_case)) then begin ; check for unsupported multiple wcs definitions
        print, "A line definition including the form 'wcs*' has been detected. ds9 functionality for multiple wcs in a header is not supported. The line is:"
        print, linevar
        print, "the routine will now stop at code point 2d2"
        stop
      endif

      ;***************************************
      ; f) deal with lines specifying the co-ordinate system
      ;***************************************

      if (strmatch(scolon_def, 'image', /fold_case)) then begin  ; IMAGE  # pixel coords of current file ; only this co-ordinate system is currently suppported
        cmode = 'image' ; switch mode ; currently does nothing
        continue ; no need to do anything else so continue
      endif else if (strmatch(scolon_def, 'fk4', /fold_case) || strmatch(scolon_def, 'B1950', /fold_case) ) then begin ; FK4, B1950 # sky coordinate systems
        f_error,["The ds9 region file contains fk4/B1950 co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif else if (strmatch(scolon_def ,'fk5', /fold_case) || strmatch(scolon_def ,'J2000' , /fold_case)) then begin       ; FK5, J2000 # sky coordinate systems
        if (image_equinox eq 2000) then cmode = 'wcs' else begin
          f_error,["The ds9 region file contains definitions with an fk5/J2000 equinox, but the image header is not in this format.","This cannot be handled by the interpreter.","The line in this format is:",linevar]
        endelse
      endif else if (strmatch(scolon_def ,'IRCS', /fold_case)) then begin ; ICRS # currently same as J2000
        f_error,["The ds9 region file contains IRCS co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif else if (strmatch(scolon_def ,'GALACTIC', /fold_case)) then begin ; GALACTIC # sky coordinate systems
        f_error,["The ds9 region file contains GALACTIC co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif else if (strmatch(scolon_def ,'ECLIPTIC', /fold_case)) then begin ; ECLIPTIC  # sky coordinate systems
        f_error,["The ds9 region file contains ECLIPTIC co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif else if (strmatch(scolon_def ,'LINEAR', /fold_case)) then begin ; LINEAR # linear wcs as defined in file
        f_error,["The ds9 region file contains LINEAR co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif else if (strmatch(scolon_def ,'AMPLIFIER', /fold_case) || strmatch(scolon_def ,'DETECTOR' , /fold_case)) then begin       ; AMPLIFIER # mosaic coords of original file using ATM/ATV ; DETECTOR  # mosaic coords of original file usingDTM/DTV
        f_error,["The ds9 region file contains definitions a mosaic world coordinate system (either DETECTOR or AMPLIFIER).","Mosaic images are not currently supported by the interpreter.","The line in this format is:",linevar]
      endif else if (strmatch(scolon_def ,'PHYSICAL', /fold_case)) then begin   ; PHYSICAL  # pixel coords of original file using LTM/LTV
        f_error,["The ds9 region file contains PHYSICAL co-ordinate specifications.","Only IMAGE co-ordinates are supported.","Please manually convert in ds9, use the supplied script 'mask_ds9_file_convert' or the /convert keyword in mask_tool to convert the region file to image co-ordinates."]
      endif

      ;***************************************
      ; g) interpret region definitions
      ;***************************************
      line_split = strtrim(strsplit(scolon_def, '(,) ',/extract),2) ; split for spaces, brackets and commas and remove all leading and trailing spaces
      line_type = line_split[0]

      if (strmatch(line_split[0],'circle', /fold_case)) then  begin ; region format: circle x y radius
        ; line_split[1] = x, line_split[2] =y, line_split[3] = radius
        circle_test = mask_inside_circle(float(line_split[1]), float(line_split[2]), float(line_split[3]), image_xcents, image_ycents)
        mask_list = where(circle_test eq 1, mask_count) ; mask within circle (including edge of circle)
        if mask_count gt 0 then ds9_mask[mask_list] = mask_val
      endif else if (strmatch(line_split[0],'ellipse', /fold_case)) then begin  ; region format:  ellipse (x, y, maj, min, angle)
        if n_elements(line_split) eq 6 then begin
          ; line_split[1] = x, line_split[2] =y, line_split[3] = radius(x), line_split[4] = radius(y), , line_split[5] = angle
          ellipse_angle = float(line_split[5])*!dtor
          ellipse_test = mask_inside_ellipse(float(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4]), ellipse_angle, $
                                             image_xcents, image_ycents)
          mask_list = where(ellipse_test eq 1, mask_count) ; mask within circle (including edge of circle)
          if mask_count gt 0 then ds9_mask[mask_list] = mask_val
        endif else f_error,"Elliptical annuli are not currently supported by the interpreter"
      endif else if (strmatch(line_split[0],'box', /fold_case)) then begin  ; region format: box x y width height angle
        if n_elements(line_split) eq 6 then begin
          box_angle = float(line_split[5])*!dtor
          mask_ds9_box_vertices, line_split[1], line_split[2], line_split[3], line_split[4], box_angle, box_xpoints, box_ypoints
          ; now just apply polynomial masking method:
          poly_roi = obj_new('IDLanROI', box_xpoints, box_ypoints) ; create polynomial Region Of Interest object
          cpcentmask = poly_roi->containspoints(image_xcents,image_ycents)  ; returns these values: 0 = Exterior. The point lies strictly out of bounds of the ROI; 1 = Interior. The point lies strictly inside the bounds of the ROI, 2 = On edge. The point lies on an edge of the ROI boundary, 3 = On vertex. The point matches a vertex of the ROI
          mask_list = where(cpcentmask gt 0, mask_count)       ; mask interior and boundary points. Set null to allow a polygon that masks nothing. E.g. one defined exterior to the image region
          if mask_count gt 0 then ds9_mask[mask_list] = mask_val
          obj_destroy, poly_roi ; destroy object to prevent memory leakage is versions of IDL pre-8.0
        endif else f_error,"Box annuli are not currently supported by the interpreter"
      endif else if (strmatch(line_split[0],'polygon', /fold_case)) then begin  ; add in other region formats
        if cmode eq 'image' then begin
          poly_xx = float(line_split[1:*:2]) ; obtain all elements of line_split with odd subscripts and convert from string to float
          poly_yy = float(line_split[2:*:2]) ; obtain all elements of line_split with even subscripts, neglecting the 0th element, which is "polygon" and convert from string to float
        endif else if cmode eq 'wcs' then begin
          ra  = ten(line_split[1:*:2])*15
          dec = ten(line_split[2:*:2])
          adxy, image_head, ra, dec, poly_xx, poly_yy
        endif else f_error,"Image co-ordinate mode is not specified"
        if (n_elements(poly_xx) gt 2) then begin ; ensure region defined is not a line as this would create spurious masking.
          poly_roi = obj_new('IDLanROI', poly_xx, poly_yy) ; create polynomial Region Of Interest object
          cpcentmask = poly_roi->containspoints(image_xcents,image_ycents)  ; returns these values: 0 = Exterior. The point lies strictly out of bounds of the ROI; 1 = Interior. The point lies strictly inside the bounds of the ROI, 2 = On edge. The point lies on an edge of the ROI boundary, 3 = On vertex. The point matches a vertex of the ROI
          mask_list = where(cpcentmask gt 0, mask_count)       ; mask interior and boundary points. Set null to allow a polygon that masks nothing. E.g. one defined exterior to the image region
          if mask_count gt 0 then ds9_mask[mask_list] = mask_val
          obj_destroy, poly_roi ; destroy object to prevent memory leakage in versions of IDL pre-8.0
        endif
      endif else if (strmatch(line_split[0],'point', /fold_case) || strmatch(line_split[0],'line', /fold_case)  || strmatch(line_split[0],'vector', /fold_case) || strmatch(line_split[0],'text', /fold_case) || strmatch(line_split[0],'ruler', /fold_case)  || strmatch(line_split[0],'compass', /fold_case) || strmatch(line_split[0],'projection', /fold_case)) then begin
        continue ; lines, compasses, text etc. are not masked regions
      endif else if (strmatch(line_split[0],'annulus', /fold_case)) then begin
        f_error,['An annulus region is defined. This type of region is not yet handled by the interpreter. The line is:',linevar]
      endif else if (strmatch(line_split[0],'panda', /fold_case) || strmatch(line_split[0],'epanda', /fold_case) || strmatch(line_split[0],'bpanda', /fold_case)) then begin
        f_error,['A panda, epanda or bpanda region is defined. This type of region is not yet handled by the interpreter. The line is:',linevar]
      endif else begin
        f_error,["Something has gone wrong! A line in the ds9 region file has not been correctly interpreted. The line is:",linevar]
      endelse
    endfor ; end semi-colon loop
  endwhile  ; end loop over file lines

  close, ds9_lun ; close ds9 region file

end
