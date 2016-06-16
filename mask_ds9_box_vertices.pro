;----------------------------------------------------------------------------------------
pro mask_ds9_box_vertices, box_x, box_y, box_width, box_height, box_angle $ ; inputs
                    , box_xpoints, box_ypoints ; output
;----------------------------------------------------------------------------------------
; *** ds9_box_vertices
;---input--------------------------------------------------------------------------------
; *** box_x       = x co-ordinate of the centre of the box (from ds9 definition)
; *** box_y       = y co-ordinate of the centre of the box (from ds9 definition)
; *** box_width   = width of the box (from ds9 definition)
; *** box_height  = height of the box (from ds9 definition)
; *** box_angle   = angle of rotation of the box with respect to the image in radians
; ***               (requires multiplcation by !dtor from ds9 definition)
;---output-------------------------------------------------------------------------------
; *** box_xpoints = ordered x co-ordinates of the box (in format neccessary to create
; ***               an 'IDLanROI' object
; *** box_ypoints = ordered y co-ordinates of the box (in format neccessary to create
; ***               an 'IDLanROI' object
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults 

  ; 1) calculate half-width and half-height
  box_wr = box_width/2.0
  box_hr = box_height/2.0

  ; 2) calculate vertices of the box
  ; top-left from box frame of reference
  xx_tl = box_x + (- box_hr*sin(box_angle)) + (-box_wr*cos(box_angle))
  yy_tl = box_y + (box_hr*cos(box_angle))   + (-box_wr*sin(box_angle))
  ; top-right from box frame of reference
  xx_tr = box_x + (-box_hr*sin(box_angle))  +  (box_wr*cos(box_angle))
  yy_tr = box_y + (box_hr*cos(box_angle))   +  (box_wr*sin(box_angle))
  ; bottom-right from box frame of reference
  xx_br = box_x + (box_hr*sin(box_angle))   +  (box_wr*cos(box_angle))
  yy_br = box_y + (-box_hr*cos(box_angle))  +  (box_wr*sin(box_angle))
  ; bottom-left from box frame of reference
  xx_bl = box_x + (box_hr*sin(box_angle))   + (-box_wr*cos(box_angle))
  yy_bl = box_y + (-box_hr*cos(box_angle))  +  (-box_wr*sin(box_angle))

  ; 3) create vectors of box x and y co-ordinates
  box_xpoints = [xx_tl, xx_tr, xx_br, xx_bl]
  box_ypoints = [yy_tl, yy_tr, yy_br, yy_bl]




end
