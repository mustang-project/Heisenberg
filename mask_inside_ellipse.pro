;----------------------------------------------------------------------------------------
function mask_inside_ellipse, ellipse_x, ellipse_y, ellipse_maj, ellipse_min, ellipse_ang $
                      , test_x, test_y $
                      , no_edge_points = no_edge_points
;----------------------------------------------------------------------------------------
; *** inside_ellipse determines whether points with x-co-ordinates (test_x) and
; *** y-co-ordinates (test_y) are within an ellipse
;---input--------------------------------------------------------------------------------
; *** ellipse_x   = The ellipse centre x-co-ordinate
; *** ellipse_y   = The ellipse centre y-co-ordinate
; *** ellipse_maj = The ellipse's major-axis
; *** ellipse_min = The ellipse's minor-axis
; *** box_angle   = angle of rotation of the ellipse with respect to the image in radians
; ***               (requires multiplcation by !dtor from ds9 definition)
; *** test_x = The ellipse centre x-co-ordinate (must be of the same dimensions as test_y)
; *** test_y = The ellipse centre x-co-ordinate (must be of the same dimensions as test_x)
;---function returns---------------------------------------------------------------------
; *** returns -> an array of the same dimensions as test_x and test_y (or simply a number
; ***            if test_x and test_y are only numbers) with 0s in positions not within
; ***            the ellipse and 1s for points within the ellipse
;---keywords-----------------------------------------------------------------------------
; *** no_edge_points = (1/0) If set the routine returns 0 for points that lie exactly on
; ***                  the ellipse's perimiter. The normal behaviour is to return 1 for
; ***                  these points.
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  if array_equal(size(test_x),size(test_y)) ne 1 then f_error,"inside_ellipse error: test_x and test_y must have the same dimensions"

  test_dist =((((cos(ellipse_ang) * (test_x -ellipse_x)) + (sin(ellipse_ang) * (test_y -ellipse_y)))^2)/(ellipse_maj*ellipse_maj)) + $
             ((((sin(ellipse_ang) * (test_x -ellipse_x)) - (cos(ellipse_ang) * (test_y -ellipse_y)))^2)/(ellipse_min*ellipse_min))

   if keyword_set(no_edge_points) then begin
    return, test_dist lt 1.0 ; don't include edge_points
   endif else begin
     return, test_dist le 1.0 ; do include edge_points
   endelse

end
