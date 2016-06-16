;----------------------------------------------------------------------------------------
function mask_inside_circle, circ_x, circ_y, circ_r, test_x, test_y, no_edge_points = no_edge_points
;----------------------------------------------------------------------------------------
; *** inside_circle determines whether points with x-co-ordinates (test_x) and
; *** y-co-ordinates (test_y) are within a circle
;---input--------------------------------------------------------------------------------
; *** circ_x = The circle centre x-co-ordinate
; *** circ_y = The circle centre y-co-ordinate
; *** circ_r = The circle's radius
; *** test_x = The circle centre x-co-ordinate (must be of the same dimensions as test_y)
; *** test_y = The circle centre x-co-ordinate (must be of the same dimensions as test_x)
;---function returns---------------------------------------------------------------------
; *** returns -> an array of the same dimensions as test_x and test_y (or simply a number
; ***            if test_x and test_y are only numbers) with 0s in positions not within
; ***            the circle and 1s for points within the circle
;---keywords-----------------------------------------------------------------------------
; *** no_edge_points = (1/0) If set the routine returns 0 for points that lie exactly on
; ***                  the circle's perimiter. The normal behaviour is to return 1 for
; ***                  these points.
;----------------------------------------------------------------------------------------
  if array_equal(size(test_x),size(test_y)) ne 1 then begin
    print, " inside_circle error: test_x and test_y must have the same dimensions"
    print, " quitting..."
    stop
  endif

  test_r = sqrt( (test_x - circ_x)^2 + (test_y - circ_y)^2)

  if keyword_set(no_edge_points) then begin
   return, test_r lt circ_r  ; don't include edge_points
  endif else begin
    return, test_r le circ_r ; do include edge_points
  endelse




end
