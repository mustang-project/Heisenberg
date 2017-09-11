;----------------------------------------------------------------------------------------
pro iteration_plot, plot_filename, iter_vec, errmin, errmax, ytitle, zero_ymin = zero_ymin
;----------------------------------------------------------------------------------------
;
;--(dependencies)------------------------------------------------------------------------
; *** To run, fourier_diffuse_fraction requires:
; *** 1) The coyote graphics library: http://www.idlcoyote.com/documents/programs.php
;--(input)-------------------------------------------------------------------------------
; *** plot_filename = the filename to output the plot image to
; *** iter_vec      = quanity to plot over the iteration
; *** errmax        = upwards error on the quanity being plotted
; *** errmin        = downwards error on the quanity being plotted
; *** ytitle        = title of the yaxis (name of the quanity being plotted with units)
;----------------------------------------------------------------------------------------

  n_els = n_elements(iter_vec)

  ; ensure error bars do not overflow the plot area
  if n_elements(errmax) gt 0 then ymax = 1.1 * max(iter_vec + errmax) else ymax = 1.1 * max(iter_vec)
  if (n_elements(zero_ymin) eq 1 && zero_ymin eq 1) then ymin = 0.0 else begin
    if n_elements(errmin) gt 0 then ymin = 0.9 * min(iter_vec - errmin) else ymin =  0.9 * min(iter_vec)
  endelse



  set_plot, 'PS' ; set device to PostScript
  device, filename=plot_filename, xsize=24,ysize=18,/color,bits_per_pixel=8,/encapsulated


    cgplot, indgen(n_els), iter_vec , err_yhigh = errmax, err_ylow = errmin, $
      xrange = [-0.2, n_els], yrange = [ymin, ymax], $
      xtitle = '!6Iteration', ytitle = ytitle , $
      font = -1, charsize = 1.3, charthick = 1.3, xthick = 1.3, ythick = 1.3, err_thick = 1.3, thick = 1.3

  device, /close ; close the PostScript file
  set_plot, 'X'  ; return plotting to the original device


end
