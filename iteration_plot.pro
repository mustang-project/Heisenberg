;----------------------------------------------------------------------------------------
pro iteration_plot, plot_filename, iter_vec, errmin, errmax, ytitle
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

  set_plot, 'PS' ; set device to PostScript
  device, filename=plot_filename, /landscape ; Use device to set some PostScript device options

    cgplot, indgen(n_els), iter_vec, err_yhigh = errmax, err_ylow = errmin, $
      xrange = [-0.1, n_els], $
      xtitle = 'Iteration', ytitle = ytitle

  device, /close ; close the PostScript file
  set_plot, 'X'  ; return plotting to the original device

end
