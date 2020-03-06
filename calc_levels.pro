function calc_levels, loglevels, logrange_x, logspacing_x, nlinlevels_x

  if loglevels then begin ;set number of logarithmic contour levels for peak identification
      nlevels_x=logrange_x/logspacing_x+1 ;number of stellar contour levels
  endif else begin ;set number of linear contour levels for peak identification
      nlevels_x=nlinlevels_x ;number of stellar contour levels
  endelse

  return, nlevels_x

end
