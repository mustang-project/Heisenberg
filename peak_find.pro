function peak_find, peak_find_tui, $
  loglevels, peakid_x,logspacing_x,logrange_x, nlevels_x, $
  peakdir, x_file2short, $
  flux_weight, npixmin, use_x2, $ ; inputs (x is gas or stars as appropriate)
  nsigma, sens_x, off_x
  ; ########################################################################
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  if loglevels eq 1 then begin
    maxval=max(alog10(peakid_x),/nan) ;maximum value
    maxlevel=(floor(maxval/logspacing_x)-1)*logspacing_x ;level below maximum value
    levels=10.^(maxlevel-logrange_x+dindgen(nlevels_x)/(nlevels_x-1)*logrange_x) ;array with levels
  endif else begin
    maxval=max(peakid_x,/nan)
    minval=min(peakid_x,/nan)
    levels=minval+(maxval-minval)*dindgen(nlevels_x)/(nlevels_x-1)
  endelse
  x_peaks=peaks2d(peakdir+x_file2short,levels=levels,/log,fluxweight=flux_weight,npixmin=npixmin-1,peak_find_tui=peak_find_tui) ;Nx4 array (N is # peaks) listing (x,y,F_tot,area) of peaks in pixel IDs

  if (size(x_peaks,/type) eq 7) then begin ; catch errors. 7 = string
    if x_peaks eq 'TOOFEWLINES' then return, 'TOOFEWLINES'
  endif

  sortpeaks=reverse(sort(x_peaks[*,2])) ;sort by decreasing intensity
  x_peaks[*,0]=x_peaks[sortpeaks,0]
  x_peaks[*,1]=x_peaks[sortpeaks,1]
  x_peaks[*,2]=x_peaks[sortpeaks,2]
  x_peaks[*,3]=x_peaks[sortpeaks,3]

  if use_x2 then i_x_max=n_elements(x_peaks[*,0])-1 else i_x_max=max(where(x_peaks[*,2] gt nsigma*sens_x+off_x)) ;last peak with intensity larger than the sensitivity limit
  x_peaks=x_peaks[0:i_x_max,*] ;reject stellar peaks with stellar flux lower than sensitivity limit

  return, x_peaks


end
