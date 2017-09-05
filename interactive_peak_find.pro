pro interactive_peak_find, $
  loglevels, peakid_x,logspacing_x,logrange_x, nlevels_x, $
  peakdir, x_file2short, $
  flux_weight, npixmin, use_x2, $ ; inputs (x is gas or stars as appropriate)
  nsigma, sens_x, off_x, $
  ; ########################################################################
  x_peaks ; outputs

  if loglevels then begin
    maxval=max(alog10(peakid_x),/nan) ;maximum value
    maxlevel=(floor(maxval/logspacing_x)-1)*logspacing_x ;level below maximum value
    levels=10.^(maxlevel-logrange_x+dindgen(nlevels_x)/(nlevels_x-1)*logrange_x) ;array with levels
  endif else begin
    maxval=max(peakid_x,/nan)
    minval=min(peakid_x,/nan)
    levels=minval+(maxval-minval)*dindgen(nlevels_x)/(nlevels_x-1)
  endelse
  x_peaks=peaks2d(peakdir+x_file2short,levels=levels,/log,fluxweight=flux_weight,npixmin=npixmin-1) ;Nx4 array (N is # peaks) listing (x,y,F_tot,area) of peaks in pixel IDs
  sortpeaks=reverse(sort(x_peaks(*,2))) ;sort by decreasing intensity


  x_peaks(*,0)=x_peaks(sortpeaks,0)
  x_peaks(*,1)=x_peaks(sortpeaks,1)
  x_peaks(*,2)=x_peaks(sortpeaks,2)
  x_peaks(*,3)=x_peaks(sortpeaks,3)


  if use_x2 then i_x_max=n_elements(x_peaks(*,0))-1 else i_x_max=max(where(x_peaks(*,2) gt nsigma*sens_x+off_x)) ;last peak with intensity larger than the sensitivity limit
  x_peaks=x_peaks(0:i_x_max,*) ;reject stellar peaks with stellar flux lower than sensitivity limit









  ; if loglevels then begin
  ;   maxval=max(alog10(peakidstar),/nan) ;maximum value
  ;   maxlevel=(floor(maxval/logspacing_s)-1)*logspacing_s ;level below maximum value
  ;     levels=10.^(maxlevel-logrange_s+dindgen(nlevels_s)/(nlevels_s-1)*logrange_s) ;array with levels
  ; endif else begin
  ;     maxval=max(peakidstar,/nan)
  ;     minval=min(peakidstar,/nan)
  ;     levels=minval+(maxval-minval)*dindgen(nlevels_s)/(nlevels_s-1)
  ; endelse
  ; starpeaks=peaks2d(peakdir+starfile2short,levels=levels,/log,fluxweight=flux_weight,npixmin=npixmin-1) ;Nx4 array (N is # peaks) listing (x,y,F_tot,area) of peaks in pixel IDs
  ; sortpeaks=reverse(sort(starpeaks(*,2))) ;sort by decreasing intensity
  ; starpeaks(*,0)=starpeaks(sortpeaks,0)
  ; starpeaks(*,1)=starpeaks(sortpeaks,1)
  ; starpeaks(*,2)=starpeaks(sortpeaks,2)
  ; starpeaks(*,3)=starpeaks(sortpeaks,3)
  ; if use_star2 then istarmax=n_elements(starpeaks(*,0))-1 else istarmax=max(where(starpeaks(*,2) gt nsigma*sensstar+offstar)) ;last peak with intensity larger than the sensitivity limit
  ; starpeaks=starpeaks(0:istarmax,*) ;reject stellar peaks with stellar flux lower than sensitivity limit




end
