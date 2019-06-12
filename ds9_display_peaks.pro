pro ds9_display_peaks, peakdir, x_file2short, window_name, x, y, no_disp = no_disp  ; , $
  ; logspacing_x,logrange_x, nlevels_x, npixmin, nsigma

  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  region_filename = strcompress(peakdir + x_file2short + '_peakid.reg',/remove_all)
  openw, peaks_lun, region_filename, width = 250, /get_lun
  printf, peaks_lun, '# Region file format: DS9 version 4.1'
  printf, peaks_lun, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
  printf, peaks_lun, 'image'

  if n_params() eq 5 then begin
    for ii = 0, n_elements(x)-1 do begin
      printf, peaks_lun, strcompress('point(' + string(x[ii]+1) + ',' + string(y[ii]+1) + ')', /remove_all) + ' # point=cross' ; convert idl to DS9 Co-ordinates with +1
    endfor
    free_lun, peaks_lun
  endif

  if n_elements(no_disp) ne 1 then begin
    image_path = strcompress(peakdir + x_file2short + '.fits', /remove_all )
    spawn_command = 'ds9 ' + image_path + ' -regions load ' + region_filename + ' -title ' + window_name + ' & '
    spawn, spawn_command
  endif


end
