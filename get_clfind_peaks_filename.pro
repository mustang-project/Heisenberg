;----------------------------------------------------------------------------------------
function get_clfind_peaks_filename, file
;----------------------------------------------------------------------------------------
; Simple routine to ensure that clfind filenames are unified between routines
;--(input)-------------------------------------------------------------------------------
; ***  file  = filename to append '_clfind.dat' to
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  return, strcompress(file+'_peaks.dat', /remove_all)
end

; get_clfind_peaks_filename
