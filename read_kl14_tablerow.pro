pro read_kl14_tablerow, filename, name_vec, value_vec, de_log = de_log
  ; de_log; [1] convert from log values [0] do nothing
  ; ; NOT pre-8.0 compatible


  split_string = ' '
  comment_string = '#'

  openr, r_lun, filename, /get_lun
  line_var = ''   ; dummy variable to hold the read line
  name_vec = []  ; NOT pre-8.0 compatible
  value_vec = [] ; NOT pre-8.0 compatible

  while ~eof(r_lun) do begin
    readf, r_lun, line_var ; read line
    line_split = strsplit(line_var, split_string ,/extract) ;split line
    if n_elements(line_split) gt 1 && strmid(line_split[0],0,1) ne comment_string then begin ;if line is not empty or comment
      name_vec = [name_vec, line_split[0]]   ; store variable name
      value_vec = [value_vec, line_split[1]] ; store variable value
    endif
  endwhile

  name_vec = strcompress(name_vec, /remove_all) ; remove any remaining spaces

  if de_log eq 1 then value_vec = 10.0^value_vec  ; this uses automatic conversion from strings to numbers

  free_lun, r_lun


end
