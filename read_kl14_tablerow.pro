pro read_kl14_tablerow, filename, name_vec, value_vec, de_log = de_log, compress_names = compress_names
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  ; de_log; [1] convert from log values [0] do nothing
  ; compress_names[1] remove white space from names_vec [0] do nothing
  
  split_string = ' '
  comment_string = '#'

  openr, r_lun, filename, /get_lun
  line_var = ''   ; dummy variable to hold the read line

  nlines = 0
  for ii = 0, 1, 1 do begin
    if ii eq 1 then begin ; create vectors on second run
      point_lun, r_lun, 0 ; rewind file to the beginning
      line_counter = -1 ; initialise line counter
      name_vec = strarr(nlines)
      value_vec = dblarr(nlines)
    endif
    while ~eof(r_lun) do begin
      readf, r_lun, line_var ; read line
      line_split = strsplit(line_var, split_string ,/extract) ;split line
      if n_elements(line_split) gt 1 && strmid(line_split[0],0,1) ne comment_string then begin ;if line is not empty or comment
        if ii eq 0 then nlines ++ else if ii eq 1 then begin ; initial run to count number of lines
          line_counter ++
          name_vec[line_counter] = line_split[0]
          value_vec[line_counter] = double(line_split[1])
          ; name_vec = [name_vec, line_split[0]]   ; store variable name
          ; value_vec = [value_vec, line_split[1]] ; store variable value
        endif
      endif
    endwhile
  endfor

  name_vec = strcompress(name_vec, /remove_all) ; remove any remaining spaces

  if de_log eq 1 then value_vec = 10.0^value_vec  ; this uses automatic conversion from strings to numbers
  if compress_names then name_vec = strcompress(name_vec, /remove_all) ; remove all whitespace from the names


  free_lun, r_lun


end
