function peak_var_read, var_name, original_value, global_var = global_var
  ; var_name : string name of the variable (e.g. 'logspacing_s')

  inp_break = 0
  while (inp_break eq 0) do begin
    print, "Do you want to change " + var_name + " ?"
    print, "Currently " +  var_name  + " = ", original_value
    if n_elements(global_var) eq 1 && global_var eq 1 then print, "Warning " + var_name + " affects peak identification in both maps"
    print, "Type 'No' to leave the variable unaltered, or enter the new value for the variable"
    var_in =''
    read, var_in

    var_in = strcompress(var_in, /remove_all)

    if strmatch(var_in, 'no', /fold_case) || strmatch(var_in, 'n', /fold_case)  then begin
      var_out = original_value
      inp_break = 1
    endif else if valid_num(var_in, val) eq 1 then begin
      var_out = val
      inp_break = 1
    endif else begin
      print, "the input: " + var_in + " . Is not in the correct format. Please enter a valid number or the word 'No'."
    endelse

  endwhile

  return, var_out

end
