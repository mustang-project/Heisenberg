;----------------------------------------------------------------------------------------
function output_var_search, var_name, name_vec, value_vec
;----------------------------------------------------------------------------------------
; returns the value of the variable "var_name" given the outputted vectors
; (name_vec and value_vec) from read_heisenberg_output.pro
;--(input)-------------------------------------------------------------------------------
; *** var_name          = the name of the variable to get a value for
; *** name_vec          = the vector of variable names
; *** value_vec         = the vector of variable values corresponding to
; ***                     name_vec
;--(output)-----------------------------------------------------------------------------
; *** return_var        = the value of valiable "var_name"
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  var_loc = where(name_vec eq var_name, var_count)
  if var_count eq 1 then return_var = value_vec[var_loc] else begin
   print, 'variable ' + var_name + ' not found' ; handle not found variable
   stop
  endelse

  return, return_var

end
