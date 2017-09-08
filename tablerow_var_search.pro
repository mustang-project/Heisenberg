function tablerow_var_search, var_name, name_vec, value_vec

  var_loc = where(name_vec eq var_name, var_count)
  if var_count eq 1 then return_var = value_vec[var_loc] else begin
   print, 'variable ' + var_name + ' not found' ; handle not foud variable
   stop
  endelse

  return, return_var

end
