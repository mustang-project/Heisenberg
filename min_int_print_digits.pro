function min_int_print_digits, number
  compile_opt idl2, strictarrsubs

  temp = max(number) ; get maximum of a set of numbers and protect input number
  ndigits = 1
  stopcrit = 0

  if temp lt 10 then return, ndigits

  while (stopcrit eq 0) do begin
    temp/=10
    if temp lt 10 then stopcrit =1
    ndigits ++
    ; print, b, num
  endwhile
  return, ndigits

end
