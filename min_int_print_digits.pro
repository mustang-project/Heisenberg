;----------------------------------------------------------------------------------------
function min_int_print_digits, number
;----------------------------------------------------------------------------------------
; calculate the minimum number of digits needed to print the integer part of a number
; or set of numbers
;--(input)-------------------------------------------------------------------------------
; *** number     = a number or set of numbers
;--(output)-----------------------------------------------------------------------------
; *** ndigits    = the number of digits needed to print the integer part of number
;----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  temp = max(number) ; get maximum of a set of numbers and protect input number
  ndigits = 1
  stopcrit = 0

  if temp lt 10 then return, ndigits

  ; divide number by 10 in a loop to get number of digits required
  while (stopcrit eq 0) do begin
    temp/=10
    if temp lt 10 then stopcrit =1
    ndigits ++
  endwhile
  return, ndigits

end
