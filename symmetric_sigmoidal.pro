;----------------------------------------------------------------------------------------
function symmetric_sigmoidal, x,a,b,c,d
;----------------------------------------------------------------------------------------
; Calculates the value of a symmetric sigmoidal function
;--(input)-------------------------------------------------------------------------------
; *** x           = the datapoint(s) to calculate for
; *** a,b,c and d = fitting parameters of the model
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  if n_params() lt 5 then begin ; check for parameters
   print, 'Syntax: symmetric_sigmoidal, x,a,b,c,d   '
   return, 'error'
  endif

  return, d + (a - d)/(1.0d0 + (x/c)^b)
end
