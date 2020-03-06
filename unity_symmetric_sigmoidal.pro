;----------------------------------------------------------------------------------------
function unity_symmetric_sigmoidal, x,a,b,c
;----------------------------------------------------------------------------------------
; Calculates the value of a symmetric sigmoidal function that tends to 1 as x
; increases
;--(input)-------------------------------------------------------------------------------
; *** x     = the datapoint(s) to calculate for
; *** a,b,c = fitting parameters of the model
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs

  if n_params() lt 4 then begin ; check for parameters
   print, 'Syntax: unity_symmetric_sigmoidal, x,a,b,c  '
   return, 'error'
  endif

  d = 1.0d0
  return, d + (a - d)/(1.0d0 + (x/c)^b)
end
