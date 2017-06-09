;----------------------------------------------------------------------------------------
function is_an_integer, variable
;----------------------------------------------------------------------------------------
; replicates isa(var, /integer, /scalar) for idl versions pre-8.0, i.e. returns 1 if
; variable is a single integer and 0 otherwise
;--(input)-------------------------------------------------------------------------------
; *** variable  = variable to be tested
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  if n_elements(variable) ne 1 then return, 0

  tn = size(variable, /type)

  ; 2 = int, 3 = long, 12 = uint, 13 = ulong, 14 = long64, 15 = ulong64
  if tn eq 2 || tn eq 3 || tn eq 12 || tn eq 13 || tn eq 14 || tn eq 15 then return, 1 else return, 0

end
