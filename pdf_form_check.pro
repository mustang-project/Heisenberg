;----------------------------------------------------------------------------------------
function pdf_form_check, cumulative_pdf
;----------------------------------------------------------------------------------------
; checks to see that cumulative pdf has been generated with the expected properties to
; avoid errors arising from spurious interpolation that results from a poor pdf
;--(input)-------------------------------------------------------------------------------
; *** cumulative_pdf  = the cumulative pdf to be tested
;--(output)------------------------------------------------------------------------------
; *** returns 1 if the pdf passes the checks, 0 otherwise
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs


  ; make a series of checks in order of computational effort:
  if cumulative_pdf[-1] ne 1.0d0 then return, 0 ; the final value in the pdf should be 1
  if total(long(cumulative_pdf lt 0.0d0)) ge 1 then return, 0 ; there should be no negative values in the cumulative pdf
  if total(long((cumulative_pdf[1:-1] - cumulative_pdf[0:-2]) lt 0.0d0)) ge 1 then return, 0 ; the pdf should be monotonically increasing

  ; if the pdf passes all checks then return 1
  return,1

end
