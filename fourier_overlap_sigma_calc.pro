;----------------------------------------------------------------------------------------
function fourier_overlap_sigma_calc, x, a, b, c, d, covariance_matrix, sigma_x
;----------------------------------------------------------------------------------------
; calculates the value of the 1 sigma unncertainty on the correction fractor, due to
; overlapping regions, applied to the compact emission fraction of an image
;--(input)-------------------------------------------------------------------------------
; *** x                  = the datapoint to calculate for
; *** a,b,c and d        = fitting parameters of the model (for indices a=0, b=1, c=2
; ***                      and d=3)
; *** covariance_matrix = the covariance matrix of the fitting parameters
; *** sigma_x           = the 1 sigma uncertainty on x
;--(output)-----------------------------------------------------------------------------
; *** sigma_final       = the 1 sigma uncertainty on the the compact fraction corrective
; ***                     factor to account for overlapping peaks
;----------------------------------------------------------------------------------------


  compile_opt idl2, strictarrsubs

  if n_params() lt 7 then begin ; check for parameters
   print, 'Syntax: fourier_overlap_sigma_calc, x, a, b, c, d, covariance_matrix'
   return, 'error'
  endif

  if n_elements(x) gt 1 then f_error,['fourier_overlap_sigma_calc: the input variable x must have only one element']

  n_variables = 5

  ; calculate the error on each parameter
  sigma_a = sqrt(covariance_matrix[0,0])
  sigma_b = sqrt(covariance_matrix[1,1])
  sigma_c = sqrt(covariance_matrix[2,2])
  sigma_d = sqrt(covariance_matrix[3,3])
  sigma_vec = [sigma_a, sigma_b, sigma_c, sigma_d, sigma_x]

 ; partial derivatives of a symmetric sigmoidal for each parameter
  aderiv = 1.0d0/(((x/c)^b) + 1.0d0)
  bderiv =  -((a - d) * ((x/c)^b) * alog(x/c))/((((x/c)^b) + 1.0d0)^2)
  cderiv = (b*(a - d) * (x/c)^b)/(c *(((x/c)^b) + 1.0d0)^2)
  dderiv = 1.0d0 - (1.0d0/(((x/c)^b) + 1.0d0))
  xderiv = (b*(d-a)*(x/c)^b)/(x*(((x/c)^b)+1)^2)
  differentials_vec = [aderiv, bderiv, cderiv, dderiv, xderiv]



  ; make the covariance_array
  covariance_arr = fltarr(n_variables,n_variables)
  covariance_arr[0:-2,0:-2] = covariance_matrix ; covariance of a,b,c,d
  covariance_arr[-1,-1] = sigma_x^2 ; sigma^2 to match rest of the matrix. The error in metallicity is correlated with itself and nothing else



  ;hughes and hase (2010) calculation
  sigma_parts_arr = fltarr(n_variables, n_variables)
  for ii = 0, n_variables - 1, 1 do begin
    for jj = 0, n_variables - 1, 1 do begin
      ; eliminate unneccesary step
      ; coorelation_coefficient = covariance_arr[ii,jj]/(sigma_vec[ii] * sigma_vec[jj])
      ; sigma_parts_arr[ii,jj] = differentials_vec[ii] * differentials_vec[jj] * coorelation_coefficient * sigma_vec[ii] * sigma_vec[jj]
      sigma_parts_arr[ii,jj] = differentials_vec[ii] * differentials_vec[jj] * covariance_arr[ii,jj]
    endfor
  endfor
  sigma_final = sqrt(total(sigma_parts_arr))

  return, sigma_final

end
