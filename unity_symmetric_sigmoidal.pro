function unity_symmetric_sigmoidal, x,a,b,c
  d = 1.0d0
  return, d + (a - d)/(1.0d0 + (x/c)^b)
end
