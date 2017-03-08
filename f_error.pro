;Standard error message

pro f_error,estring
    nstring=n_elements(estring)
    print, '; ---------- KL14 ERROR MESSAGE ----------'
    for i=0,nstring-1 do print,'; ' + estring(i)
    print, '; ---------- KL14 ERROR MESSAGE ----------'
    MESSAGE, 'KL14 Error encountered', /NoName
end

