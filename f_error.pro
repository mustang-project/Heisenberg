;Standard error message

pro f_error,estring
    nstring=n_elements(estring)
    print, '; ---------- HEISENBERG ERROR MESSAGE ----------'
    for i=0,nstring-1 do print,'; ' + estring(i)
    print, '; ---------- HEISENBERG ERROR MESSAGE ----------'
    MESSAGE, 'Heisenberg Error encountered', /NoName
end
