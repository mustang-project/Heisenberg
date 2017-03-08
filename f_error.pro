;Standard error message

pro f_error,estring
    nstring=n_elements(estring)
    estring(0)='error: '+estring(0)
    for i=0,nstring-1 do print,' '+estring(i)
    stop,' quitting...'
end

