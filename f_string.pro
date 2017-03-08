;String conversion function

function f_string,numorig,ndec ;converts any real number to string with up to four decimals
    if ndec gt 4 then stop ;only works with floats and up to four decimals
    arraylen=n_elements(numorig)
    num=fltarr(arraylen)
    numstr=strarr(arraylen)
    for i=0L,arraylen-1 do begin
        num(i)=float(10.^(-ndec)*round(10.^ndec*numorig(i)))
        sign=num(i) ge 0.
        lsign=abs(num(i)) ge 1.
        nlog=fix(alog10(abs(num(i))+1.d-5))+lsign-1
        n1=6+(num(i) ne 0.)*((nlog lt 0)*nlog+sign-1)
        n2=ndec+2-(ndec eq 0)+(num(i) ne 0.)*((nlog gt 0)*nlog-sign+1)
        numstr(i)=strmid(num(i),n1,n2)
    endfor
    if arraylen eq 1 then return,strmid(num(0),n1(0),n2(0)) else return,numstr
end
