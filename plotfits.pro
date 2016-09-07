;;;;;;;;;;;;;;;
;PLOT FITS MAP;
;;;;;;;;;;;;;;;

function logticks,axis,index,value
    exponent=long(alog10(value))
    tickmark='10!U'+strtrim(string(exponent),2)+'!N'
    return,tickmark
end

pro plotfits,file,loff,roff,ict,log=log
    
    if n_elements(loff) eq 0 then loff=0
    if n_elements(roff) eq 0 then roff=0
    if n_elements(ict) eq 0 then ict=0
    
    map = readfits(file+'.fits', hdr)

    fasthist, map

    loadct, ict
    reversect

    make_axes, hdr, raxis=raxis, daxis=daxis, rimg=rimg, dimg=dimg, vaxis=vaxis    
    
    nx=n_elements(map(*,0))
    ny=n_elements(map(0,*))
    npl=max([nx,ny])
    
    plmin=min(alog10(map),/nan)+loff
    plmax=max(alog10(map),/nan)-roff
    plrange=plmax-plmin
    
    if log then disp, alog10(map), min=plmin, max=plmax, /squarepix $
           else disp, map, min=10.^plmin, max=10.^plmax, /squarepix

    loadct,0
    
end



