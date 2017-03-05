;;;;;;;;;;;;;;;;;;;;;;
;IDENTIFY PEAKS IN 2D;
;;;;;;;;;;;;;;;;;;;;;;

function peaks2d,map,levels=levels,log=log,fluxweight=fluxweight,npixmin=npixmin
    
    if n_elements(log) eq 0 then log=0
    if n_elements(fluxweight) eq 0 then fluxweight=0
    
    resolve_routine,'clfind2d'
    if log then clfind2d,file=map,levels=levels,/log,npixmin=npixmin else clfind2d,file=map,levels=levels,npixmin=npixmin
    
    clstats2d,file=map,log=log,/silent,fluxweight=fluxweight
    file=map+'_peaks.dat'
    nlines=min(file_lines(file)-8)
    if nlines le 2 then f_error,'insufficient number of peaks identified ('+f_string(nlines,0)+')' ;analysis fundamentally does not work if number of peaks is 2 or less
    n=dblarr(nlines)
    x=dblarr(nlines)
    y=dblarr(nlines)
    maxflux=dblarr(nlines)
    fwhmx=dblarr(nlines)
    fwhmy=dblarr(nlines)
    r=dblarr(nlines)
    totflux=dblarr(nlines)
    npix=dblarr(nlines)
    nval=0
    openr,lun,file,/get_lun ;open the file
    skip_lun,lun,7,/lines ;skip the 7-line header
    line='' ;initialise line
    for i=0,nlines-1 do begin ;read lines until end string is reached
        readf,lun,n0,x0,y0,maxflux0,fwhmx0,fwhmy0,r0,totflux0,npix0 ;read line containing variable
        n(i)=n0
        x(i)=x0
        y(i)=y0
        maxflux(i)=maxflux0
        fwhmx(i)=fwhmx0
        fwhmy(i)=fwhmy0
        r(i)=r0
        totflux(i)=totflux0
        npix(i)=npix0
    endfor
    close,lun ;close the logical unit
    free_lun,lun ;make the logical unit available again
        
    coordinates=[[x],[y],[totflux],[npix]]
    
    return,coordinates
    
end



