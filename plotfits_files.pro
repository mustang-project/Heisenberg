;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOT FITS MAP FOR MORE THAN ONE FILE;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function logticks,axis,index,value
    exponent=long(alog10(value))
    tickmark='10!U'+strtrim(string(exponent),2)+'!N'
    return,tickmark
end

pro plotfits_files,files,galaxy,name,loffs,roffs,ict,log=log
    
    if n_elements(loffs) eq 0 then loffs=0
    if n_elements(roffs) eq 0 then roffs=0
    if n_elements(ict) eq 0 then ict=0
    
set_plot,'ps'
device,filename='galaxies/'+galaxy+'/figures/'+name+'maps.ps'

nfiles=n_elements(files)
!p.multi=[0,nfiles,1,0,0]
for i=0,nfiles-1 do begin
    multiplot
    map = readfits(files(i)+'.fits', hdr)

    fasthist, map

    loadct, ict
    reversect

    make_axes, hdr, raxis=raxis, daxis=daxis, rimg=rimg, dimg=dimg, vaxis=vaxis    
    
    nx=n_elements(map(*,0))
    ny=n_elements(map(0,*))
    npl=max([nx,ny])
    
    plmin=min(alog10(map),/nan)+loffs(i)
    plmax=max(alog10(map),/nan)-roffs(i)
    plrange=plmax-plmin
        
    if log then disp, alog10(map), min=plmin, max=plmax, /squarepix,/noerase, xstyle=4, ystyle=4 $
           else disp, map, min=10.^plmin, max=10.^plmax, /squarepix,/noerase, xstyle=4, ystyle=4
endfor
multiplot,/reset
device,/close
set_plot,'x'
!p.multi=0
    loadct,0
    
end



