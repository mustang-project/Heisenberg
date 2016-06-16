;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CUT MAP AREA TO BE COVERED BY ALL MAPS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function cutarea,map1,map2,map3,edge
    
    if n_elements(edge) eq 0 then edge=0
    map1=double(map1)
    map2=double(map2)
    map3=double(map3)
    if mean([size(map1) eq size(map2),size(map1) eq size(map3)]) ne 1. then begin
        print,' CUTAREA: maps have different pixel dimensions, quitting ...'
        stop
    endif
    
    sz=size(map1)
    nx=sz(1)
    ny=sz(2)
    x1=dblarr(nx)
    y1=dblarr(ny)
    x2=dblarr(nx)
    y2=dblarr(ny)
    x3=dblarr(nx)
    y3=dblarr(ny)
    for i=0,nx-1 do begin
        x1(i)=total(map1(i,*),/nan)
        x2(i)=total(map2(i,*),/nan)
        x3(i)=total(map3(i,*),/nan)
    endfor
    for i=0,ny-1 do begin
        y1(i)=total(map1(*,i),/nan)
        y2(i)=total(map2(*,i),/nan)
        y3(i)=total(map3(*,i),/nan)
    endfor
    
    xmin1=min(where(x1 ne 0.))
    xmax1=max(where(x1 ne 0.))
    ymin1=min(where(y1 ne 0.))
    ymax1=max(where(y1 ne 0.))
    xmin2=min(where(x2 ne 0.))
    xmax2=max(where(x2 ne 0.))
    ymin2=min(where(y2 ne 0.))
    ymax2=max(where(y2 ne 0.))
    xmin3=min(where(x3 ne 0.))
    xmax3=max(where(x3 ne 0.))
    ymin3=min(where(y3 ne 0.))
    ymax3=max(where(y3 ne 0.))
    
    xmin=round(min([xmin1,xmin2,xmin3])+edge)
    xmax=round(max([xmax1,xmax2,xmax3])-edge)
    ymin=round(min([ymin1,ymin2,ymin3])+edge)
    ymax=round(max([ymax1,ymax2,ymax3])-edge)
    
    incl=dblarr(nx,ny)
    incl(xmin:xmax,ymin:ymax)=1 ;set indices to be included to unity

    return,incl
    
end





