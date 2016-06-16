;;;;;;;;;;;;;;;;;;;;;;
;GENERATE TOP HAT PSF;
;;;;;;;;;;;;;;;;;;;;;;

function psf_tophat,width
    
    npix=2*round(width/2.)+1 ;dimension of psf array (always an odd number to have a pixel at the centre)
    frac=dblarr(npix,npix) ;fractional coverage array
    psf=dblarr(npix,npix) ;psf array
    index=findgen(npix) ;array with x and y pixel indices
    middle=0.5*(npix-1.) ;pixel value of top hat centre
    radius=0.5*width
    for i=0,npix-1 do begin
        for j=0,npix-1 do begin
            x1=index(i)-0.5-middle ;left x coord (with respect to middle)
            x2=index(i)+0.5-middle ;right x coord (with respect to middle)
            y1=index(j)-0.5-middle ;low y coord (with respect to middle)
            y2=index(j)+0.5-middle ;high y coord (with respect to middle)
            x=[x1,x1,x2,x2] ;x coordinates of pixel corners
            y=[y1,y2,y1,y2] ;y coordinates of pixel corners
            rad=sqrt(x^2.+y^2.) ;radii of pixel corners
            sortrad=sort(rad) ;index list of sorted radii
            x=x(sortrad) ;x sorted by distance
            y=y(sortrad) ;y sorted by distance
            rad=rad(sortrad) ;radii sorted by distance
            incircle=where(rad le radius,nin) ;pixel corners within circle
            if nin eq 0 then begin
                frac(i,j)=0. ;if entirely outside, then pixel value is zero
                if max(abs(x))/max(abs(y)) gt 1. then rmin=min(abs(x)) else rmin=min(abs(y)) ;minimum radius if part of central column or row
                if (i eq middle or j eq middle) and rmin le radius then $ ;if only part of a single edge overlaps with circle
                    frac(i,j)=radius^2.*acos(rmin/radius)-rmin*sqrt(2.*radius*(radius-rmin)-(radius-rmin)^2.)
                if i eq middle and j eq middle then $
                    frac(i,j)=!pi*radius^2.-4.*(radius^2.*acos(0.5/radius)-0.5*sqrt(2.*radius*(radius-0.5)-(radius-0.5)^2.)) ;if central, then fraction of central pixel
            endif
            if nin eq 4 then frac(i,j)=1. ;if entirely within, then fractional coverage is unity
            if nin ne 0 and nin ne 4 then begin
                if nin eq 1 then begin ;if one is within circle
                    lx=sqrt(radius^2.-y(0)^2.)-abs(x(0)) ;triangle x leg
                    ly=sqrt(radius^2.-x(0)^2.)-abs(y(0)) ;triangle y leg
                    a=sqrt(lx^2.+ly^2.) ;circular segment base length
                    area=0.5*lx*ly ;triangle area
                endif
                if nin eq 2 then begin ;if two are within circle
                    if x(0) eq x(1) then begin
                        ccst=x ;coordinate with same index within circle 
                        cvar=y ;coordinate with different index within circle
                    endif else begin
                        ccst=y ;coordinate with same index within circle 
                        cvar=x ;coordinate with different index within circle
                    endelse
                    l1=sqrt(radius^2.-cvar(0)^2.)-abs(ccst(0)) ;long edge of trapezoid
                    l2=sqrt(radius^2.-cvar(1)^2.)-abs(ccst(1)) ;short edge of trapezoid
                    dl=l1-l2 ;edge length difference
                    a=sqrt(1.+dl^2.) ;circular segment base length
                    rmax=max(abs([x,y])) ;minimum radius if part of central column or row
                    if (i eq middle or j eq middle) and rmax le radius then $
                        minarea=radius^2.*acos(rmax/radius)-rmax*sqrt(2.*radius*(radius-rmax)-(radius-rmax)^2.) $ ;part of circlular segment outside pixel
                        else minarea=0.
                    area=l2+0.5*dl-minarea ;trapezoid area (width is unity)
                endif
                if nin eq 3 then begin ;if three are within circle
                    lx=abs(x(3))-sqrt(radius^2.-y(3)^2.)
                    ly=abs(y(3))-sqrt(radius^2.-x(3)^2.)
                    a=sqrt(lx^2.+ly^2.) ;circular segment base length
                    area=1.-0.5*lx*ly ;square minus trangle area
                endif
                rcs=0.5*sqrt(4.*radius^2.-a^2.) ;distance to middle of circular segment base
                hcs=radius-rcs ;height of circular segment base
                areacs=radius^2.*acos(rcs/radius)-rcs*sqrt(2.*radius*hcs-hcs^2.) ;area of circular segment
                frac(i,j)=(area+areacs) ;fractional coverage equals area covered within pixel by top hat
            endif
            psf(i,j)=frac(i,j)*(!pi*radius^2.)^(-1.) ;divide fractional coverage by area to get normalised psf
;stop
        endfor
    endfor

    return,psf
    
end





