;-----------------------------------------------------------------------------
; CLPLOT2d.PRO
;-----------------------------------------------------------------------------
; plot the clumps found by clfind2d
;
; 6/10/04  jpw
;
; needs: clplot2d.cb      ; common block
; usage: .run,clplot2d    ; to compile all procedures
;-----------------------------------------------------------------------------
pro clplot2d,file=file
@clplot2d.cb
@header.cb

if NOT keyword_set(file) then begin
  print,'PROCEDURE clplot2d,file=filename'
  print,'-------------------------------------------------------'
  print,'filename = root name of the FITS data cube (in quotes)'
  print,'           assumes .fits and .fits.clf extensions'
  print,'-------------------------------------------------------'
  return
endif

print,'-------------------------------------------------------'
print,"CLPLOT2d: ",systime()
dfile=file+".fits"
afile=file+".fits.clf"
print,format='("Reading ",a0)',dfile
data=readfits(dfile,head)
readhd,head
nx=long(naxis1)
ny=long(naxis2)
x0=crval1
y0=crval2
if(n_elements(cdelt1) gt 0) then dx=cdelt1
if(n_elements(cd1_1) gt 0) then dx=cd1_1
if(n_elements(cdelt2) gt 0) then dy=cdelt2
if(n_elements(cd2_2) gt 0) then dy=cd2_2
i0=crpix1
j0=crpix2

print,format='("Reading ",a0)',afile
assign=readfits(afile,head)

; get clfind parameters
s=size(head)
nlines=s(1)
levels=fltarr(99)
nlevs=0
for i=0,nlines-1 do begin
  if (strpos(head(i),'contour') gt 0) then begin
    s=strmid(head(i),28,10)
    nlevs=nlevs+1
    levels(nlevs)=float(strtrim(s))
  endif
endfor
levels=levels(1:nlevs)
print,format='("Contour levels at =",i3)',levels
nlevs=n_elements(levels)
print,'-------------------------------------------------------'

; MAIN LOOP
ncl=1
while (ncl gt 0) do begin
  read,str_ncl,prompt="Enter clump number (RETURN or 0 to end): ",format='(a0)'
  if(strlen(str_ncl) gt 0) then begin
    reads,str_ncl,ncl
    makeplot,ncl,levels
  endif else begin
    ncl=0
  endelse
end

print,'-------------------------------------------------------'
print,"CLPLOT2d exits sucessfully"


return
end
;---------------------------------------------------------------------------
pro makeplot,ncl,levs
@clplot2d.cb
if (ncl le 0) then return

; map coordinates (offsets) and limits
x=(findgen(nx)-i0)*dx
y=(findgen(ny)-j0)*dy
xmin=min(x)
xmax=max(x)
ymin=min(y)
ymax=max(y)

; convert to arcmin
x=60*x
y=60*y
xmin=60*xmin
xmax=60*xmax
ymin=60*ymin
ymax=60*ymax

; aspect ratio of map and display
asp=(ymax-ymin)/(xmax-xmin)
if(asp gt 1) then begin
  xs=600.0
  ys=776.0
endif else begin
  xs=776.0
  ys=600.0
endelse
asp0=ys/xs

; set up screen device
device,decomposed=0
window,0,xsize=xs,ysize=ys
loadct,39

; set up plot position on the screen
; so that the angular units are the same on each axis
sx=0.75
sy=sx*asp/asp0
if(sy gt 0.88) then begin
  sy=0.88
  sx=sy*asp0/asp
endif
px1=0.5*(1-sx)+0.04
px2=0.5*(1+sx)+0.04
py1=0.5*(1-sy)+0.04
py2=0.5*(1+sy)+0.04
pos=[px1,py1,px2,py2]


bmap=bytscl(data,min=0.5*min(levs),max=max(data),/NaN,top=254)
tvimage,bmap,position=pos,/erase,/minus_one

tic=0.035
contour,data,x,y $
       ,position=pos $
       ,levels=levs,c_thick=1,c_colors=255 $
       ,xra=[xmax,xmin],xsty=1 $
       ,yra=[ymin,ymax],ysty=1 $
       ,xtitle="!7Da !6(arcmin)" $
       ,ytitle="!7Dd !6(arcmin)" $
       ,xticklen=tic/asp,yticklen=tic $
       ,xminor=4,yminor=4 $
       ,charsize=1.4,charthick=2 $
       ,/noerase

; plot clump
clump=fltarr(nx,ny)
pix=where(assign eq ncl,npix)
if(npix eq 0) then begin
  print,'No clump!'
  print,'Max clump number =',max(assign)
  return
endif
clump(pix)=data(pix)
contour,clump,x,y $
       ,position=pos $
       ,levels=levs,c_thick=3,c_colors=255 $
       ,xra=[xmax,xmin],xsty=5 $
       ,yra=[ymin,ymax],ysty=5 $
       ,/noerase,color=fsc_color('red')

return
end
;---------------------------------------------------------------------------
