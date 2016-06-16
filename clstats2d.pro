;-----------------------------------------------------------------------------
; CLSTATS2d.PRO
;-----------------------------------------------------------------------------
; Find clump statistics from data cube and clump
; assignment cube as created by program clfind
;
; The output of the program is very generic for maximum flexibility.
; Because of the different units found in FITS headers, the output
; positions and sizes are in pixels. You can convert these to header
; variables using, e.g. x=crval1+(i-crpxi1)*cdelt1.
; Similarly, only the clump integrated intensity is given.
;
; 6/10/04  jpw   modified from dark ages version
;
; needs: clstats2d.cb      ; common block
; usage: .run,clstats2d    ; to compile all procedures
;-----------------------------------------------------------------------------
pro clstats2d,file=file,log=log,silent=silent,fluxweight=fluxweight
@clstats2d.cb
@header.cb

if NOT keyword_set(fluxweight) then fluxweight=0
if NOT keyword_set(file) then begin
  print,'PROCEDURE clstats2d,file=filename,[/log]'
  print,'--------------------------------------------------'
  print,'filename = name of the FITS data cube (in quotes)'
  print,'           assumes .fits extension'
  print,'/log       for screen output copied to clstats2d.log'
  print,'--------------------------------------------------'
  return
endif

print,"-----------------------------------------------------------------------------------------------"
print,"CLSTATS2d: ",systime()
print,"-----------------------------------------------------------------------------------------------"
print,format='("Filename = ",a0)',file
if keyword_set(log) then begin
  openw,1,file+'_peaks.dat'
  printf,1,"#-----------------------------------------------------------------------------------------------"
  printf,1,"#CLSTATS2d: ",systime()
  printf,1,"#-----------------------------------------------------------------------------------------------"
  printf,1,format='("#Filename = ",a0)',file
endif

; read in data + header
data=readfits(file+'.fits',head)
readhd,head
nx=long(naxis1)
ny=long(naxis2)

; mask out bad data (NaN replaced by -999.9)
bad=where(finite(data) eq 0,count)
if (count gt 0) then data(bad)=-999.9

; read clump assign fle
assign=readfits(file+'.fits.clf',head)

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

if not(keyword_set(silent)) then begin
  print,"-----------------------------------------------------------------------------------------------"
  print," Ncl     x     y            MaxFlux     FWHMx     FWHMy       R            TotFlux    Npix"
  print,"-----------------------------------------------------------------------------------------------"
endif
if keyword_set(log) then begin
  printf,1,"#-----------------------------------------------------------------------------------------------"
  printf,1,"# Ncl    x     y            MaxFlux     FWHMx     FWHMy       R            TotFlux    Npix"
  printf,1,"#-----------------------------------------------------------------------------------------------"
endif

; MAIN LOOP
ncl_tot=max(assign)
for ncl=1,ncl_tot do begin
  dostats,ncl,log=log,silent=silent,fluxweight=fluxweight
endfor
if not(keyword_set(silent)) then $
  print,"-----------------------------------------------------------------------------------------------"
if keyword_set(log) then begin
  printf,1,"#-----------------------------------------------------------------------------------------------"
  close,1
endif

return
end
;-----------------------------------------------------------------------------
pro dostats,ncl,log=log,silent=silent,fluxweight=fluxweight
; note that the sizes are no longer corrected for the beam size
; (on account of possible mismatch in units).
; Best to do this yourself later...
@clstats2d.cb

clump_pix=where(assign eq ncl,npix)
j=clump_pix/nx
i=clump_pix-j*nx
flux=data(clump_pix)

; check to see if any pixel in the clump
; lie at the limits of the data cube
edge=0
ind=where(i eq nx-1 OR i eq 0,count)
if(count gt 0) then edge=1
ind=where(j eq ny-1 OR j eq 0,count)
if(count gt 0) then edge=edge+2


; peak position
if fluxweight then begin
    imin=min(i)
    imax=max(i)
    ilen=imax-imin+1
    jmin=min(j)
    jmax=max(j)
    jlen=jmax-jmin+1
    iarr=dindgen(ilen)+imin
    jarr=dindgen(jlen)+jmin
    iflux=dblarr(ilen)
    jflux=dblarr(jlen)
    for k=0,imax-imin do begin
        ipos=where(i eq iarr(k),ct)
        if ct gt 0 then iflux(k)=total(flux(ipos))
    endfor
    for k=0,jmax-jmin do begin
        jpos=where(j eq jarr(k),ct)
        if ct gt 0 then jflux(k)=total(flux(jpos))
    endfor
    i0=round(total(iarr*iflux)/total(iflux))
    j0=round(total(jarr*jflux)/total(jflux))
    pos=where(i eq i0 and j eq j0,ct)
    if ct gt 0 then peak=flux(pos) else peak=0.
endif else begin
    peak=max(flux)
    pos=where(data eq peak AND assign eq ncl)
    j0=pos(0)/nx
    i0=pos(0)-j0*nx
endelse

; calculate clump integrated intensity, size, and velocity dispersion
sumflux=total(flux)
ibar=total(i*flux)/sumflux
jbar=total(j*flux)/sumflux
sigi2=total(i*i*flux)/sumflux-ibar^2
sigj2=total(j*j*flux)/sumflux-jbar^2

; get clump size based on area on xy plane
radius=sqrt(npix/!pi)

; make dispersions into FWHM
sx=2.355*sqrt(sigi2)
sy=2.355*sqrt(sigj2)

; FLAG clumps if they lie on an edge
flag=""
if(edge eq 1) then flag="X"
if(edge eq 2) then flag="Y"
if(edge eq 3) then flag="XY"

if not(keyword_set(silent)) then $
  print,format='(i6,2(2x,i5),3x,f16.8,3(2x,f8.4),x,f16.4,2x,i6,x,a4)',ncl,i0,j0,peak,sx,sy,radius,sumflux,npix,flag
if keyword_set(log) then $
  printf,1,format='(i6,2(2x,i5),3x,f16.8,3(2x,f8.4),x,f16.4,2x,i6,x,a4)',ncl,i0,j0,peak,sx,sy,radius,sumflux,npix,flag

return
end
;-----------------------------------------------------------------------------
