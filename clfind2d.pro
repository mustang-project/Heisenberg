;-----------------------------------------------------------------------------
; CLFIND2d.PRO
;-----------------------------------------------------------------------------
; Find clumps in a x-y fits dataset
; based on the algorithm described in
; Williams, de Geus, & Blitz 1994, ApJ, 428, 693
; (please cite this paper in publications that use this procedure)
;
; 6/10/04   jpw
;
; needs: clfind2d.cb      ; common block
; usage: .run,clfind2d    ; to compile all procedures
;-----------------------------------------------------------------------------
pro clfind2d,file=file,levels=levels,log=log,npixmin=npixmin
; read in fits cube and call clump finding routines
@clfind2d.cb
@header.cb

if NOT keyword_set(file) then begin
  print,'PROCEDURE clfind2d,file=filename,levels=[l1,l2,l3...],[/log]'
  print,'------------------------------------------------------------'
  print,'filename = root name of the FITS data cube (in quotes)'
  print,'           assumes .fits extension'
  print,'levels   = vector of contour levels'
  print,'/log       for screen output copied to clfind2d.log'
  print,'------------------------------------------------------------'
  return
endif 

print,"----------------------------------------------------------------"
print,"CLFIND2d: ",systime()
print,"----------------------------------------------------------------"
print,format='("Filename = ",a0)',file
print,format='("Contour levels at =",f10.6)',levels
print,"----------------------------------------------------------------"
if keyword_set(log) then begin
  openw,1,file+'_clfind.dat'
  printf,1,"#----------------------------------------------------------------"
  printf,1,"#CLFIND2d: ",systime()
  printf,1,"#----------------------------------------------------------------"
  printf,1,format='("#Filename = ",a0)',file
  printf,1,format='("#Contour levels at =",f10.6)',levels
  printf,1,"#----------------------------------------------------------------"
endif

; read in data + header
data=readfits(file+'.fits',header,/silent)
readhd,header
nx=long(naxis1)
ny=long(naxis2)

; save for common block
infile=file
levs=levels

; mask out bad data (NaN replaced by -999.9)
bad=where(finite(data) eq 0,count)
if (count gt 0) then data(bad)=-999.9

; initialize clump assignment arrays
max_clumps=99999
assign=intarr(nx,ny)
clump_peak=intarr(max_clumps,2)


; MAIN CLUMP FINDING LOOP
ncl=0
t0=systime(1)
nlev=n_elements(levs)
levs=[levs,99999]

for nwork=nlev-1L,0,-1 do begin
  defreg,nwork,npix,reg,nreg
  
  print,format='("Contour level ",I3,": ",I7," pixels ",I7," regions ",$)'$
       ,levs(nwork),npix,nreg
  if keyword_set(log) then $
    printf,1,format='("#Contour level ",I3,": ",I7," pixels ",I7," regions ",$)'$
            ,levs(nwork),npix,nreg
  print,''
  defclump,nwork,reg,nreg,nnew
  print,format='(i7," new clumps")',nnew
  if keyword_set(log) then printf,1,format='(i7,"# new clumps")',nnew
  ncl=ncl+nnew

endfor

; reject clumps with =< npixmin pixels
if n_elements(npixmin) eq 0 then npixmin=20
testbad,npixmin,nbad
ncl=ncl-nbad
print,format='(i7," clumps found (",i7," rejected)")',ncl,nbad
print,"================================================================="
if keyword_set(log) then begin
  printf,1,format='("#",i7," clumps found (",i7," rejected)")',ncl,nbad
  printf,1,"#================================================================="
endif

; write assign file to a fits file
outfile=infile+".fits.clf"
mk_hdr,header,outheader
writefits,outfile,assign,outheader
print,"Writing output file: ",outfile
print,""
print,"CLFIND2d exits sucessfully"

delta_t=(systime(1)-t0)/60.0
print,format='(f5.1," minutes elapsed")',delta_t

if keyword_set(log) then begin
  printf,1,format='("#Writing output file: ",a0)',outfile
  printf,1,format='("#",f5.1," minutes elapsed")',delta_t
  close,1
endif

return
end
;-----------------------------------------------------------------------------
pro defreg,nwork,npix0,reg,nreg
; define regions at contour level nwork
@clfind2d.cb

; get all points at this contour level
levmin=levs(nwork)
levmax=levs(nwork+1)
levpix=where(data ge levmin AND data lt levmax,npix0)
if(npix0 eq 0) then return

; create region array for working level
reg=intarr(nx,ny)
nreg=0

; set all pixels not in level to -1
reg(levpix)=1
reg=reg-1

; loop through finding sucessive unassigned peaks -> new regions
pix1=where(reg eq 0,npix1)
npix1init=npix1
while (npix1 gt 0) do begin
  nreg=nreg+1
  dpeak=max(data(pix1),peak)
  jpeak=pix1(peak)/nx
  ipeak=pix1(peak)-jpeak*nx
  pix=search2d(data,ipeak,jpeak,levmin,levmax,/diagonal)
  reg(pix)=nreg
  pix1=where(reg eq 0,npix1)
  progress,' defreg',npix1init-npix1,npix1init
endwhile

return
end
;-----------------------------------------------------------------------------
pro defclump,nwork,reg,nreg,nnew
; define clumps at working contour level
@clfind2d.cb

levmin=levs(nwork)
nnew=0
if (nreg eq 0) then return

; extend each defined region upwards and see
; if it merges with previously defined clumps
for nr=1L,nreg do begin
  pix1=where(reg eq nr,npix1)
  dpeak=max(data(pix1),peak)
  jpeak=pix1(peak)/nx
  ipeak=pix1(peak)-jpeak*nx
  pix2=search2d(data,ipeak,jpeak,levmin,99999,/diagonal)
  apix2=assign(pix2)

  if(max(apix2) eq 0) then begin
;    print,"Found new clump!"
    nnew=nnew+1
    assign(pix2)=ncl+nnew
    dpeak=max(data(pix2),peak)
    jpeak=pix2(peak)/nx
    ipeak=pix2(peak)-jpeak*nx
    clump_peak(ncl+nnew,0)=ipeak
    clump_peak(ncl+nnew,1)=jpeak

  endif else begin
;   define list of merged clumps
;   (note that includes unassigned pixels at the working level)
    nc=apix2(uniq(apix2,sort(apix2)))

    if(min(nc) eq 0) then begin
      if(n_elements(nc) eq 2) then begin
;     if just one clump above working level, then extend it
;        print,"Extending clump ",nc(1)
        assign(pix2)=nc(1)

      endif else begin
;       if more than one clump, this region is a merger
        clump_merge=nc(1:*)
;        print,"Merging clumps ",clump_merge
        ic=clump_peak(clump_merge,0)
        jc=clump_peak(clump_merge,1)

;       go through pixel by pixel and assign to nearest clump
        for nr1=0L,npix1-1 do begin
          j=pix1(nr1)/nx
          i=pix1(nr1)-j*nx
          d=(i-ic)^2 + (j-jc)^2
          dmin=min(d,m)
          assign(pix1(nr1))=clump_merge(m)
        endfor
      endelse
    endif

  endelse
  progress,' defclump',nr-1,nreg-1
endfor


return
end
;-----------------------------------------------------------------------------
pro testbad,nmin,nbad
; sorts clumps in order of peak flux and
; rejects those with number of pixels <= nmin
@clfind2d.cb

nc=indgen(ncl)+1
dmax=data(clump_peak(nc,0),clump_peak(nc,1))
new_order=reverse(sort(dmax))

ncl_new=0
nbad=0
assign0=assign

for n1=1L,ncl do begin
  n0=new_order(n1-1)+1
  iclp=where(assign0 eq n0,count)

  if (count lt nmin) then begin
    nbad=nbad+1
    assign(iclp)=0
  endif else begin
    ncl_new=ncl_new+1
    assign(iclp)=ncl_new
  endelse

endfor


return
end
;-----------------------------------------------------------------------------
pro mk_hdr,header,header_out
;------------------------------------
; create fits header for output file
; (basically a copy of input header
; with some program information)
;------------------------------------
@clfind2d.cb

s=size(header)
nlines=s(1)
header_out=strarr(99)
header_out(0)="SIMPLE  =                    T /"
header_out(1)="BITPIX  =                   16 /"
header_out(2)="NAXIS   =                    2 /"
n=2
for i=0L,nlines-1 do begin
  if (strpos(header(i),'NAXIS1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'NAXIS2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRVAL1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRVAL2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRPIX1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRPIX2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CDELT1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CDELT2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CTYPE1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CTYPE2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
end

nl=n_elements(levs)-1
for m=0L,nl-1 do header_out(n+m+1)=$
  "HIERARCH CLFIND2d contour = "+string(levs(m),format='(f12.8)')+"    /"


; Add final END statement
for i=0L,nlines-1 do begin
  if (strpos(header(i),'END') eq 0) then header_out(n+nl+1)=header(i)
endfor

; get rid of blank lines at end
header_out=header_out(0:n+nl+1)

return
end
;-----------------------------------------------------------------------------
