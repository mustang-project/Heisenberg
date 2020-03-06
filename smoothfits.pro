function smoothfits,distance,incl,res_pc,filename,origdir,newdir,max_sample,astr_tolerance,tophat=tophat,regrid=regrid

if keyword_set(regrid) eq 0 then regrid=0
; motivate clfind levels

; TARGET RESOLUTION
  res_as=res_pc/distance/!dtor*3600.*sqrt(cos(incl))

if regrid then begin
;    NYQUIST RESOLUTION
     nyquist=res_as/max_sample

;    for resolutions of 0 - 5 round to 0.5 values
     ind = where(nyquist lt 5, ct)
     if ct ne 0 then nyquist = round(2*nyquist)/2.

;    for resolutions of 5 - 10 round to 1 values
     ind = where(nyquist ge 5 and nyquist lt 10, ct)
     if ct ne 0 then nyquist = round(nyquist)

;    for resolutions of 10 - 20 round to 2 values
     ind = where(nyquist ge 10 and nyquist lt 20, ct)
     if ct ne 0 then nyquist = 2*round(nyquist/2)

;    for resolutions of 20+ round to 5 values
     ind = where(nyquist ge 20, ct)
     if ct ne 0 then nyquist = 5*round(nyquist/5)
endif

;    LOAD ORIGINAL FILE
     file=origdir+filename
     data = readfits(file,hdr,/silent)

;    ORIGINAL BEAM SIZE
     bmaj = sxpar(hdr, 'BMAJ', count=ct)
     bmin = sxpar(hdr, 'BMIN', count=ct)
     if ct ne 1 then STOP
     beam = sqrt(bmaj*bmin)

;       CHECK IF DATA HAS COARSER RESOLUTION THAN TARGET RESOLUTION
        doconvolve=1
        if ~tophat && res_as le 1.05*beam*3600. then begin
            print,' Unable to process '+file+' at res '+strtrim(res_pc)+'pc'
            print,' Native res [arcsec] =',beam*3600.
            print,' Target res [arcsec] =',res_as
            print,' Native resolution is >95% of target resolution'
            print,' Output written at native resolution'
            doconvolve=0
        endif
        if tophat && res_as le 1.05*beam*3600. then print,' WARNING: tophat convolution at width < resolution FWHM -- proceeding ...'

        print,' Process '+file+' at res '+strtrim(res_pc)+'pc'

;       Construct Convolution Kernel
        sz = size(data)
        if doconvolve then begin
            cdelt = get_platescale(hdr, astr_tolerance)
            if ~tophat then begin
                psf_fwhm  = sqrt(res_as^2-(beam*3600.)^2) / 3600. / cdelt ; pixel units
                psf_width = floor(5.*psf_fwhm) < sxpar(hdr,'NAXIS1') < sxpar(hdr,'NAXIS2') ; NAXIS(psf)<=NAXIS(map)
                psf_width = 2.*ceil(psf_width/2.)-1. > 0. ;ensure odd number of pixels
                psf = psf_gaussian(npixel=psf_width, fwhm=[psf_fwhm,psf_fwhm], /normal) ;Gaussian kernel
            endif else psf = psf_tophat(res_as/3600./cdelt) ;tophat kernel

;           Check if map is 2D or 3D
            nans = where(finite(data, /nan), n_nans)
            if n_nans ne 0 then data(nans)=0.

            if sz[0] eq 2 then begin
               conv = convolve(data, psf) ;convolve in 2D
            endif else begin
               conv = fltarr(sz[1:3])
               for k = 0, sz[3]-1 do begin
                  counter, k+1, sz[3], 'Convolve channel '
                  conv[*,*,k] = convolve(data[*,*,k], psf, /no_ft) ;convolve in 3D
               endfor
            endelse

            data = conv ;replace data by convolved data

;           Set back Blanked Pixel
            if n_nans ne 0 then data[nans] = !values.f_nan

;           Update Header
            sxaddpar, hdr, 'BMAJ', res_as/3600. ;new beam major axis
            sxaddpar, hdr, 'BMIN', res_as/3600. ;new beam minor axis
        endif

if regrid then begin
;       NYQUIST SAMPLE MAP
        cdelt = 3600.*get_platescale(hdr, astr_tolerance)
        if nyquist ge 1.2*cdelt then begin
           naxis1 = sxpar(hdr, 'NAXIS1')
           naxis2 = sxpar(hdr, 'NAXIS2')

           cdelt1 = get_platescale(hdr, astr_tolerance, /xplatescale) ; get x-platescale with sign
           cdelt2 = get_platescale(hdr, astr_tolerance, /yplatescale) ; get y-platescale with sign
           crpix1 = sxpar(hdr, 'CRPIX1')
           crpix2 = sxpar(hdr, 'CRPIX2')
           crval1 = sxpar(hdr, 'CRVAL1')
           crval2 = sxpar(hdr, 'CRVAL2')

           scale = nyquist/cdelt
           hdr_new = hdr
           sxaddpar, hdr_new, 'NAXIS1', round(naxis1 / scale)
           sxaddpar, hdr_new, 'NAXIS2', round(naxis2 / scale)
           sxaddpar, hdr_new, 'CRPIX1', round(crpix1 / scale)
           sxaddpar, hdr_new, 'CRPIX2', round(crpix2 / scale)
           sxaddpar, hdr_new, 'CDELT1', cdelt1*scale
           sxaddpar, hdr_new, 'CDELT2', cdelt2*scale
           sxaddpar, hdr_new, 'CD1_1',  cdelt1*scale
           sxaddpar, hdr_new, 'CD2_2',  cdelt2*scale

;          Check if map is 2D or 3D
           if sz[0] eq 2 then begin
              hastrom, data, hdr, hdr_new $
                       , cubic=-0.5, interp=2, missing=!values.f_nan
           endif else begin
              for k = 0, sz[3]-1 do begin
                 counter, k+1, sz[3], 'Align channel '
                 chan = data[*,*,k]
                 hdr_2d = twod_head(hdr)
                 hastrom, chan, hdr_2d, hdr_new $
                          , cubic=-0.5, interp=2, missing=!values.f_nan
                 sz_chan = size(chan)
                 if k eq 0 then new = fltarr(sz_chan[1], sz_chan[2], sz[3])
                 new[*,*,k] = chan
              endfor
              sxaddpar, hdr_new, 'NAXIS',  3
              sxaddpar, hdr_new, 'NAXIS3', sxpar(hdr,'NAXIS3'), after='NAXIS2'
              hdr  = hdr_new
              data = new
           endelse

        endif
endif

;       SAVE RESULTS
        ;write fits
        dir=newdir
        dummy=file_search(dir,count=ct)
        if ct eq 0 then spawn,'mkdir '+dir
        writefits,newdir+filename,data,hdr

    return,data

end
