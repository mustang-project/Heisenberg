;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DERIVE PHYSICAL QUANTITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function f_getdvar,var ;Note: input must be positive and real
    COMMON numbers
    nvar=n_elements(var)
    dvar=dblarr(nvar)
    for i=0L,nvar-1 do begin
        if i eq 0 then dvar(i)=min([2.*var(i),var(i+1)-var(i)]) ;make sure that var-0.5*dvar cannot be negative
        if i eq nvar-1 then dvar(i)=var(i)-var(i-1)
        if i gt 0 and i lt nvar-1 then dvar(i)=.5*(var(i+1)-var(i-1))
        if dvar(i) eq 0. then dvar(i)=tiny
    endfor
    return,dvar
end

function f_createpdf,array,darray,cumpdf,nnew
    COMMON numbers
    norig=n_elements(array) ;originial number of elements
    minval=min(array)-.5*darray(0) ;assume linearly-symmetric bins as input
    maxval=max(array)+.5*darray(norig-1) ;assume linearly-symmetric bins as input
    pdf=dblarr(nnew)
    if minval le 0. then logrange=0 else logrange=1
    if minval ne maxval then begin
        if logrange then begin
            newarray=minval*(maxval/minval)^((dindgen(nnew)+.5)/nnew)
            left=newarray*(maxval/minval)^(-1./(2.*nnew)) ;left bin edge
            right=newarray*(maxval/minval)^(1./(2.*nnew)) ;right bin edge
            newdarray=right-left
            if max(newdarray) eq 0 then logrange=0
        endif
        if ~logrange then begin
            newarray=minval+(maxval-minval)*((dindgen(nnew)+.5)/nnew)
            left=newarray-(maxval-minval)/(2.*nnew) ;left bin edge
            right=newarray+(maxval-minval)/(2.*nnew) ;right bin edge
            newdarray=0.*left+(maxval-minval)/nnew
        endif
        minlimit=min(where(array gt min(array),nmin))
        maxlimit=max(where(array lt max(array),nmax))
        if nmin eq 0 then interpolmin=0 else interpolmin=minlimit ;avoid interpolation over saturated maximum values (crashes interpol function) and just include number once
        if nmax eq 0 then interpolmax=n_elements(cumpdf)-1 else interpolmax=maxlimit ;avoid interpolation over saturated maximum values (crashes interpol function) and just include number once
        interpolarr=[0,findgen(interpolmax-interpolmin+1)+interpolmin] ;index array for interpolation
        if interpolmax lt n_elements(cumpdf)-1 then interpolarr=[interpolarr,n_elements(cumpdf)-1] ;update index array for interpolation
        for i=0L,nnew-1 do begin
            if left(i) le min(array) then cumpdfleft=max([0.,interpol(cumpdf(interpolarr),array(interpolarr),left(i))]) ;left cumpdf value
            if left(i) gt max(array) then cumpdfleft=1. ;by definition
            if left(i) gt min(array) && left(i) le max(array) then cumpdfleft=interpol(cumpdf,array,left(i)) ;left cumpdf value
            if left(i) eq max(array) then cumpdfleft=min([interpol(cumpdf(interpolarr),array(interpolarr),left(i)),1.]) ;left cumpdf value
            if right(i) lt min(array) then cumpdfright=0. ;by definition
            if right(i) ge max(array) then cumpdfright=min([interpol(cumpdf(interpolarr),array(interpolarr),right(i)),1.]) ;right cumpdf value
            if right(i) lt max(array) && right(i) ge min(array) then cumpdfright=interpol(cumpdf,array,right(i)) ;right cumpdf value
            pdf(i)=(cumpdfright-cumpdfleft)/newdarray(i) ;pdf value in bin
        endfor
    endif else begin
        newarray=dblarr(nnew)
        newdarray=dblarr(nnew)+tiny
        for i=0L,nnew-1 do newarray(i)=array(0)+(i-nnew)/2*tiny
        pdf(*)=(nnew*tiny)^(-1.)
    endelse
    if finite(total(pdf)) eq 0 then stop
    return,[[newarray],[newdarray],[pdf]]
end

function f_writecorr,matrix,complete,galaxy,outputdir
    openw,lun,outputdir+galaxy+'_correlation_matrix.dat',/get_lun
    printf,lun,'# 2D correlation coefficient matrix between all constrained quantities for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how this matrix was calculated'
    printf,lun,'# This symmetric array lists the correlation coefficients for quantities that are, from left to right AND from top to bottom, in the order:'
    if complete then begin
        printf,lun,'# tgas, tover, lambda, tstar, ttotal, betastar, betagas, surfglobalstar, surfglobalgas, surfconstar, surfcongas, rpeakstar, rpeakgas, zetastar, zetagas, vfb, surfsfr, surfgas, tdepl, esf, mdotsf, mdotfb, etainst, etaavg, chie, chip, fcl, fgmc'
    endif else begin
        printf,lun,'# tgas, tover, lambda, tstar, ttotal, betastar, betagas, surfglobalstar, surfglobalgas, surfconstar, surfcongas, rpeakstar, rpeakgas, zetastar, zetagas, vfb, fcl, fgmc'
    endelse
    printf,lun,''
    n=n_elements(matrix(0,*))
    nstr=f_string(n-1,0)
    for i=0L,n-1 do begin
        printf,lun,format='(f12.5,'+nstr+'(3x,f12.5))',reform(matrix(i,*))
    endfor
    close,lun
    free_lun,lun
    report='         2D correlation coefficient matrix written'
    return,report
end

function f_writecov,matrix,complete,galaxy,outputdir
    openw,lun,outputdir+galaxy+'_covariance_matrix.dat',/get_lun
    printf,lun,'# 2D covariance matrix between all constrained quantities for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how this matrix was calculated'
    printf,lun,'# This symmetric array lists the covariances for quantities that are, from left to right AND from top to bottom, in the order:'
    if complete then begin
        printf,lun,'# tgas, tover, lambda, tstar, ttotal, betastar, betagas, surfglobalstar, surfglobalgas, surfconstar, surfcongas, rpeakstar, rpeakgas, zetastar, zetagas, vfb, surfsfr, surfgas, tdepl, esf, mdotsf, mdotfb, etainst, etaavg, chie, chip, fcl, fgmc'
    endif else begin
        printf,lun,'# tgas, tover, lambda, tstar, ttotal, betastar, betagas, surfglobalstar, surfglobalgas, surfconstar, surfcongas, rpeakstar, rpeakgas, zetastar, zetagas, vfb, fcl, fgmc'
    endelse
    printf,lun,''
    n=n_elements(matrix(0,*))
    nstr=f_string(n-1,0)
    for i=0L,n-1 do begin
        printf,lun,format='(f12.5,'+nstr+'(3x,f12.5))',reform(matrix(i,*))
    endfor
    close,lun
    free_lun,lun
    report='         2D covariance matrix written'
    return,report
end

function derivephys,surfsfr,surfsfr_err,surfgas,surfgas_err,area,tgas,tover,lambda,beta_star,beta_gas,fstarover,fgasover,fcl,fgmc,tstariso,tstariso_rerrmin,tstariso_rerrmax,tstar_incl, $
                    surfglobals,surfglobalg,surfcontrasts,surfcontrastg,apertures_star,apertures_gas,psie,psip,peak_prof,ntry,nphysmc,galaxy,outputdir,arrdir,figdir,map_units, $
                    emfrac_cor_mode, filter_choice, filter_len_conv, $
                    rpeak_cor_mode, rpeaks_cor_val, rpeaks_cor_emin, rpeaks_cor_emax, rpeakg_cor_val, rpeakg_cor_emin, rpeakg_cor_emax

    if map_units eq 1 then complete=1 else complete=0

    ;derived quantities for best-fitting model
    if tstar_incl eq 0 then tstar=tstariso+tover else tstar=tstariso
    ttotal=tgas+tstar-tover
    betastar=interpol(beta_star,fstarover,tover/tstar)
    betagas=interpol(beta_gas,fgasover,tover/tgas)
    surfglobals_fit=10.^interpol(alog10(surfglobals),alog10(apertures_star),alog10(lambda))
    surfglobalg_fit=10.^interpol(alog10(surfglobalg),alog10(apertures_gas),alog10(lambda))
    surfs_fit=10.^interpol(alog10(surfcontrasts),alog10(apertures_star),alog10(lambda))
    surfg_fit=10.^interpol(alog10(surfcontrastg),alog10(apertures_gas),alog10(lambda))
    zetastar=f_zeta(surfs_fit,ttotal,tstar,peak_prof)
    zetagas=f_zeta(surfg_fit,ttotal,tgas,peak_prof)
    rpeakstar=f_rpeak(zetastar,lambda)
    rpeakgas=f_rpeak(zetagas,lambda)
    vfb=f_vfb(tover,lambda)
    if complete then begin
        tdepl=surfgas/surfsfr/1.d9
        esf=f_esf(tgas,tdepl,fcl,fgmc)
        mdotsf=f_mdotsf(tgas,lambda,surfgas*surfglobalg_fit,fgmc,esf)
        mdotfb=f_mdotfb(tover,lambda,surfgas*surfglobalg_fit,fgmc,esf)
        etainst=f_etainst(tgas,tover,esf)
        etaavg=f_etaavg(esf)
        chie=f_chie(tover,esf,vfb,psie)
        chip=f_chip(tover,esf,vfb,psip)
    endif
    if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
      cut_len_temp = filter_len_conv * lambda

      if rpeak_cor_mode eq 0 then begin ; use measured rpeak value
        fwhm_star_temp = rpeakstar * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
        fwhm_gas_temp = rpeakgas * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
      endif else if rpeak_cor_mode eq 1 then begin ; use supplied rpeak value
        fwhm_star_temp = rpeaks_cor_val * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
        fwhm_gas_temp = rpeakg_cor_val * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
      endif

      qconstar = fourier_flux_loss_correction(filter_choice, cut_len_temp/fwhm_star_temp)
      qcongas = fourier_flux_loss_correction(filter_choice, cut_len_temp/fwhm_gas_temp)
      fcl /= qconstar ; flux loss correction
      fgmc /= qcongas  ; flux loss correction
    endif else begin
      qconstar = 1.0d0
      qcongas = 1.0d0
    endelse
    if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin

      if rpeak_cor_mode eq 0 then begin ; use measured rpeak value
        etastar = sqrt(tstar/ttotal) * zetastar ; zeta corrected for timeline. i.e. peak density in image
        etagas = sqrt(tgas/ttotal) * zetagas ; zeta corrected for timeline. i.e. peak density in image
      endif else if rpeak_cor_mode eq 1 then begin ; use supplied rpeak value
        ; old method
        ; ; etastar = (sqrt(tstar/ttotal) * zetastar) / (rpeakstar/rpeaks_cor_val) ; eta with correction for flux loss from signal regions
        ; ; etagas = (sqrt(tgas/ttotal) * zetagas) / (rpeakgas/rpeakg_cor_val) ; eta with correction for flux loss from signal regions
        ; new method
        etastar = (sqrt(tstar/ttotal) * 2.0* (rpeaks_cor_val/lambda))
        etagas = (sqrt(tgas/ttotal) * 2.0* (rpeakg_cor_val/lambda))



      endif

      qzetastar = fourier_zeta_correction(filter_choice, etastar)
      qzetagas = fourier_zeta_correction(filter_choice, etagas)

      fcl  /=  qzetastar ; zeta correction
      fgmc  /= qzetagas ; zeta correction
    endif else begin
      qzetatar = 1.0d0
      qzetagas = 1.0d0
    endelse

    restore,filename=arrdir+'probnorm.sav'
    restore,filename=arrdir+'probtgastover.sav'
    restore,filename=arrdir+'probtgas.sav'
    restore,filename=arrdir+'tgasarr.sav'
    restore,filename=arrdir+'toverarr.sav'
    restore,filename=arrdir+'lambdaarr.sav'
    restore,filename=arrdir+'dtgas.sav'
    restore,filename=arrdir+'dtover.sav'
    restore,filename=arrdir+'dtlambda.sav'
    ; diffuse fraction
    restore,filename=arrdir+'fclarr.sav'
    restore,filename=arrdir+'fgmcarr.sav' ; ; add
    ; fcl_test=10.^interpol(alog10(star_flux_frac_arr),alog10(lambdaarr),alog10(lambda))
    ; fgmc_test=10.^interpol(alog10(gas_flux_frac_arr),alog10(lambdaarr),alog10(lambda))
    ;

    ;get cumulative PDF for use as 1D generation function
    probtgascum=total(probtgas*dtgas,/cumulative)

    ;create Monte Carlo arrays
    tgasmc=dblarr(nphysmc)
    tovermc=dblarr(nphysmc)
    lambdamc=dblarr(nphysmc)
    tstarisomc=dblarr(nphysmc)
    if complete then begin
        surfsfrmc=dblarr(nphysmc)
        surfgasmc=dblarr(nphysmc)
    endif

    tstarmc=dblarr(nphysmc)+tstar
    ttotalmc=dblarr(nphysmc)+ttotal
    betastarmc=dblarr(nphysmc)+betastar
    betagasmc=dblarr(nphysmc)+betagas
    surfglobalsmc=dblarr(nphysmc)+surfglobals_fit
    surfglobalgmc=dblarr(nphysmc)+surfglobalg_fit
    surfsmc=dblarr(nphysmc)+surfs_fit
    surfgmc=dblarr(nphysmc)+surfg_fit
    zetastarmc=dblarr(nphysmc)+zetastar
    zetagasmc=dblarr(nphysmc)+zetagas
    rpeakstarmc=dblarr(nphysmc)+rpeakstar
    rpeakgasmc=dblarr(nphysmc)+rpeakgas
    vfbmc=dblarr(nphysmc)+vfb
    if complete then begin
        tdeplmc=dblarr(nphysmc)+tdepl
        esfmc=dblarr(nphysmc)+esf
        mdotsfmc=dblarr(nphysmc)+mdotsf
        mdotfbmc=dblarr(nphysmc)+mdotfb
        etainstmc=dblarr(nphysmc)+etainst
        etaavgmc=dblarr(nphysmc)+etaavg
        chiemc=dblarr(nphysmc)+chie
        chipmc=dblarr(nphysmc)+chip
    endif
    ; diffuse fraction
    fclmc=dblarr(nphysmc)+fcl   ; monte-carlo array for fcl
    fgmcmc=dblarr(nphysmc)+fgmc ; monte-carlo array for fgmc
    ; diffuse fraction corrections
    if rpeak_cor_mode eq 1 then begin
      rpeakscormc=dblarr(nphysmc)
      rpeakgcormc=dblarr(nphysmc)
    endif
    qconstarmc = dblarr(nphysmc)
    qcongasmc = dblarr(nphysmc)
    etastarmc = dblarr(nphysmc)
    etagasmc = dblarr(nphysmc)
    qzetastarmc = dblarr(nphysmc)
    qzetagasmc = dblarr(nphysmc)


    nredo=0
    nredomax=floor(nphysmc/3.)
    ; nredomax=floor(nphysmc*10.)

    nrnd=3*(nphysmc+nredomax)
    randomarr=randomu(systime(1),nrnd)
    gaussarr=randomn(systime(1),nrnd)
    irnd=long(0)
    igauss=long(0)

    tstariso_rerr=0.5*(tstariso_rerrmin+tstariso_rerrmax)
    if rpeak_cor_mode eq 1 then begin
      rpeaks_cor_remin = rpeaks_cor_emin / rpeaks_cor_val
      rpeaks_cor_remax = rpeaks_cor_emax / rpeaks_cor_val
       rpeakscor_rerr=0.5*(rpeaks_cor_remin+rpeaks_cor_remax)


       rpeakg_cor_remin = rpeakg_cor_emin / rpeakg_cor_val
       rpeakg_cor_remax = rpeakg_cor_emax / rpeakg_cor_val
       rpeakgcor_rerr=0.5*(rpeakg_cor_remin+rpeakg_cor_remax)
     endif
    for i=0L,nphysmc-1 do begin ;do Monte Carlo error propagation
        redo=0
        tstarisomc(i)=tstariso+gaussarr(igauss)*tstariso_rerr*tstariso
        xi=tstarisomc(i)/tstariso
        igauss+=1
        tgasmc(i)=10.^interpol(alog10(tgasarr),probtgascum,randomarr(irnd))*xi
        i1=max(where(tgasarr le tgasmc(i),ct))
        if ct eq 0 then i1=0
        if i1 eq ntry-1 then i1=ntry-2
        i2=i1+1
        frac12=alog10(tgasmc(i)/tgasarr(i1))/alog10(tgasarr(i2)/tgasarr(i1))
        probtovertemp=probtgastover(i1,*)+frac12*(probtgastover(i2,*)-probtgastover(i1,*))
        probtovercumtemp=total(probtovertemp*dtover,/cumulative)
        probtovercumtemp=probtovercumtemp/max(probtovercumtemp)
        irnd+=1
        tovermc(i)=10.^interpol(alog10(toverarr),probtovercumtemp,randomarr(irnd))*xi
        i3=max(where(toverarr le tovermc(i),ct))
        if ct eq 0 then i3=0
        if i3 eq ntry-1 then i3=ntry-2
        i4=i3+1
        frac34=alog10(tovermc(i)/toverarr(i3))/alog10(toverarr(i4)/toverarr(i3))
        problambdatemp1=probnorm(i1,i3,*)+frac34*(probnorm(i1,i4,*)-probnorm(i1,i3,*))
        problambdatemp2=probnorm(i2,i3,*)+frac34*(probnorm(i2,i4,*)-probnorm(i2,i3,*))
        problambdatemp=problambdatemp1+frac12*(problambdatemp2-problambdatemp1)
        problambdacumtemp=total(problambdatemp*dlambda,/cumulative)
        problambdacumtemp=problambdacumtemp/max(problambdacumtemp)
        irnd+=1
        lambdamc(i)=10.^interpol(alog10(lambdaarr),problambdacumtemp,randomarr(irnd))
        irnd+=1
        if complete then begin
            surfsfrmc(i)=surfsfr+gaussarr(igauss)*surfsfr_err
            igauss+=1
            surfgasmc(i)=surfgas+gaussarr(igauss)*surfgas_err
            igauss+=1
        endif

        if tstar_incl eq 0 then tstarmc(i)=tstariso*xi+tovermc(i) else tstarmc(i)=tstariso*xi
        ttotalmc(i)=tgasmc(i)+tstarmc(i)-tovermc(i)
        betastarmc(i)=interpol(beta_star,fstarover,tovermc(i)/tstarmc(i))
        betagasmc(i)=interpol(beta_gas,fgasover,tovermc(i)/tgasmc(i))
        surfglobalsmc(i)=10.^interpol(alog10(surfglobals),alog10(apertures_star),alog10(lambdamc(i)))
        surfglobalgmc(i)=10.^interpol(alog10(surfglobalg),alog10(apertures_gas),alog10(lambdamc(i)))
        surfsmc(i)=10.^interpol(alog10(surfcontrasts),alog10(apertures_star),alog10(lambdamc(i)))
        surfgmc(i)=10.^interpol(alog10(surfcontrastg),alog10(apertures_gas),alog10(lambdamc(i)))
        zetastarmc(i)=f_zeta(surfsmc(i),ttotalmc(i),tstarmc(i),peak_prof)
        zetagasmc(i)=f_zeta(surfgmc(i),ttotalmc(i),tgasmc(i),peak_prof)
        rpeakstarmc(i)=f_rpeak(zetastarmc(i),lambdamc(i))
        rpeakgasmc(i)=f_rpeak(zetagasmc(i),lambdamc(i))
        vfbmc(i)=f_vfb(tovermc(i),lambdamc(i))
        if min(finite([tstarisomc(i),tgasmc(i),tovermc(i),lambdamc(i),tstarmc(i),ttotalmc(i),betastarmc(i),betagasmc(i),rpeakstarmc(i),rpeakgasmc(i),zetastarmc(i),zetagasmc(i),vfbmc(i)])) eq 0 then redo=1
        if complete then begin
            tdeplmc(i)=surfgasmc(i)/surfsfrmc(i)/1.d9
            esfmc(i)=f_esf(tgasmc(i),tdeplmc(i),fcl,fgmc)
            mdotsfmc(i)=f_mdotsf(tgasmc(i),lambda,surfgasmc(i)*surfglobalgmc(i),fgmc,esfmc(i))
            mdotfbmc(i)=f_mdotfb(tovermc(i),lambda,surfgasmc(i)*surfglobalgmc(i),fgmc,esfmc(i))
            etainstmc(i)=f_etainst(tgasmc(i),tovermc(i),esfmc(i))
            etaavgmc(i)=f_etaavg(esfmc(i))
            chiemc(i)=f_chie(tovermc(i),esfmc(i),vfbmc(i),psie)
            chipmc(i)=f_chip(tovermc(i),esfmc(i),vfbmc(i),psip)
            if min(finite([surfsfrmc(i),surfgasmc(i),tdeplmc(i),esfmc(i),mdotsfmc(i),mdotfbmc(i),etainstmc(i),etaavgmc(i),chiemc(i),chipmc(i)])) eq 0 then redo=1
        endif
        ; diffuse fraction
        ; fclmc[i]=interpol(star_flux_frac_arr,lambdaarr,lambdamc[i]) ; linear interpolation
        ; fgmcmc[i]=interpol(gas_flux_frac_arr,lambdaarr,lambdamc[i])
        fclmc[i]=10.^interpol(alog10(star_flux_frac_arr),alog10(lambdaarr),alog10(lambdamc[i]))  ; logspace interpolation
        fgmcmc[i]=10.^interpol(alog10(gas_flux_frac_arr),alog10(lambdaarr),alog10(lambdamc[i]))

        if rpeak_cor_mode eq 1 then begin
          rpeakscormc[i]=rpeaks_cor_val+gaussarr(igauss)*rpeakscor_rerr*rpeaks_cor_val
          rpeakgcormc[i]=rpeakg_cor_val+gaussarr(igauss)*rpeakgcor_rerr*rpeakg_cor_val
        endif

        if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
          cut_len_temp = filter_len_conv * lambdamc[i]

          ; if rpeak_cor_mode eq 0 then begin ; use measured rpeak value
          ;   fwhm_star_temp =  rpeakstar * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          ;   fwhm_gas_temp =  rpeakgas * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          ; endif else if rpeak_cor_mode eq 1 then begin ; use supplied rpeak value
          ;   fwhm_star_temp =  rpeaks_cor_val * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          ;   fwhm_gas_temp =  rpeakg_cor_val * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          ; endif

          if rpeak_cor_mode eq 0 then begin ; use measured rpeak value
            fwhm_star_temp = rpeakstarmc[i] * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
            fwhm_gas_temp = rpeakgasmc[i] * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          endif else if rpeak_cor_mode eq 1 then begin ; use supplied rpeak value
            fwhm_star_temp = rpeakscormc[i] * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
            fwhm_gas_temp = rpeakgcormc[i] * (2.0d0 *sqrt(2.0d0 * alog(2.0d0)) )
          endif

          qconstarmc[i] = fourier_flux_loss_correction(filter_choice, cut_len_temp/fwhm_star_temp)
          qcongasmc[i] = fourier_flux_loss_correction(filter_choice, cut_len_temp/fwhm_gas_temp)

          fclmc[i] /= qconstarmc[i] ; flux loss correction
          fgmcmc[i] /= qcongasmc[i]  ; flux loss correction
        endif else begin
          qconstarmc[i] = 1.0d0
          qcongasmc[i] = 1.0d0
        endelse
        if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin

          if rpeak_cor_mode eq 0 then begin ; use measured rpeak value
            etastarmc[i] = sqrt(tstarmc[i]/ttotalmc[i]) * zetastarmc[i] ; corrected zeta for timeline
            etagasmc[i] = sqrt(tgasmc[i]/ttotalmc [i]) * zetagasmc[i]
          endif else if rpeak_cor_mode eq 1 then begin ; use supplied rpeak value
            ; old model:
            ; etastarmc[i] = (sqrt(tstarmc[i]/ttotalmc[i]) * zetastarmc[i]) / (rpeakstarmc[i]/rpeakscormc[i]) ; eta with correction for flux loss from signal regions
            ; etagasmc[i] = (sqrt(tgasmc[i]/ttotalmc[i]) * zetagasmc[i])  / (rpeakgasmc[i]/rpeakgcormc[i]) ; eta with correction for flux loss from signal regions

            ; new model:
            ; etastar = (sqrt(tstar/ttotal) * 2rpeak/lambda)


            etastarmc[i] = (sqrt(tstarmc[i]/ttotalmc[i]) * 2.0 * (rpeakstarmc[i]/lambdamc[i]) ) ; eta with correction for flux loss from signal regions
            etagasmc[i] = (sqrt(tgasmc[i]/ttotalmc[i]) * 2.0 * (rpeakgasmc[i]/lambdamc[i])) ; eta with correction for flux loss from signal regions

          endif

          qzetastarmc[i] = fourier_zeta_correction(filter_choice, etastarmc[i])
          qzetagasmc[i] = fourier_zeta_correction(filter_choice, etagasmc[i])

          fclmc[i]  /=  qzetastarmc[i] ; zeta correction
          fgmcmc[i]  /= qzetagasmc[i] ; zeta correction
        endif else begin
          qzetastarmc[i] = 1.0d0
          qzetagasmc[i] = 1.0d0
        endelse

        if min(finite([fclmc[i], fgmcmc[i], qconstarmc[i], qcongasmc[i], etastarmc[i], etagasmc[i], qzetastarmc[i], qzetagasmc[i] ])) eq 0 then redo=1


        if redo eq 1 then begin ;if any Infinity or NaN has been found, redo draw
            i=i-1
            nredo+=1
            if nredo eq nredomax then f_error,['too many Infinity or NaN values encountered in Monte-Carlo error propagation','chi^2 landscape is too sparsely sampled to converge, please verify quality of best fit']
        endif

        progress,'     ==> derived quantity error propagation progress',i,nphysmc-1
    endfor

    if complete then begin
        quantmc=[[tgasmc],[tovermc],[lambdamc], $
                 [tstarmc],[ttotalmc],[betastarmc],[betagasmc], $
                 [surfglobalsmc],[surfglobalgmc],[surfsmc],[surfgmc],[rpeakstarmc],[rpeakgasmc],[zetastarmc],[zetagasmc],[vfbmc], $
                 [surfsfrmc],[surfgasmc],[tdeplmc],[esfmc],[mdotsfmc],[mdotfbmc],[etainstmc],[etaavgmc],[chiemc],[chipmc], $
                 [fclmc],[fgmcmc],[qconstarmc],[qcongasmc],[etastarmc],[etagasmc],[qzetastarmc],[qzetagasmc]]

    endif else begin
        quantmc=[[tgasmc],[tovermc],[lambdamc], $
                 [tstarmc],[ttotalmc],[betastarmc],[betagasmc], $
                 [surfglobalsmc],[surfglobalgmc],[surfsmc],[surfgmc],[rpeakstarmc],[rpeakgasmc],[zetastarmc],[zetagasmc],[vfbmc], $
                 [fclmc],[fgmcmc],[qconstarmc],[qcongasmc],[etastarmc],[etagasmc],[qzetastarmc],[qzetagasmc]]
    endelse
    ; if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then quantmc=[[quantmc], [qconstarmc],[qcongasmc]]
    ; if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then  quantmc=[[quantmc], [qzetastarmc],[qzetagasmc]]
    ;





    nquant=n_elements(quantmc(0,*))
    corrquant=dblarr(nquant,nquant)
    covquant=dblarr(nquant,nquant)
    for i=0L,nquant-1 do begin
        for j=0L,nquant-1 do begin
            corrquant(i,j)=correlate(quantmc(*,i),quantmc(*,j))
            covquant(i,j)=correlate(quantmc(*,i),quantmc(*,j),/covariance)
        endfor
        progress,'     ==> correlation and covariance matrix progress',i,nquant-1
    endfor
    report=f_writecorr(corrquant,complete,galaxy,outputdir)
    report=f_writecov(covquant,complete,galaxy,outputdir)
    ;sort all arrays
    tgasmc=tgasmc(sort(tgasmc))
    tovermc=tovermc(sort(tovermc))
    lambdamc=lambdamc(sort(lambdamc))
    tstarisomc=tstarisomc(sort(tstarisomc))
    if complete then begin
        surfsfrmc=surfsfrmc(sort(surfsfrmc))
        surfgasmc=surfgasmc(sort(surfgasmc))
    endif

    print,'     ==> calculating derived quantity PDFs and uncertainties'
    tstarmc=tstarmc(sort(tstarmc))
    dtstarmc=f_getdvar(tstarmc)
    ttotalmc=ttotalmc(sort(ttotalmc))
    dttotalmc=f_getdvar(ttotalmc)
    betastarmc=betastarmc(sort(betastarmc))
    dbetastarmc=f_getdvar(betastarmc)
    betagasmc=betagasmc(sort(betagasmc))
    dbetagasmc=f_getdvar(betagasmc)
    surfglobalsmc=surfglobalsmc(sort(surfglobalsmc))
    dsurfglobalsmc=f_getdvar(surfglobalsmc)
    surfglobalgmc=surfglobalgmc(sort(surfglobalgmc))
    dsurfglobalgmc=f_getdvar(surfglobalgmc)
    surfsmc=surfsmc(sort(surfsmc))
    dsurfsmc=f_getdvar(surfsmc)
    surfgmc=surfgmc(sort(surfgmc))
    dsurfgmc=f_getdvar(surfgmc)
    zetastarmc=zetastarmc(sort(zetastarmc))
    dzetastarmc=f_getdvar(zetastarmc)
    zetagasmc=zetagasmc(sort(zetagasmc))
    dzetagasmc=f_getdvar(zetagasmc)
    rpeakstarmc=rpeakstarmc(sort(rpeakstarmc))
    drpeakstarmc=f_getdvar(rpeakstarmc)
    rpeakgasmc=rpeakgasmc(sort(rpeakgasmc))
    drpeakgasmc=f_getdvar(rpeakgasmc)
    vfbmc=vfbmc(sort(vfbmc))
    dvfbmc=f_getdvar(vfbmc)
    if complete then begin
        tdeplmc=tdeplmc(sort(tdeplmc))
        dtdeplmc=f_getdvar(tdeplmc)
        esfmc=esfmc(sort(esfmc))
        desfmc=f_getdvar(esfmc)
        mdotsfmc=mdotsfmc(sort(mdotsfmc))
        dmdotsfmc=f_getdvar(mdotsfmc)
        mdotfbmc=mdotfbmc(sort(mdotfbmc))
        dmdotfbmc=f_getdvar(mdotfbmc)
        etainstmc=etainstmc(sort(etainstmc))
        detainstmc=f_getdvar(etainstmc)
        etaavgmc=etaavgmc(sort(etaavgmc))
        detaavgmc=f_getdvar(etaavgmc)
        chiemc=chiemc(sort(chiemc))
        dchiemc=f_getdvar(chiemc)
        chipmc=chipmc(sort(chipmc))
        dchipmc=f_getdvar(chipmc)
    endif
    ; diffuse fraction
    fclmc=fclmc[sort(fclmc)]
    dfclmc=f_getdvar(fclmc)
    fgmcmc=fgmcmc[sort(fgmcmc)]
    dfgmcmc=f_getdvar(fgmcmc)
    ;
    ; if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
    qconstarmc = qconstarmc[sort(qconstarmc)]
    dqconstarmc = f_getdvar(qconstarmc)
    qcongasmc = qcongasmc[sort(qcongasmc)]
    dqcongasmc = f_getdvar(qcongasmc)
    ; endif
    ; if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin
    etastarmc = etastarmc[sort(etastarmc)]
    detastarmc = f_getdvar(etastarmc)
    etagasmc = etagasmc[sort(etagasmc)]
    detagasmc = f_getdvar(etagasmc)


    qzetastarmc = qzetastarmc[sort(qzetastarmc)]
    dqzetastarmc = f_getdvar(qzetastarmc)
    qzetagasmc = qzetagasmc[sort(qzetagasmc)]
    dqzetagasmc = f_getdvar(qzetagasmc)
    ; endif
    ;get error bars on best-fitting values
    cumdistr=findgen(nphysmc)/(nphysmc-1.)
    dummy=f_pdftovalues(tstarmc,dtstarmc,cumdistr,tstar)
    tstar_errmin=dummy(1)
    tstar_errmax=dummy(2)
    dummy=f_pdftovalues(ttotalmc,dttotalmc,cumdistr,ttotal)
    ttotal_errmin=dummy(1)
    ttotal_errmax=dummy(2)
    dummy=f_pdftovalues(betastarmc,dbetastarmc,cumdistr,betastar)
    betastar_errmin=dummy(1)
    betastar_errmax=dummy(2)
    dummy=f_pdftovalues(betagasmc,dbetagasmc,cumdistr,betagas)
    betagas_errmin=dummy(1)
    betagas_errmax=dummy(2)
    dummy=f_pdftovalues(surfglobalsmc,dsurfglobalsmc,cumdistr,surfglobals_fit)
    surfglobals_errmin=dummy(1)
    surfglobals_errmax=dummy(2)
    dummy=f_pdftovalues(surfglobalgmc,dsurfglobalgmc,cumdistr,surfglobalg_fit)
    surfglobalg_errmin=dummy(1)
    surfglobalg_errmax=dummy(2)
    dummy=f_pdftovalues(surfsmc,dsurfsmc,cumdistr,surfs_fit)
    surfs_errmin=dummy(1)
    surfs_errmax=dummy(2)
    dummy=f_pdftovalues(surfgmc,dsurfgmc,cumdistr,surfg_fit)
    surfg_errmin=dummy(1)
    surfg_errmax=dummy(2)
    dummy=f_pdftovalues(zetastarmc,dzetastarmc,cumdistr,zetastar)
    zetastar_errmin=dummy(1)
    zetastar_errmax=dummy(2)
    dummy=f_pdftovalues(zetagasmc,dzetagasmc,cumdistr,zetagas)
    zetagas_errmin=dummy(1)
    zetagas_errmax=dummy(2)
    dummy=f_pdftovalues(rpeakstarmc,drpeakstarmc,cumdistr,rpeakstar)
    rpeakstar_errmin=dummy(1)
    rpeakstar_errmax=dummy(2)
    dummy=f_pdftovalues(rpeakgasmc,drpeakgasmc,cumdistr,rpeakgas)
    rpeakgas_errmin=dummy(1)
    rpeakgas_errmax=dummy(2)
    dummy=f_pdftovalues(vfbmc,dvfbmc,cumdistr,vfb)
    vfb_errmin=dummy(1)
    vfb_errmax=dummy(2)
    if complete then begin
        dummy=f_pdftovalues(tdeplmc,dtdeplmc,cumdistr,tdepl)
        tdepl_errmin=dummy(1)
        tdepl_errmax=dummy(2)
        dummy=f_pdftovalues(esfmc,desfmc,cumdistr,esf)
        esf_errmin=dummy(1)
        esf_errmax=dummy(2)
        dummy=f_pdftovalues(mdotsfmc,dmdotsfmc,cumdistr,mdotsf)
        mdotsf_errmin=dummy(1)
        mdotsf_errmax=dummy(2)
        dummy=f_pdftovalues(mdotfbmc,dmdotfbmc,cumdistr,mdotfb)
        mdotfb_errmin=dummy(1)
        mdotfb_errmax=dummy(2)
        dummy=f_pdftovalues(etainstmc,detainstmc,cumdistr,etainst)
        etainst_errmin=dummy(1)
        etainst_errmax=dummy(2)
        dummy=f_pdftovalues(etaavgmc,detaavgmc,cumdistr,etaavg)
        etaavg_errmin=dummy(1)
        etaavg_errmax=dummy(2)
        dummy=f_pdftovalues(chiemc,dchiemc,cumdistr,chie)
        chie_errmin=dummy(1)
        chie_errmax=dummy(2)
        dummy=f_pdftovalues(chipmc,dchipmc,cumdistr,chip)
        chip_errmin=dummy(1)
        chip_errmax=dummy(2)
    endif
    ; diffuse fraction
    dummy=f_pdftovalues(fclmc,dfclmc,cumdistr,fcl) ; get errors for fcl
    fcl_errmin=dummy[1]
    fcl_errmax=dummy[2]
    dummy=f_pdftovalues(fgmcmc,dfgmcmc,cumdistr,fgmc)
    fgmc_errmin=dummy[1]
    fgmc_errmax=dummy[2]

    ;
    ; if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
    dummy=f_pdftovalues(qconstarmc,dqconstarmc,cumdistr,qconstar) ; get errors for qconstar
    qconstar_errmin=dummy[1]
    qconstar_errmax=dummy[2]
    dummy=f_pdftovalues(qcongasmc,dqcongasmc,cumdistr,qcongas) ; get errors for qcongas
    qcongas_errmin=dummy[1]
    qcongas_errmax=dummy[2]
    ; endif
    ; if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin
    dummy=f_pdftovalues(etastarmc,detastarmc,cumdistr,etastar) ; get errors for qzetastar
    etastar_errmin=dummy[1]
    etastar_errmax=dummy[2]
    dummy=f_pdftovalues(etagasmc,detagasmc,cumdistr,etagas) ; get errors for qzetagas
    etagas_errmin=dummy[1]
    etagas_errmax=dummy[2]

    dummy=f_pdftovalues(qzetastarmc,dqzetastarmc,cumdistr,qzetastar) ; get errors for qzetastar
    qzetastar_errmin=dummy[1]
    qzetastar_errmax=dummy[2]
    dummy=f_pdftovalues(qzetagasmc,dqzetagasmc,cumdistr,qzetagas) ; get errors for qzetagas
    qzetagas_errmin=dummy[1]
    qzetagas_errmax=dummy[2]
    ; endif



    ;plot PDFs and write tables
    print,'     ==> plotting PDFs and writing them to output directory'
    dummy=f_createpdf(tstarmc,dtstarmc,cumdistr,ntry)
    tstararr=dummy(*,0)
    dtstar=dummy(*,1)
    probtstar=dummy(*,2)
    report=f_plotdistr(tstararr,dtstar,probtstar,tstar,tstar_errmin,tstar_errmax,galaxy,figdir,'tstar','!8t!6!Dstar!N','!6Myr',0)
    report=f_writepdf(alog10(tstararr),alog10(dtstar),alog10(probtstar),galaxy,outputdir,'tstar','# log10(tstar[Myr]), log10(dtstar[Myr]), log10(PDF[Myr^-1])')
    dummy=f_createpdf(ttotalmc,dttotalmc,cumdistr,ntry)
    ttotalarr=dummy(*,0)
    dttotal=dummy(*,1)
    probttotal=dummy(*,2)
    report=f_plotdistr(ttotalarr,dttotal,probttotal,ttotal,ttotal_errmin,ttotal_errmax,galaxy,figdir,'ttotal','!7s!6','!6Myr',0)
    report=f_writepdf(alog10(ttotalarr),alog10(dttotal),alog10(probttotal),galaxy,outputdir,'ttotal','# log10(ttotal[Myr]), log10(dttotal[Myr]), log10(PDF[Myr^-1])')
    dummy=f_createpdf(betastarmc,dbetastarmc,cumdistr,ntry)
    betastararr=dummy(*,0)
    dbetastar=dummy(*,1)
    probbetastar=dummy(*,2)
    report=f_plotdistr(betastararr,dbetastar,probbetastar,betastar,betastar_errmin,betastar_errmax,galaxy,figdir,'betastar','!7b!6!Dstar!N','',0)
    report=f_writepdf(alog10(betastararr),alog10(dbetastar),alog10(probbetastar),galaxy,outputdir,'betastar','# log10(betastar), log10(dbetastar), log10(PDF)')
    dummy=f_createpdf(betagasmc,dbetagasmc,cumdistr,ntry)
    betagasarr=dummy(*,0)
    dbetagas=dummy(*,1)
    probbetagas=dummy(*,2)
    report=f_plotdistr(betagasarr,dbetagas,probbetagas,betagas,betagas_errmin,betagas_errmax,galaxy,figdir,'betagas','!7b!6!Dgas!N','',0)
    report=f_writepdf(alog10(betagasarr),alog10(dbetagas),alog10(probbetagas),galaxy,outputdir,'betagas','# log10(betagas), log10(dbetagas), log10(PDF)')
    dummy=f_createpdf(surfglobalsmc,dsurfglobalsmc,cumdistr,ntry)
    surfglobalsarr=dummy(*,0)
    dsurfglobals=dummy(*,1)
    probsurfglobals=dummy(*,2)
    report=f_plotdistr(surfglobalsarr,dsurfglobals,probsurfglobals,surfglobals_fit,surfglobals_errmin,surfglobals_errmax,galaxy,figdir,'surfglobalstar','!7E!6!Dstar,global!N','',0)
    report=f_writepdf(alog10(surfglobalsarr),alog10(dsurfglobals),alog10(probsurfglobals),galaxy,outputdir,'surfglobalstar','# log10(surfglobalstar), log10(dsurfglobalstar), log10(PDF)')
    dummy=f_createpdf(surfglobalgmc,dsurfglobalgmc,cumdistr,ntry)
    surfglobalgarr=dummy(*,0)
    dsurfglobalg=dummy(*,1)
    probsurfglobalg=dummy(*,2)
    report=f_plotdistr(surfglobalgarr,dsurfglobalg,probsurfglobalg,surfglobalg_fit,surfglobalg_errmin,surfglobalg_errmax,galaxy,figdir,'surfglobalgas','!7E!6!Dgas,global!N','',0)
    report=f_writepdf(alog10(surfglobalgarr),alog10(dsurfglobalg),alog10(probsurfglobalg),galaxy,outputdir,'surfglobalgas','# log10(surfglobalgas), log10(dsurfglobalgas), log10(PDF)')
    dummy=f_createpdf(surfsmc,dsurfsmc,cumdistr,ntry)
    surfsarr=dummy(*,0)
    dsurfs=dummy(*,1)
    probsurfs=dummy(*,2)
    report=f_plotdistr(surfsarr,dsurfs,probsurfs,surfs_fit,surfs_errmin,surfs_errmax,galaxy,figdir,'surfconstar','!7E!6!Dstar!N','',0)
    report=f_writepdf(alog10(surfsarr),alog10(dsurfs),alog10(probsurfs),galaxy,outputdir,'surfconstar','# log10(surfconstar), log10(dsurfconstar), log10(PDF)')
    dummy=f_createpdf(surfgmc,dsurfgmc,cumdistr,ntry)
    surfgarr=dummy(*,0)
    dsurfg=dummy(*,1)
    probsurfg=dummy(*,2)
    report=f_plotdistr(surfgarr,dsurfg,probsurfg,surfg_fit,surfg_errmin,surfg_errmax,galaxy,figdir,'surfcongas','!7E!6!Dgas!N','',0)
    report=f_writepdf(alog10(surfgarr),alog10(dsurfg),alog10(probsurfg),galaxy,outputdir,'surfcongas','# log10(surfcongas), log10(dsurfcongas), log10(PDF)')
    dummy=f_createpdf(zetastarmc,dzetastarmc,cumdistr,ntry)
    zetastararr=dummy(*,0)
    dzetastar=dummy(*,1)
    probzetastar=dummy(*,2)
    report=f_plotdistr(zetastararr,dzetastar,probzetastar,zetastar,zetastar_errmin,zetastar_errmax,galaxy,figdir,'zetastar','!7f!6!Dstar!N','',0)
    report=f_writepdf(alog10(zetastararr),alog10(dzetastar),alog10(probzetastar),galaxy,outputdir,'zetastar','# log10(zetastar), log10(dzetastar), log10(PDF)')
    dummy=f_createpdf(zetagasmc,dzetagasmc,cumdistr,ntry)
    zetagasarr=dummy(*,0)
    dzetagas=dummy(*,1)
    probzetagas=dummy(*,2)
    report=f_plotdistr(zetagasarr,dzetagas,probzetagas,zetagas,zetagas_errmin,zetagas_errmax,galaxy,figdir,'zetagas','!7f!6!Dgas!N','',0)
    report=f_writepdf(alog10(zetagasarr),alog10(dzetagas),alog10(probzetagas),galaxy,outputdir,'zetagas','# log10(zetagas), log10(dzetagas), log10(PDF)')
    dummy=f_createpdf(rpeakstarmc,drpeakstarmc,cumdistr,ntry)
    rpeakstararr=dummy(*,0)
    drpeakstar=dummy(*,1)
    probrpeakstar=dummy(*,2)
    report=f_plotdistr(rpeakstararr,drpeakstar,probrpeakstar,rpeakstar,rpeakstar_errmin,rpeakstar_errmax,galaxy,figdir,'rpeakstar','!8r!6!Dstar!N','!6pc',0)
    report=f_writepdf(alog10(rpeakstararr),alog10(drpeakstar),alog10(probrpeakstar),galaxy,outputdir,'rpeakstar','# log10(rpeakstar[pc]), log10(drpeakstar[pc]), log10(PDF[pc^-1])')
    dummy=f_createpdf(rpeakgasmc,drpeakgasmc,cumdistr,ntry)
    rpeakgasarr=dummy(*,0)
    drpeakgas=dummy(*,1)
    probrpeakgas=dummy(*,2)
    report=f_plotdistr(rpeakgasarr,drpeakgas,probrpeakgas,rpeakgas,rpeakgas_errmin,rpeakgas_errmax,galaxy,figdir,'rpeakgas','!8r!6!Dgas!N','!6pc',0)
    report=f_writepdf(alog10(rpeakgasarr),alog10(drpeakgas),alog10(probrpeakgas),galaxy,outputdir,'rpeakgas','# log10(rpeakgas[pc]), log10(drpeakgas[pc]), log10(PDF[pc^-1])')
    dummy=f_createpdf(vfbmc,dvfbmc,cumdistr,ntry)
    vfbarr=dummy(*,0)
    dvfb=dummy(*,1)
    probvfb=dummy(*,2)
    report=f_plotdistr(vfbarr,dvfb,probvfb,vfb,vfb_errmin,vfb_errmax,galaxy,figdir,'vfb','!8v!6!Dfb!N','!6km s!U-1!N',0)
    report=f_writepdf(alog10(vfbarr),alog10(dvfb),alog10(probvfb),galaxy,outputdir,'vfb','# log10(vfb[km s^-1]), log10(dvfb[km s^-1]), log10(PDF[km^-1 s])')
    if complete then begin
        dummy=f_createpdf(tdeplmc,dtdeplmc,cumdistr,ntry)
        tdeplarr=dummy(*,0)
        dtdepl=dummy(*,1)
        probtdepl=dummy(*,2)
        report=f_plotdistr(tdeplarr,dtdepl,probtdepl,tdepl,tdepl_errmin,tdepl_errmax,galaxy,figdir,'tdepl','!8t!6!Ddepl!N','!6Gyr',0)
        report=f_writepdf(alog10(tdeplarr),alog10(dtdepl),alog10(probtdepl),galaxy,outputdir,'tdepl','# log10(tdepl[Gyr]), log10(dtdepl[Gyr]), log10(PDF[Gyr^-1])')
        dummy=f_createpdf(esfmc,desfmc,cumdistr,ntry)
        esfarr=dummy(*,0)
        desf=dummy(*,1)
        probesf=dummy(*,2)
        report=f_plotdistr(esfarr,desf,probesf,esf,esf_errmin,esf_errmax,galaxy,figdir,'esf','!7e!6!Dsf!N','',0)
        report=f_writepdf(alog10(esfarr),alog10(desf),alog10(probesf),galaxy,outputdir,'esf','# log10(esf), log10(desf), log10(PDF)')
        dummy=f_createpdf(mdotsfmc,dmdotsfmc,cumdistr,ntry)
        mdotsfarr=dummy(*,0)
        dmdotsf=dummy(*,1)
        probmdotsf=dummy(*,2)
        report=f_plotdistr(mdotsfarr,dmdotsf,probmdotsf,mdotsf,mdotsf_errmin,mdotsf_errmax,galaxy,figdir,'mdotsf','!6(d!8M!6/d!8t!6)!Dsf!N','!6M!D!9n!6!N yr!U-1!N',0)
        report=f_writepdf(alog10(mdotsfarr),alog10(dmdotsf),alog10(probmdotsf),galaxy,outputdir,'mdotsf','# log10(mdotsf[Msun yr^-1]), log10(dmdotsf[Msun yr^-1]), log10(PDF[Msun^-1 yr])')
        dummy=f_createpdf(mdotfbmc,dmdotfbmc,cumdistr,ntry)
        mdotfbarr=dummy(*,0)
        dmdotfb=dummy(*,1)
        probmdotfb=dummy(*,2)
        report=f_plotdistr(mdotfbarr,dmdotfb,probmdotfb,mdotfb,mdotfb_errmin,mdotfb_errmax,galaxy,figdir,'mdotfb','!6(d!8M!6/d!8t!6)!Dfb!N','!6M!D!9n!6!N yr!U-1!N',0)
        report=f_writepdf(alog10(mdotfbarr),alog10(dmdotfb),alog10(probmdotfb),galaxy,outputdir,'mdotfb','# log10(mdotfb[Msun yr^-1]), log10(dmdotfb[Msun yr^-1]), log10(PDF[Msun^-1 yr])')
        dummy=f_createpdf(etainstmc,detainstmc,cumdistr,ntry)
        etainstarr=dummy(*,0)
        detainst=dummy(*,1)
        probetainst=dummy(*,2)
        report=f_plotdistr(etainstarr,detainst,probetainst,etainst,etainst_errmin,etainst_errmax,galaxy,figdir,'etainst','!7g!6!Dfb!N','',0)
        report=f_writepdf(alog10(etainstarr),alog10(detainst),alog10(probetainst),galaxy,outputdir,'etainst','# log10(etainst), log10(detainst), log10(PDF)')
        dummy=f_createpdf(etaavgmc,detaavgmc,cumdistr,ntry)
        etaavgarr=dummy(*,0)
        detaavg=dummy(*,1)
        probetaavg=dummy(*,2)
        report=f_plotdistr(etaavgarr,detaavg,probetaavg,etaavg,etaavg_errmin,etaavg_errmax,galaxy,figdir,'etaavg','!6<!7g!6!Dfb!N>','',0)
        report=f_writepdf(alog10(etaavgarr),alog10(detaavg),alog10(probetaavg),galaxy,outputdir,'etaavg','# log10(etaavg), log10(detaavg), log10(PDF)')
        dummy=f_createpdf(chiemc,dchiemc,cumdistr,ntry)
        chiearr=dummy(*,0)
        dchie=dummy(*,1)
        probchie=dummy(*,2)
        report=f_plotdistr(chiearr,dchie,probchie,chie,chie_errmin,chie_errmax,galaxy,figdir,'chie','!7v!6!Dfb,!8E!6!N','',0)
        report=f_writepdf(alog10(chiearr),alog10(dchie),alog10(probchie),galaxy,outputdir,'chie','# log10(chie), log10(dchie), log10(PDF)')
        dummy=f_createpdf(chipmc,dchipmc,cumdistr,ntry)
        chiparr=dummy(*,0)
        dchip=dummy(*,1)
        probchip=dummy(*,2)
        report=f_plotdistr(chiparr,dchip,probchip,chip,chip_errmin,chip_errmax,galaxy,figdir,'chip','!7v!6!Dfb,!8p!6!N','',0)
        report=f_writepdf(alog10(chiparr),alog10(dchip),alog10(probchip),galaxy,outputdir,'chip','# log10(chip), log10(dchip), log10(PDF)')
    endif
    ; diffuse fraction
    dummy=f_createpdf(fclmc,dfclmc,cumdistr,ntry)
    fclarr=dummy(*,0)
    dfcl=dummy(*,1)
    probfcl=dummy(*,2)
    report=f_plotdistr(fclarr,dfcl,probfcl,fcl,fcl_errmin,fcl_errmax,galaxy,figdir,'fcl','!8f!6!Dcl!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(fclarr),alog10(dfcl),alog10(probfcl),galaxy,outputdir,'fcl','# log10(fcl), log10(dfcl), log10(PDF)')

    dummy=f_createpdf(fgmcmc,dfgmcmc,cumdistr,ntry)
    fgmcarr=dummy(*,0)
    dfgmc=dummy(*,1)
    probfgmc=dummy(*,2)
    report=f_plotdistr(fgmcarr,dfgmc,probfgmc,fgmc,fgmc_errmin,fgmc_errmax,galaxy,figdir,'fgmc','!8f!6!Dgmc!N','',0);'!7f!6!Dgmc!N' = f_gmc
    report=f_writepdf(alog10(fgmcarr),alog10(dfgmc),alog10(probfgmc),galaxy,outputdir,'fgmc','# log10(fgmc), log10(dfgmc), log10(PDF)')
    ;


    ; if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
    dummy=f_createpdf(qconstarmc,dqconstarmc,cumdistr,ntry)
    qconstararr=dummy(*,0)
    dqconstar=dummy(*,1)
    probqconstar=dummy(*,2)
    report=f_plotdistr(qconstararr,dqconstar,probqconstar,qconstar,qconstar_errmin,qconstar_errmax,galaxy,figdir,'qconstar','!6Q!Dcon,star!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(qconstararr),alog10(dqconstar),alog10(probqconstar),galaxy,outputdir,'qconstar','# log10(qconstar), log10(dqconstar), log10(PDF)')


    dummy=f_createpdf(qcongasmc,dqcongasmc,cumdistr,ntry)
    qcongasarr=dummy(*,0)
    dqcongas=dummy(*,1)
    probqcongas=dummy(*,2)
    report=f_plotdistr(qcongasarr,dqcongas,probqcongas,qcongas,qcongas_errmin,qcongas_errmax,galaxy,figdir,'qcongas','!6Q!Dcon,gas!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(qcongasarr),alog10(dqcongas),alog10(probqcongas),galaxy,outputdir,'qcongas','# log10(qcongas), log10(dqcongas), log10(PDF)')
    ; endif

    ; if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin
    dummy=f_createpdf(etastarmc,detastarmc,cumdistr,ntry)
    etastararr=dummy(*,0)
    detastar=dummy(*,1)
    probetastar=dummy(*,2)
    report=f_plotdistr(etastararr,detastar,probetastar,etastar,etastar_errmin,etastar_errmax,galaxy,figdir,'etastar','!7g!6!Dstar!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(etastararr),alog10(detastar),alog10(probetastar),galaxy,outputdir,'etastar','# log10(etastar), log10(detastar), log10(PDF)')


    dummy=f_createpdf(etagasmc,detagasmc,cumdistr,ntry)
    etagasarr=dummy(*,0)
    detagas=dummy(*,1)
    probetagas=dummy(*,2)
    report=f_plotdistr(etagasarr,detagas,probetagas,etagas,etagas_errmin,etagas_errmax,galaxy,figdir,'etagas','!7g!6!Dgas!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(etagasarr),alog10(detagas),alog10(probetagas),galaxy,outputdir,'etagas','# log10(etagas), log10(detagas), log10(PDF)')

    dummy=f_createpdf(qzetastarmc,dqzetastarmc,cumdistr,ntry)
    qzetastararr=dummy(*,0)
    dqzetastar=dummy(*,1)
    probqzetastar=dummy(*,2)
    report=f_plotdistr(qzetastararr,dqzetastar,probqzetastar,qzetastar,qzetastar_errmin,qzetastar_errmax,galaxy,figdir,'qzetastar','!6Q!D!7f!6,star!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(qzetastararr),alog10(dqzetastar),alog10(probqzetastar),galaxy,outputdir,'qzetastar','# log10(qzetastar), log10(dqzetastar), log10(PDF)')


    dummy=f_createpdf(qzetagasmc,dqzetagasmc,cumdistr,ntry)
    qzetagasarr=dummy(*,0)
    dqzetagas=dummy(*,1)
    probqzetagas=dummy(*,2)
    report=f_plotdistr(qzetagasarr,dqzetagas,probqzetagas,qzetagas,qzetagas_errmin,qzetagas_errmax,galaxy,figdir,'qzetagas','!6Q!D!7f!6,gas!N','',0)  ;'!7f!6!Dcl!N' = f_cl
    report=f_writepdf(alog10(qzetagasarr),alog10(dqzetagas),alog10(probqzetagas),galaxy,outputdir,'qzetagas','# log10(qzetagas), log10(dqzetagas), log10(PDF)')
    ; endif




    returnarr = [tstar,tstar_errmin,tstar_errmax, $
                ttotal,ttotal_errmin,ttotal_errmax, $
                betastar,betastar_errmin,betastar_errmax, $
                betagas,betagas_errmin,betagas_errmax, $
                surfglobals_fit,surfglobals_errmin,surfglobals_errmax, $
                surfglobalg_fit,surfglobalg_errmin,surfglobalg_errmax, $
                surfs_fit,surfs_errmin,surfs_errmax, $
                surfg_fit,surfg_errmin,surfg_errmax, $
                fcl, fcl_errmin, fcl_errmax, $
                fgmc, fgmc_errmin, fgmc_errmax, $
                qconstar,qconstar_errmin,qconstar_errmax, $
                qcongas,qcongas_errmin,qcongas_errmax, $
                etastar,etastar_errmin,etastar_errmax, $
                etagas,etagas_errmin,etagas_errmax, $
                qzetastar,qzetastar_errmin,qzetastar_errmax, $
                qzetagas,qzetagas_errmin,qzetagas_errmax, $
                zetastar,zetastar_errmin,zetastar_errmax, $
                zetagas,zetagas_errmin,zetagas_errmax, $
                rpeakstar,rpeakstar_errmin,rpeakstar_errmax, $
                rpeakgas,rpeakgas_errmin,rpeakgas_errmax, $
                vfb,vfb_errmin,vfb_errmax]

    if complete then returnarr = [returnarr, $
                                  tdepl,tdepl_errmin,tdepl_errmax, $
                                  esf,esf_errmin,esf_errmax, $
                                  mdotsf,mdotsf_errmin,mdotsf_errmax, $
                                  mdotfb,mdotfb_errmin,mdotfb_errmax, $
                                  etainst,etainst_errmin,etainst_errmax, $
                                  etaavg,etaavg_errmin,etaavg_errmax, $
                                  chie,chie_errmin,chie_errmax, $
                                  chip,chip_errmin,chip_errmax]

    ; if emfrac_cor_mode eq 1 || emfrac_cor_mode eq 3  then begin
    ;   returnarr = [returnarr, $
    ;                qconstar,qconstar_errmin,qconstar_errmax, $
    ;                qcongas,qcongas_errmin,qcongas_errmax]
    ; endif
    ; if emfrac_cor_mode eq 2 || emfrac_cor_mode eq 3  then begin
    ;   returnarr = [returnarr, $
    ;                qzetastar,qzetastar_errmin,qzetastar_errmax, $
    ;                qzetagas,qzetagas_errmin,qzetagas_errmax]
    ; endif


    return,returnarr



end
