;;;;;;;;;;;;;;;;;;;;
;FIT KL14 PRINCIPLE;
;;;;;;;;;;;;;;;;;;;;

function f_string,numorig,ndec ;converts any real number to string with up to four decimals
    if ndec gt 4 then stop ;only works with floats and up to four decimals
    arraylen=n_elements(numorig)
    num=fltarr(arraylen)
    numstr=strarr(arraylen)
    for i=0,arraylen-1 do begin
        num(i)=float(10.^(-ndec)*round(10.^ndec*numorig(i)))
        sign=num(i) ge 0.
        lsign=abs(num(i)) ge 1.
        nlog=fix(alog10(abs(num(i))+1.d-5))+lsign-1
        n1=6+(num(i) ne 0.)*((nlog lt 0)*nlog+sign-1)
        n2=ndec+2-(ndec eq 0)+(num(i) ne 0.)*((nlog gt 0)*nlog-sign+1)
        numstr(i)=strmid(num(i),n1,n2)
    endfor
    if arraylen eq 1 then return,strmid(num(0),n1(0),n2(0)) else return,numstr
end

function f_intfrac,var,r1,r2
    a0=var(0) ;Gaussian normalisation
    a1=var(1) ;Gaussian x offset
    a2=var(2) ;Gaussian dispersion
    sigma1=gaussian(r1,[a0,a1,a2]) ;surface density at inner edge
    sigma2=gaussian(r2,[a0,a1,a2]) ;surface density at outer edge
    erf1=erf((r1-a1)/(sqrt(2.)*a2)) ;error function at inner edge
    erf2=erf((r2-a1)/(sqrt(2.)*a2)) ;error function at outer edge
    int1=a2^2.*(sigma1-sigma2)+sqrt(!pi/2.)*a0*a1*a2*(erf2-erf1) ;integral of surface density between inner and outer edge
    int2=a2^2.*sigma1+sqrt(!pi/2.)*a0*a1*a2*(1.-erf1) ;integral of surface density between inner edge and infinity
    int=int1/int2 ;fraction covered within outer edge
    return,int
end

function f_fluxratiostar,tgas,tstar,tover,laps,lambda,beta_gas,surfcontrasts,surfcontrastg,peak_prof
    forward_function F_ZETA
    forward_function F_RPEAK
    ttotal=tgas+tstar-tover
    zetas=f_zeta(surfcontrasts,ttotal,tstar,peak_prof)
    zetag=f_zeta(surfcontrastg,ttotal,tgas,peak_prof)
    rpeaks=f_rpeak(zetas,lambda)
    rpeakg=f_rpeak(zetag,lambda)
    nl=n_elements(laps) ;number of aperture sizes
    surfratios=dblarr(nl)
    surfratiog=dblarr(nl)
    for i=0,nl-1 do begin
        if peak_prof le 1 then begin
            if rpeaks lt 0.5*laps(i) then surfratios(i)=ttotal/tstar*(lambda/laps(i))^2. else surfratios(i)=surfcontrasts
            if rpeakg lt 0.5*laps(i) then surfratiog(i)=ttotal/tgas*(lambda/laps(i))^2. else surfratiog(i)=surfcontrastg
        endif
        if peak_prof eq 2 then begin
            surfratios(i)=ttotal/tstar*(lambda/laps(i))^2.*(1.-exp(-laps(i)^2./8./rpeaks^2.))
            surfratiog(i)=ttotal/tgas*(lambda/laps(i))^2.*(1.-exp(-laps(i)^2./8./rpeakg^2.))
        endif
    endfor
    others=1.
    gas=(beta_gas*tover/tstar/(1.+(beta_gas-1.)*tover/tgas))*surfratiog
    star=surfratios
    tdeplstar=(gas+others)/(star+others)
    return,tdeplstar
end

function f_fluxratiogas,tgas,tstar,tover,lapg,lambda,beta_star,surfcontrasts,surfcontrastg,peak_prof
    forward_function F_ZETA
    forward_function F_RPEAK
    ttotal=tgas+tstar-tover
    zetas=f_zeta(surfcontrasts,ttotal,tstar,peak_prof)
    zetag=f_zeta(surfcontrastg,ttotal,tgas,peak_prof)
    rpeaks=f_rpeak(zetas,lambda)
    rpeakg=f_rpeak(zetag,lambda)
    nl=n_elements(lapg) ;number of aperture sizes
    surfratios=dblarr(nl)
    surfratiog=dblarr(nl)
    for i=0,nl-1 do begin
        if peak_prof le 1 then begin
            if rpeaks lt 0.5*lapg(i) then surfratios(i)=ttotal/tstar*(lambda/lapg(i))^2. else surfratios(i)=surfcontrasts
            if rpeakg lt 0.5*lapg(i) then surfratiog(i)=ttotal/tgas*(lambda/lapg(i))^2. else surfratiog(i)=surfcontrastg
        endif
        if peak_prof eq 2 then begin
            surfratios(i)=ttotal/tstar*(lambda/lapg(i))^2.*(1.-exp(-lapg(i)^2./8./rpeaks^2.))
            surfratiog(i)=ttotal/tgas*(lambda/lapg(i))^2.*(1.-exp(-lapg(i)^2./8./rpeakg^2.))
        endif
    endfor
    others=1.
    gas=surfratiog
    star=(beta_star*tover/tgas/(1.+(beta_star-1.)*tover/tstar))*surfratios
    tdeplgas=(gas+others)/(star+others)
    return,tdeplgas
end

function f_pdftovalues,variable,dvariable,cumpdf,refvalue ;convert PDF to refvalue^+errmax_-errmin, where the errors are set by the 1sigma percentiles relative to refvalue
    COMMON numbers
    errminfrac_gaussian=1.-.5*(1.+erf(1./sqrt(2.)))
    errmaxfrac_gaussian=1.-errminfrac_gaussian
    refprob=interpol(cumpdf,variable,refvalue)
    if refprob lt 0. then refprob=tiny
    if refprob gt 1. then refprob=1.-tiny
    if finite(refprob) eq 0. then refprob=0.5
    errminfrac=2.*errminfrac_gaussian*refprob ;(.5-errmin)/.5=(1-errmin2/p) => errmin2 = 2*errmin*p (32nd percentile below best-fitting value)
    errmaxfrac=2.*(1.-errmaxfrac_gaussian)*refprob+2.*(errmaxfrac_gaussian-.5) ;(errmax-.5)/.5=(errmax2-p)/(1-p) => errmax2 = (2*errmax-1)*(1-p)+p = 2*(1-errmax)*p+2*(errmax-.5) (68th percentile above best-fitting value)
    nvar=n_elements(variable)
    pdf=dblarr(nvar)
    pdf(0)=cumpdf(0)/dvariable(0)
    for i=1,nvar-1 do pdf(i)=(cumpdf(i)-cumpdf(i-1))/dvariable(i)
    evalue=total(variable*pdf*dvariable)
    minperc=interpol(variable,cumpdf,errminfrac)
    maxperc=interpol(variable,cumpdf,errmaxfrac)
    if minperc ne maxperc then begin
        errmin=max([0.,refvalue-minperc])
        errmax=max([0.,maxperc-refvalue])
    endif else begin
        errmin=tiny
        errmax=tiny
    endelse
    return,[evalue,errmin,errmax]
end

function f_writepdf,array,darray,probarray,galaxy,outputdir,varstring,commentstring ;write table with PDF to file
    openw,lun,outputdir+galaxy+'_PDF_'+varstring+'.dat',/get_lun
    printf,lun,'# 1D Probability distribution function of '+varstring+' for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how this PDF was calculated'
    printf,lun,'# IMPORTANT: this PDF uses a coarser grid than the PDF used to calculate the best-fitting values in the PDF figures'
    printf,lun,'# IMPORTANT: rounding errors may therefore cause small differences when using this file to reproduce the PDF figures'
    printf,lun,commentstring
    printf,lun,''
    n=n_elements(array)
    for i=0,n-1 do begin
        ar=array(i)
        dar=darray(i)
        pdf=probarray(i)
        printf,lun,format='(f12.5,2(3x,f12.5))',ar,dar,pdf
    endfor
    close,lun
    free_lun,lun
    report='         1D Probability distribution function of '+varstring+' written'
    return,report
end

function f_writetf_model,xarray,yarraystar,yarraygas,galaxy,outputdir,commentstring ;write table with tuningfork model line to file
    openw,lun,outputdir+galaxy+'_tuningfork_model.dat',/get_lun
    printf,lun,'# Best-fitting model output for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how this model was calculated'
    printf,lun,commentstring
    printf,lun,''
    n=n_elements(xarray)
    for i=0,n-1 do begin
        v1=xarray(i)
        v2=yarraystar(i)
        v3=yarraygas(i)
        printf,lun,format='(f12.5,2(3x,f12.5))',v1,v2,v3
    endfor
    close,lun
    free_lun,lun
    report='         Best-fitting model output written'
    return,report
end

function f_writetf_obs,xarraystar,xarraygas,yarraystar,yarraygas,errbarstar,errbargas,errfitstar,errfitgas,galaxy,outputdir,commentstring ;write table with tuningfork observational data points to file
    openw,lun,outputdir+galaxy+'_tuningfork_obs.dat',/get_lun
    printf,lun,'# Observational data point output for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how these data points were calculated'
    printf,lun,commentstring
    printf,lun,''
    n=n_elements(xarraystar)
    for i=0,n-1 do begin
        v1=xarraystar(i)
        v2=xarraygas(i)
        v3=yarraystar(i)
        v4=yarraygas(i)
        v5=errbarstar(i)
        v6=errbargas(i)
        v7=errfitstar(i)
        v8=errfitgas(i)
        printf,lun,format='(f12.5,7(3x,f12.5))',v1,v2,v3,v4,v5,v6,v7,v8
    endfor
    close,lun
    free_lun,lun
    report='         Observational data point output written'
    return,report
end

function f_plotdistr,array,darray,pdf,value,errmin,errmax,galaxy,figdir,varstring,symstring,unitstring,verbose ;plot PDFs to postscript
    arraymin=min(array)
    arraymax=max(array)
    cumpdf=total(pdf*darray,/cumulative)
    dummy=f_pdftovalues(array,darray,cumpdf,value)
    value_frompdf=dummy(0)
    errmin_frompdf=dummy(1)
    errmax_frompdf=dummy(2)
    if verbose then begin
        print,'         value specified: '+f_string(value,2)
        print,'         errmin specified: '+f_string(errmin,2)
        print,'         errmax specified: '+f_string(errmax,2)
        print,'         value calculated from PDF: '+f_string(value_frompdf,2)
        print,'         errmin calculated from PDF: '+f_string(errmin_frompdf,2)
        print,'         errmax calculated from PDF: '+f_string(errmax_frompdf,2)
        print,'         the above specified and calculated values may differ due to rounding errors from using differently-sized PDF grids'
    endif
    logdiffarray=0.5*(2.-alog10(arraymax/arraymin)) ;is the range at least two orders of magnitude?
    if logdiffarray gt 0 then arrange=[arraymin/10.^logdiffarray,arraymax*10.^logdiffarray] else arrange=[arraymin,arraymax] ;if not, extend for plotting
    if strlen(unitstring) gt 0 then begin
        xlab=symstring+' ['+unitstring+']'
        ylab='!6d!8p!6/d'+symstring+' [('+unitstring+')!U-1!N]'
    endif else begin
        xlab=symstring
        ylab='!6d!8p!6/d'+symstring
    endelse
    origdevice=!D.name
    set_plot,'ps'
    device,filename=figdir+galaxy+'_distr_'+varstring+'.ps',xsize=12,ysize=9,/color,bits_per_pixel=8,/encapsulated
    ymax=max(pdf)
    plot,array,pdf,/nodata,/xlog,xtitle=xlab,ytitle=ylab,xr=arrange,xstyle=1
    oplot,array,pdf
    oplot,[1,1]*value,[0,10.*ymax],linestyle=2
    oplot,[1,1]*value-errmin,[0,10.*ymax],linestyle=1
    oplot,[1,1]*value+errmax,[0,10.*ymax],linestyle=1
    ndec=max([2,2-fix(alog10(value))])
    if ndec gt 4 then ndec=4
    xyouts,.215,.85,symstring+'!6 = '+textoidl(f_string(value,ndec)+'^{+'+f_string(errmax,ndec)+'}_{-'+f_string(errmin,ndec)+'}')+'!N '+unitstring,/normal,charsize=.8
    device,/close
    cumpdf=total(pdf*darray,/cumulative)
    ylabcum='!8p!6(<'+symstring+')'
    device,filename=figdir+galaxy+'_cumdistr_'+varstring+'.ps',xsize=12,ysize=9,/color,bits_per_pixel=8,/encapsulated
    plot,array,cumpdf,/nodata,/xlog,xtitle=xlab,ytitle=ylab,xr=arrange,xstyle=1
    oplot,array,cumpdf
    oplot,[1,1]*value,[0,1],linestyle=2
    oplot,[1,1]*value-errmin,[0,1],linestyle=1
    oplot,[1,1]*value+errmax,[0,1],linestyle=1
    ndec=max([2,2-fix(alog10(value))])
    if ndec gt 4 then ndec=4
    xyouts,.215,.85,symstring+'!6 = '+textoidl(f_string(value,ndec)+'^{+'+f_string(errmax,ndec)+'}_{-'+f_string(errmin,ndec)+'}')+'!N '+unitstring,/normal,charsize=.8
    device,/close
    set_plot,origdevice
    report='         1D Probability distribution functions of '+varstring+' plotted'
    return,report
end

function fitKL14,fluxratio_star,fluxratio_gas,err_star_log,err_gas_log,tstariso,beta_star,beta_gas,fstarover,fgasover,apertures_star,apertures_gas, $
                 surfcontrasts,surfcontrastg,peak_prof,tstar_incl,tgasmini,tgasmaxi,tovermini,tovermaxi,nfitstar,nfitgas,ndepth,ntry,galaxy,figdir,genplot,outputdir,arrdir,window_plot

    COMMON numbers

    if n_elements(apertures_star) ne n_elements(apertures_gas) then begin ;ensure that there are the same number of star and gas apertures
        print,' error: apertures_star and apertures_gas must have the same number of elements'
        print,' quitting...'
        stop
    endif
    napertures=n_elements(apertures_star)

    use=where(finite(alog10(fluxratio_star(0:napertures-1))) and finite(alog10(fluxratio_gas(0:napertures-1))))
    nvar=3
    nfit=total(nfitstar(use))+total(nfitgas(use))
    ndeg=nfit-nvar-1.
    weights=[nfitstar(use),nfitgas(use)]
    weights/=mean(weights)


    lambdamini=.3*min([apertures_star[use],apertures_gas[use]])
    lambdamaxi=3.*max([apertures_star[use],apertures_gas[use]])

    tgasmin=tgasmini
    tgasmax=tgasmaxi
    tovermin=tovermini
    tovermax=tovermaxi
    lambdamin=lambdamini
    lambdamax=lambdamaxi

    go=1
    n=0
    while go eq 1 and n lt ndepth do begin
        tgasarr=tgasmin*(tgasmax/tgasmin)^((dindgen(ntry)+.5)/ntry)
        toverarr=tovermin*(tovermax/tovermin)^((dindgen(ntry)+.5)/ntry)
        lambdaarr=lambdamin*(lambdamax/lambdamin)^((dindgen(ntry)+.5)/ntry)
        tgasarr3d=dblarr(ntry,ntry,ntry)
        toverarr3d=dblarr(ntry,ntry,ntry)
        lambdaarr3d=dblarr(ntry,ntry,ntry)
        dtgas=tgasarr*((tgasmax/tgasmin)^(1./(2.*ntry))-(tgasmax/tgasmin)^(-1./(2.*ntry)))
        dtover=toverarr*((tovermax/tovermin)^(1./(2.*ntry))-(tovermax/tovermin)^(-1./(2.*ntry)))
        dlambda=lambdaarr*((lambdamax/lambdamin)^(1./(2.*ntry))-(lambdamax/lambdamin)^(-1./(2.*ntry)))
        dtgasarr=dblarr(ntry,ntry,ntry)
        dtoverarr=dblarr(ntry,ntry,ntry)
        dlambdaarr=dblarr(ntry,ntry,ntry)
        redchi2=dblarr(ntry,ntry,ntry)
        chi2=dblarr(ntry,ntry,ntry)
        fluxratio_star_th=dblarr(ntry,ntry,ntry,napertures)
        fluxratio_gas_th=dblarr(ntry,ntry,ntry,napertures)
        timebefore=systime(1)
        for i=0,ntry-1 do begin
            tgasarr3d(i,*,*)=tgasarr(i)
            toverarr3d(*,i,*)=toverarr(i)
            lambdaarr3d(*,*,i)=lambdaarr(i)
            dtgasarr(i,*,*)=dtgas(i)
            dtoverarr(*,i,*)=dtover(i)
            dlambdaarr(*,*,i)=dlambda(i)
            for j=0,ntry-1 do begin
                if tstar_incl eq 0 then tstaruse=tstariso+toverarr(j) else tstaruse=tstariso
                beta_star_use=interpol(beta_star,fstarover,toverarr(j)/tstaruse)
                beta_gas_use=interpol(beta_gas,fgasover,toverarr(j)/tgasarr(i))
                for k=0,ntry-1 do begin
                    if toverarr(j) le min(tgasarr(i)) then begin
                        fluxratio_star_th(i,j,k,*)=f_fluxratiostar(tgasarr(i),tstaruse,toverarr(j),apertures_star,lambdaarr(k),beta_gas_use,surfcontrasts,surfcontrastg,peak_prof)
                        fluxratio_gas_th(i,j,k,*)=f_fluxratiogas(tgasarr(i),tstaruse,toverarr(j),apertures_gas,lambdaarr(k),beta_star_use,surfcontrasts,surfcontrastg,peak_prof)
                        diff=[(alog10(fluxratio_star(use)/fluxratio_star_th(i,j,k,use)))^2./err_star_log(use)^2., $
                              (alog10(fluxratio_gas(use)/fluxratio_gas_th(i,j,k,use)))^2./err_gas_log(use)^2.]
                    endif else diff=huge
                    chi2(i,j,k)=total(diff*weights)
                endfor
            endfor
            progress,'     ==> fitting progress '+f_string(n+1,0)+' (maximum '+f_string(ndepth,0)+')',i,ntry-1
        endfor
        redchi2=chi2/ndeg
        timeafter=systime(1)
        timeloop=timeafter-timebefore
        minchi2=min(redchi2)
        index=min(where(redchi2 eq minchi2))
        kbest=long(index/ntry/ntry)
        jbest=long((index-ntry*ntry*kbest)/ntry)
        ibest=long(index-ntry*jbest-ntry*ntry*kbest)
        prob=exp(-chi2/2.)
        probcst=1./total(prob*dtgasarr*dtoverarr*dlambdaarr)
        probnorm=prob*probcst

        probtgas=dblarr(ntry)
        probtover=dblarr(ntry)
        problambda=dblarr(ntry)
        probtgastover=dblarr(ntry,ntry)
        probtgaslambda=dblarr(ntry,ntry)
        probtoverlambda=dblarr(ntry,ntry)
        for i=0,ntry-1 do begin
            probtgas(i)=total(probnorm(i,*,*)*dtoverarr(i,*,*)*dlambdaarr(i,*,*))
            probtover(i)=total(probnorm(*,i,*)*dtgasarr(*,i,*)*dlambdaarr(*,i,*))
            problambda(i)=total(probnorm(*,*,i)*dtgasarr(*,*,i)*dtoverarr(*,*,i))
            for j=0,ntry-1 do begin
                probtgastover(i,j)=total(probnorm(i,j,*)*dlambdaarr(i,j,*))
                probtgaslambda(i,j)=total(probnorm(i,*,j)*dtoverarr(i,*,j))
                probtoverlambda(i,j)=total(probnorm(*,i,j)*dtgasarr(*,i,j))
            endfor
        endfor

        meantgas=total(tgasarr*probtgas*dtgas)
        meantover=total(toverarr*probtover*dtover)
        meanlambda=total(lambdaarr*problambda*dlambda)
        stdevtgas=sqrt(total((tgasarr-meantgas)^2.*probtgas*dtgas))
        stdevtover=sqrt(total((toverarr-meantover)^2.*probtover*dtover))
        stdevlambda=sqrt(total((lambdaarr-meanlambda)^2.*problambda*dlambda))

        meantgastover=total(tgasarr3d*toverarr3d*probnorm*dtgasarr*dtoverarr*dlambdaarr)
        meantgaslambda=total(tgasarr3d*lambdaarr3d*probnorm*dtgasarr*dtoverarr*dlambdaarr)
        meantoverlambda=total(toverarr3d*lambdaarr3d*probnorm*dtgasarr*dtoverarr*dlambdaarr)
        corrtgastover=(meantgastover-meantgas*meantover)/stdevtgas/stdevtover
        corrtgaslambda=(meantgaslambda-meantgas*meanlambda)/stdevtgas/stdevlambda
        corrtoverlambda=(meantoverlambda-meantover*meanlambda)/stdevtover/stdevlambda

        tgas=tgasarr(ibest)
        tover=toverarr(jbest)
        lambda=lambdaarr(kbest)
        if tstar_incl eq 0 then tstar=tstariso+tover else tstar=tstariso
        beta_star_best=interpol(beta_star,fstarover,tover/tstar)
        beta_gas_best=interpol(beta_gas,fgasover,tover/tgas)
        probtgascum=total(probtgas*dtgas,/cumulative)
        probtovercum=total(probtover*dtover,/cumulative)
        problambdacum=total(problambda*dlambda,/cumulative)
        dummy=f_pdftovalues(tgasarr,dtgas,probtgascum,tgas)
        tgas_errmin=dummy(1)
        tgas_errmax=dummy(2)
        if abs(tgas-dummy(0)) gt max([tgas_errmin,tgas_errmax]) then print,' WARNING: tgas PDF is highly asymmetric, recommend using PDF rather than calculated value and error bars'
        dummy=f_pdftovalues(toverarr,dtover,probtovercum,tover)
        tover_errmin=dummy(1)
        tover_errmax=dummy(2)
        if abs(tover-dummy(0)) gt max([tover_errmin,tover_errmax]) then print,' WARNING: tover PDF is highly asymmetric, recommend using PDF rather than calculated value and error bars'
        dummy=f_pdftovalues(lambdaarr,dlambda,problambdacum,lambda)
        lambda_errmin=dummy(1)
        lambda_errmax=dummy(2)
        if abs(lambda-dummy(0)) gt max([lambda_errmin,lambda_errmax]) then print,' WARNING: lambda PDF is highly asymmetric, recommend using PDF rather than calculated value and error bars'

        ;if best fit gets too close to the edge of the fitted range, stop iteratively shrinking the fitting range
        errminfrac=1.-.5*(1.+erf(1./sqrt(2.)))
        errmaxfrac=1.-errminfrac
        errminfractgas=2.*errminfrac*probtgascum(ibest) ;(.5-errmin)/.5=(1-errmin2/p) => errmin2 = 2*errmin*p (32nd percentile below best-fitting value)
        errmaxfractgas=2.*(1.-errmaxfrac)*probtgascum(ibest)+2.*(errmaxfrac-.5) ;(errmax-.5)/.5=(errmax2-p)/(1-p) => errmax2 = (2*errmax-1)*(1-p)+p = 2*(1-errmax)*p+2*(errmax-.5) (68th percentile above best-fitting value)
        imin=max(where(probtgascum le errminfractgas))
        if imin(0) eq -1 then imin=0
        imax=min(where(probtgascum ge errmaxfractgas))
        if imax(0) eq -1 then imax=ntry-1
        if imax(0) eq imin then imax=imin+3
        errminfractover=2.*errminfrac*probtovercum(jbest)
        errmaxfractover=2.*(1.-errmaxfrac)*probtovercum(jbest)+2.*(errmaxfrac-.5)
        jmin=max(where(probtovercum le errminfractover))
        if jmin(0) eq -1 then jmin=0
        jmax=min(where(probtovercum ge errmaxfractover))
        if jmax(0) eq -1 then jmax=ntry-1
        if jmax(0) eq jmin then jmax=jmin+3
        errminfraclambda=2.*errminfrac*problambdacum(kbest)
        errmaxfraclambda=2.*(1.-errmaxfrac)*problambdacum(kbest)+2.*(errmaxfrac-.5)
        kmin=max(where(problambdacum le errminfraclambda))
        if kmin(0) eq -1 then kmin=0
        kmax=min(where(problambdacum ge errmaxfraclambda))
        if kmax(0) eq -1 then kmax=ntry-1
        if kmax(0) eq kmin then kmax=kmin+3

        edge=4 ;stop refining if all fitting limits are within 4 elements from the array edge
        tgasminold=tgasmin
        tgasmaxold=tgasmax
        toverminold=tovermin
        tovermaxold=tovermax
        lambdaminold=lambdamin
        lambdamaxold=lambdamax
        if imin lt edge and imax gt ntry-edge-1 and jmin lt edge and jmax gt ntry-edge-1 and kmin lt edge and kmax gt ntry-edge-1 then go=0. else begin ;set fitting range for next loop
            iredchi=dblarr(ntry)
            jredchi=dblarr(ntry)
            kredchi=dblarr(ntry)
            for i=0,ntry-1 do iredchi(i)=min(redchi2(i,*,*))
            for j=0,ntry-1 do jredchi(j)=min(redchi2(*,j,*))
            for k=0,ntry-1 do kredchi(k)=min(redchi2(*,*,k))
            imin=max([0,min(where(iredchi le minchi2+3.))-1])
            imax=min([long(ntry-1),max(where(iredchi le minchi2+3.))+1])
            jmin=max([0,min(where(jredchi le minchi2+3.))-1])
            jmax=min([long(ntry-1),max(where(jredchi le minchi2+3.))+1])
            kmin=max([0,min(where(kredchi le minchi2+3.))-1])
            kmax=min([long(ntry-1),max(where(kredchi le minchi2+3.))+1])
            tgasmin=tgasarr(imin)-.5*dtgas(imin)
            tgasmax=tgasarr(imax)+.5*dtgas(imax)
            tovermin=min([toverarr(jmin)-.5*dtover(jmin),tgasmin])
            tovermax=toverarr(jmax)+.5*dtover(jmax)
            lambdamin=lambdaarr(kmin)-.5*dlambda(kmin)
            lambdamax=lambdaarr(kmax)+.5*dlambda(kmax)
        endelse
        ;avoid edge traps
        if ibest eq 0 then tgasmin=max([tgasmin^2./tgasmax,tgasmini])
        if ibest eq ntry-1 then tgasmax=min([tgasmax^2./tgasmin,tgasmaxi])
        if jbest eq 0 then tovermin=max([2.*tovermin-tovermax,tovermini])
        if jbest eq ntry-1 then tovermax=min([2.*tovermax-tovermin,tovermaxi])
        if kbest eq 0 then lambdamin=max([lambdamin^2./lambdamax,lambdamini])
        if kbest eq ntry-1 then lambdamax=min([lambdamax^2./lambdamin,lambdamaxi])
        ;if the change in the fitting range becomes smaller than 5% for all variables, then stop iteratively shrinking it
        tol=0.05
        if abs(tgasminold/tgasmin-1.) le tol and abs(tgasmaxold/tgasmax-1.) le tol and $
            abs(toverminold/tovermin-1.) le tol and abs(tovermaxold/tovermax-1.) le tol and $
            abs(lambdaminold/lambdamin-1.) le tol and abs(lambdamaxold/lambdamax-1.) le tol then go=0.
           
        logdifftgas=0.5*(2.-alog10(tgasmax/tgasmin)) ;is the range at least two orders of magnitude?
        if logdifftgas gt 0 then tgr=[tgasmin/10.^logdifftgas,tgasmax*10.^logdifftgas] else tgr=[tgasmin,tgasmax] ;if not, extend for plotting
        logdifftover=0.5*(2.-alog10(tovermax/tovermin)) ;is the range at least two orders of magnitude?
        if logdifftover gt 0 then tor=[tovermin/10.^logdifftover,tovermax*10.^logdifftover] else tor=[tovermin,tovermax] ;if not, extend for plotting
        logdifflambda=0.5*(2.-alog10(lambdamax/lambdamin)) ;is the range at least two orders of magnitude?
        if logdifflambda gt 0 then lar=[lambdamin/10.^logdifflambda,lambdamax*10.^logdifflambda] else lar=[lambdamin,lambdamax] ;if not, extend for plotting

        xmin=min([apertures_star(0),apertures_gas(0)])^2./min([apertures_star(1),apertures_gas(1)])
        xmax=max([apertures_star(napertures-1),apertures_gas(napertures-1)])^2./max([apertures_star(napertures-2),apertures_gas(napertures-2)])
        ymin=min([fluxratio_star(use),1./fluxratio_gas(use)])
        ymax=max([fluxratio_gas(use),1./fluxratio_star(use)])
        set_plot,window_plot
        nr=1001
        r=min([apertures_star,apertures_gas])*10.^(dindgen(nr)/(nr-1)*alog10(max([apertures_star,apertures_gas])/min([apertures_star,apertures_gas])))
        fg=f_fluxratiogas(tgasarr(ibest),tstar,toverarr(jbest),r,lambdaarr(kbest),beta_star_best,surfcontrasts,surfcontrastg,peak_prof)
        fs=f_fluxratiostar(tgasarr(ibest),tstar,toverarr(jbest),r,lambdaarr(kbest),beta_gas_best,surfcontrasts,surfcontrastg,peak_prof)
        plot,apertures_star,fluxratio_star_th(ibest,jbest,kbest,*),/xlog,/ylog,xr=[xmin,xmax],yr=[ymin,ymax],/nodata,xstyle=1
        oplot,r,fg,thick=5
        oplot,r,fs,thick=5
        for i=0,ntry-1,fix(ntry/2.) do begin
            beta_gas_use=interpol(beta_gas,fgasover,toverarr(jbest)/tgasarr(i))
            fg=f_fluxratiogas(tgasarr(i),tstar,toverarr(jbest),r,lambdaarr(kbest),beta_star_best,surfcontrasts,surfcontrastg,peak_prof)
            fs=f_fluxratiostar(tgasarr(i),tstar,toverarr(jbest),r,lambdaarr(kbest),beta_gas_use,surfcontrasts,surfcontrastg,peak_prof)
            oplot,r,fg,linestyle=2
            oplot,r,fs,linestyle=2
            beta_star_use=interpol(beta_star,fstarover,toverarr(i)/tstar)
            beta_gas_use=interpol(beta_gas,fgasover,toverarr(i)/tgasarr(ibest))
            fg=f_fluxratiogas(tgasarr(ibest),tstar,toverarr(i),r,lambdaarr(kbest),beta_star_use,surfcontrasts,surfcontrastg,peak_prof)
            fs=f_fluxratiostar(tgasarr(ibest),tstar,toverarr(i),r,lambdaarr(kbest),beta_gas_use,surfcontrasts,surfcontrastg,peak_prof)
            oplot,r,fg,linestyle=1
            oplot,r,fs,linestyle=1
            fg=f_fluxratiogas(tgasarr(ibest),tstar,toverarr(jbest),r,lambdaarr(i),beta_star_best,surfcontrasts,surfcontrastg,peak_prof)
            fs=f_fluxratiostar(tgasarr(ibest),tstar,toverarr(jbest),r,lambdaarr(i),beta_gas_best,surfcontrasts,surfcontrastg,peak_prof)
            oplot,r,fg,linestyle=3
            oplot,r,fs,linestyle=3
        endfor
        oplot,apertures_star,fluxratio_star,psym=4
        oplot,apertures_gas,fluxratio_gas,psym=4
        oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)+err_star_log)-fluxratio_star,psym=3,/hibar
        oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)-err_star_log)-fluxratio_star,psym=3,/lobar
        oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)+err_gas_log)-fluxratio_gas,psym=3,/hibar
        oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)-err_gas_log)-fluxratio_gas,psym=3,/lobar
        set_plot,'x'
        if timeloop lt 1. then wait,1.

        if genplot eq 1 then begin
            a=.7
            x=a*cos(findgen(33)*!pi/16.)
            y=a*sin(findgen(33)*!pi/16.)
            USERSYM,x,y,/FILL
            
            set_plot,'ps'
            xmin=0
            xmax=1
            ymin=0
            ymax=3
            device,filename=figdir+galaxy+'_beta_star.ps',xsize=12,ysize=9,/color,bits_per_pixel=8,/encapsulated
            plot,fstarover,beta_star,xr=[xmin,xmax],yr=[ymin,ymax],/nodata,xstyle=1,ystyle=1,xtitle='!8f!6!Dstar,over!N=!8t!6!Dover!N/!8t!6!Dstar!N',ytitle='!7b!6!Dstar!N',charsize=1.3
            oplot,fstarover,beta_star
            oplot,[1,1]*tover/tstar,[ymin,ymax],linestyle=1
            oplot,[xmin,xmax],[1,1]*beta_star_best,linestyle=1
            oplot,[tover/tstar],[beta_star_best],psym=8
            device,/close
            
            device,filename=figdir+galaxy+'_beta_gas.ps',xsize=12,ysize=9,/color,bits_per_pixel=8,/encapsulated
            plot,fgasover,beta_gas,xr=[0,1],yr=[0,3],/nodata,xstyle=1,ystyle=1,xtitle='!8f!6!Dgas,over!N=!8t!6!Dover!N/!8t!6!Dgas!N',ytitle='!7b!6!Dgas!N',charsize=1.3
            oplot,fgasover,beta_gas
            oplot,[1,1]*tover/tgas,[ymin,ymax],linestyle=1
            oplot,[xmin,xmax],[1,1]*beta_gas_best,linestyle=1
            oplot,[tover/tgas],[beta_gas_best],psym=8
            device,/close

            nr=1001
            r=min([apertures_star,apertures_gas])*10.^(dindgen(nr)/(nr-1)*alog10(max([apertures_star,apertures_gas])/min([apertures_star,apertures_gas])))
            fg=f_fluxratiogas(tgas,tstar,tover,r,lambda,beta_star_best,surfcontrasts,surfcontrastg,peak_prof)
            fs=f_fluxratiostar(tgas,tstar,tover,r,lambda,beta_gas_best,surfcontrasts,surfcontrastg,peak_prof)
            device,filename=figdir+galaxy+'_fit.ps',xsize=12,ysize=11,/color,bits_per_pixel=8,/encapsulated
            xmin=10.*round(0.5*min(r)/10.)
            xmax=10.*round(2.*max(r)/10.)
            decmin=alog10(min([fs,10.^(alog10(fluxratio_star)-err_star_log)])/2.)
            decminfix=fix(decmin)-1.
            ymin=10.^decminfix*fix(10.^(decmin-decminfix))
            decmax=alog10(max([fg,10.^(alog10(fluxratio_gas)+err_gas_log)])*2.)
            decmaxfix=fix(decmax)
            ymax=10.^decmaxfix*(fix(10.^(decmax-decmaxfix))+1.)
            plot,r,fg,/xlog,/ylog,yr=[ymin,ymax],/nodata,xr=[xmin,xmax],xstyle=1,ystyle=1,xtitle='!8l!6!Dap!N [pc]',ytitle='!8t!6!Ddepl,peaks!N/!8t!6!Ddepl,tot!N',charsize=1.3,xtickformat='logticks',ytickformat='logticks'
            oplot,[xmin,xmax],[1,1]*fg(0),linestyle=2,color=fsc_color('black'),thick=3
            oplot,[xmin,xmax],[1,1]*fs(0),linestyle=2,color=fsc_color('black'),thick=3
            xdyn=alog10(xmax/xmin)
            ydyn=alog10(ymax/ymin)
            err_star_log_grey=err_star_log/sqrt(nfitstar)*mean(nfitstar)*sqrt(total(nfitgas+nfitstar)/2./n_elements(use))
            err_gas_log_grey=err_gas_log/sqrt(nfitgas)*mean(nfitgas)*sqrt(total(nfitgas+nfitstar)/2./n_elements(use))
            oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)+err_star_log_grey)-fluxratio_star,psym=3,/hibar,errcolor=180,errthick=20,/nohat
            oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)-err_star_log_grey)-fluxratio_star,psym=3,/lobar,errcolor=180,errthick=20,/nohat
            oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)+err_gas_log_grey)-fluxratio_gas,psym=3,/hibar,errcolor=180,errthick=20,/nohat
            oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)-err_gas_log_grey)-fluxratio_gas,psym=3,/lobar,errcolor=180,errthick=20,/nohat
            arrowlen=.1*ydyn
            arrowoff=0.25*arrowlen
            oplot,[1,1]*lambda,ymin*[10.^(arrowoff),10.^(arrowoff+arrowlen)],thick=3
            oplot,[1,10.^(0.7*arrowoff*xdyn/ydyn)]*lambda,ymin*[10.^(arrowoff),10.^(2.*arrowoff)],thick=3
            oplot,[1,10.^(-0.7*arrowoff*xdyn/ydyn)]*lambda,ymin*[10.^(arrowoff),10.^(2.*arrowoff)],thick=3
            oplot,r,fs,thick=5,color=fsc_color('green')
            oplot,r,fg,thick=5,color=fsc_color('green')
            oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)+err_star_log)-fluxratio_star,psym=3,/hibar,errcolor=fsc_color('blue'),errthick=3
            oploterror,apertures_star,fluxratio_star,10.^(alog10(fluxratio_star)-err_star_log)-fluxratio_star,psym=3,/lobar,errcolor=fsc_color('blue'),errthick=3
            oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)+err_gas_log)-fluxratio_gas,psym=3,/hibar,errcolor=fsc_color('red'),errthick=3
            oploterror,apertures_gas,fluxratio_gas,10.^(alog10(fluxratio_gas)-err_gas_log)-fluxratio_gas,psym=3,/lobar,errcolor=fsc_color('red'),errthick=3
            oplot,apertures_star,fluxratio_star,psym=8,color=fsc_color('blue'),thick=2
            oplot,apertures_gas,fluxratio_gas,psym=8,color=fsc_color('red'),thick=2
            xyouts,.805,.863,'!6Gas',color=fsc_color('red'),charthick=3,charsize=1.3,/normal
            xyouts,.805,.803,'!6Stars',color=fsc_color('blue'),charthick=3,charsize=1.3,/normal
            xyouts,.675-(strlen(galaxy)-3)*.023,.8204,'!6'+strmid(galaxy,0,strlen(galaxy))+':',charthick=3,charsize=1.3,/normal
            device,/close

            fracsig1=(1.+erf(1./sqrt(2.)))-1.
            fracsig2=(1.+erf(2./sqrt(2.)))-1.
            fracsig3=(1.+erf(3./sqrt(2.)))-1.

            maxprob=max(probtgastover)
            nc=1001
            contours=dindgen(nc)/(nc-1.)*maxprob
            fracarr=dblarr(nc)
            for i=0,nc-1 do begin
                percentile=where(probtgastover ge contours(i))
                dtgasarr2d=reform(dtgasarr(*,*,0))
                dtoverarr2d=reform(dtoverarr(*,*,0))
                fracarr(i)=total(probtgastover(percentile)*dtgasarr2d(percentile)*dtoverarr2d(percentile))
            endfor
            isig1=max(where(fracarr ge fracsig1))
            if isig1(0) eq -1 then isig1=0
            if isig1 eq nc-1 then isig1=isig1-1
            frac=(fracsig1-fracarr(isig1))/(fracarr(isig1+1)-fracarr(isig1))
            contoursig1=contours(isig1)+frac*(contours(isig1+1)-contours(isig1))
            isig2=max(where(fracarr ge fracsig2))
            if isig2(0) eq -1 then isig2=0
            if isig2 eq nc-1 then isig2=isig2-1
            frac=(fracsig2-fracarr(isig2))/(fracarr(isig2+1)-fracarr(isig2))
            contoursig2=contours(isig2)+frac*(contours(isig2+1)-contours(isig2))
            isig3=max(where(fracarr ge fracsig3))
            if isig3(0) eq -1 then isig3=0
            if isig3 eq nc-1 then isig3=isig3-1
            frac=(fracsig3-fracarr(isig3))/(fracarr(isig3+1)-fracarr(isig3))
            contoursig3=contours(isig3)+frac*(contours(isig3+1)-contours(isig3))
            contourlvl=[contoursig3,contoursig2,contoursig1,maxprob]
            finelvl=(dindgen(256)/255.)^3.*maxprob
            show=[2,3]
            inter=[0,1]
            labels=['3!7r!6','2!7r!6','1!7r!6',textoidl('\chi^2_{min}')]
            device,filename=figdir+galaxy+'_prob_tgas_tover.ps',xsize=11,ysize=11,/color,bits_per_pixel=8,/encapsulated
            contour,probtgastover,tgasarr,toverarr,/nodata,/xlog,/ylog,xtitle='!8t!6!Dgas!N [Myr]',ytitle='!8t!6!Dover!N [Myr]',xr=tgr,yr=tor,xstyle=5,ystyle=5,charsize=chsize,position=aspect(1.)
            ctload,22
            contour,probtgastover,tgasarr,toverarr,levels=finelvl,/fill,/overplot,charsize=chsize
            ctload,0
            contour,probtgastover,tgasarr,toverarr,/nodata,/xlog,/ylog,/noerase,xtitle='!8t!6!Dgas!N [Myr]',ytitle='!8t!6!Dover!N [Myr]',xr=tgr,yr=tor,xstyle=1,ystyle=1,charsize=chsize,position=aspect(1.)
            contour,probtgastover,tgasarr,toverarr,levels=contourlvl(show),c_annotation=labels(show),/overplot,c_linestyle=2,charsize=chsize
            contour,probtgastover,tgasarr,toverarr,levels=contourlvl(inter),c_annotation=labels(inter),/overplot,c_linestyle=1,charsize=chsize
            oplot,[tgasarr(ibest)],[toverarr(jbest)],psym=1
            oplot,[tgasarr(ibest)-tgas_errmin,tgasarr(ibest)+tgas_errmax],[1,1]*toverarr(jbest),linestyle=0
            oplot,[1,1]*tgasarr(ibest),[toverarr(jbest)-tover_errmin,toverarr(jbest)+tover_errmax],linestyle=0
            xyouts,.2,.77,'!7q!6 = '+f_string(corrtgastover,2),/normal,charsize=chsize
            device,/close

            maxprob=max(probtgaslambda)
            nc=1001
            contours=dindgen(nc)/(nc-1.)*maxprob
            fracarr=dblarr(nc)
            for i=0,nc-1 do begin
                percentile=where(probtgaslambda ge contours(i))
                dtgasarr2d=reform(dtgasarr(*,0,*))
                dlambdaarr2d=reform(dlambdaarr(*,0,*))
                fracarr(i)=total(probtgaslambda(percentile)*dtgasarr2d(percentile)*dlambdaarr2d(percentile))
            endfor
            isig1=max(where(fracarr ge fracsig1))
            if isig1(0) eq -1 then isig1=0
            if isig1 eq nc-1 then isig1=isig1-1
            frac=(fracsig1-fracarr(isig1))/(fracarr(isig1+1)-fracarr(isig1))
            contoursig1=contours(isig1)+frac*(contours(isig1+1)-contours(isig1))
            isig2=max(where(fracarr ge fracsig2))
            if isig2(0) eq -1 then isig2=0
            if isig2 eq nc-1 then isig2=isig2-1
            frac=(fracsig2-fracarr(isig2))/(fracarr(isig2+1)-fracarr(isig2))
            contoursig2=contours(isig2)+frac*(contours(isig2+1)-contours(isig2))
            isig3=max(where(fracarr ge fracsig3))
            if isig3(0) eq -1 then isig3=0
            if isig3 eq nc-1 then isig3=isig3-1
            frac=(fracsig3-fracarr(isig3))/(fracarr(isig3+1)-fracarr(isig3))
            contoursig3=contours(isig3)+frac*(contours(isig3+1)-contours(isig3))
            contourlvl=[contoursig3,contoursig2,contoursig1,maxprob]
            finelvl=(dindgen(256)/255.)^3.*maxprob
            show=[2,3]
            inter=[0,1]
            labels=['3!7r!6','2!7r!6','1!7r!6',textoidl('\chi^2_{min}')]
            device,filename=figdir+galaxy+'_prob_tgas_lambda.ps',xsize=11,ysize=11,/color,bits_per_pixel=8,/encapsulated
            contour,probtgaslambda,tgasarr,lambdaarr,/nodata,/xlog,/ylog,xtitle='!8t!6!Dgas!N [Myr]',ytitle='!7k!6 [pc]',xr=tgr,yr=lar,xstyle=5,ystyle=5,charsize=chsize,position=aspect(1.)
            ctload,22
            contour,probtgaslambda,tgasarr,lambdaarr,levels=finelvl,/fill,/overplot,charsize=chsize
            ctload,0
            contour,probtgaslambda,tgasarr,lambdaarr,/nodata,/xlog,/ylog,/noerase,xtitle='!8t!6!Dgas!N [Myr]',ytitle='!7k!6 [pc]',xr=tgr,yr=lar,xstyle=1,ystyle=1,charsize=chsize,position=aspect(1.)
            contour,probtgaslambda,tgasarr,lambdaarr,levels=contourlvl(show),c_annotation=labels(show),/overplot,c_linestyle=2,charsize=chsize
            contour,probtgaslambda,tgasarr,lambdaarr,levels=contourlvl(inter),c_annotation=labels(inter),/overplot,c_linestyle=1,charsize=chsize
            oplot,[tgasarr(ibest)],[lambdaarr(kbest)],psym=1
            oplot,[tgasarr(ibest)-tgas_errmin,tgasarr(ibest)+tgas_errmax],[1,1]*lambdaarr(kbest),linestyle=0
            oplot,[1,1]*tgasarr(ibest),[lambdaarr(kbest)-lambda_errmin,lambdaarr(kbest)+lambda_errmax],linestyle=0
            xyouts,.2,.77,'!7q!6 = '+f_string(corrtgaslambda,2),/normal,charsize=chsize
            device,/close

            maxprob=max(probtoverlambda)
            nc=1001
            contours=dindgen(nc)/(nc-1.)*maxprob
            fracarr=dblarr(nc)
            for i=0,nc-1 do begin
                percentile=where(probtoverlambda ge contours(i))
                dtoverarr2d=reform(dtoverarr(0,*,*))
                dlambdaarr2d=reform(dlambdaarr(0,*,*))
                fracarr(i)=total(probtoverlambda(percentile)*dtoverarr2d(percentile)*dlambdaarr2d(percentile))
            endfor
            isig1=max(where(fracarr ge fracsig1))
            if isig1(0) eq -1 then isig1=0
            if isig1 eq nc-1 then isig1=isig1-1
            frac=(fracsig1-fracarr(isig1))/(fracarr(isig1+1)-fracarr(isig1))
            contoursig1=contours(isig1)+frac*(contours(isig1+1)-contours(isig1))
            isig2=max(where(fracarr ge fracsig2))
            if isig2(0) eq -1 then isig2=0
            if isig2 eq nc-1 then isig2=isig2-1
            frac=(fracsig2-fracarr(isig2))/(fracarr(isig2+1)-fracarr(isig2))
            contoursig2=contours(isig2)+frac*(contours(isig2+1)-contours(isig2))
            isig3=max(where(fracarr ge fracsig3))
            if isig3(0) eq -1 then isig3=0
            if isig3 eq nc-1 then isig3=isig3-1
            frac=(fracsig3-fracarr(isig3))/(fracarr(isig3+1)-fracarr(isig3))
            contoursig3=contours(isig3)+frac*(contours(isig3+1)-contours(isig3))
            contourlvl=[contoursig3,contoursig2,contoursig1,maxprob]
            finelvl=(dindgen(256)/255.)^3.*maxprob
            show=[2,3]
            inter=[0,1]
            labels=['3!7r!6','2!7r!6','1!7r!6',textoidl('\chi^2_{min}')]
            device,filename=figdir+galaxy+'_prob_tover_lambda.ps',xsize=11,ysize=11,/color,bits_per_pixel=8,/encapsulated
            contour,probtoverlambda,toverarr,lambdaarr,/nodata,/xlog,/ylog,xtitle='!8t!6!Dover!N [Myr]',ytitle='!7k!6 [pc]',xr=tor,yr=lar,xstyle=5,ystyle=5,charsize=chsize,position=aspect(1.)
            ctload,22
            contour,probtoverlambda,toverarr,lambdaarr,levels=finelvl,/fill,/overplot,charsize=chsize
            ctload,0
            contour,probtoverlambda,toverarr,lambdaarr,/nodata,/xlog,/ylog,/noerase,xtitle='!8t!6!Dover!N [Myr]',ytitle='!7k!6 [pc]',xr=tor,yr=lar,xstyle=1,ystyle=1,charsize=chsize,position=aspect(1.)
            contour,probtoverlambda,toverarr,lambdaarr,levels=contourlvl(show),c_annotation=labels(show),/overplot,c_linestyle=2,charsize=chsize
            contour,probtoverlambda,toverarr,lambdaarr,levels=contourlvl(inter),c_annotation=labels(inter),/overplot,c_linestyle=1,charsize=chsize
            oplot,[toverarr(jbest)],[lambdaarr(kbest)],psym=1
            oplot,[toverarr(jbest)-tover_errmin,toverarr(jbest)+tover_errmax],[1,1]*lambdaarr(kbest),linestyle=0
            oplot,[1,1]*toverarr(jbest),[lambdaarr(kbest)-lambda_errmin,lambdaarr(kbest)+lambda_errmax],linestyle=0
            xyouts,.2,.77,'!7q!6 = '+f_string(corrtoverlambda,2),/normal,charsize=chsize
            device,/close

            report=f_plotdistr(tgasarr,dtgas,probtgas,tgas,tgas_errmin,tgas_errmax,galaxy,figdir,'tgas','!8t!6!Dgas!N','!6Myr',0)
            report=f_plotdistr(toverarr,dtover,probtover,tover,tover_errmin,tover_errmax,galaxy,figdir,'tover','!8t!6!Dover!N','!6Myr',0)
            report=f_plotdistr(lambdaarr,dlambda,problambda,lambda,lambda_errmin,lambda_errmax,galaxy,figdir,'lambda','!7k!6','!6pc',0)

            set_plot,'x'
        endif

        n=n+1
    endwhile
        
    print,'     ==> writing tuningfork plot data to output directory'
    report=f_writetf_model(alog10(r),alog10(fs),alog10(fg),galaxy,outputdir,'# log10(apertures), log10(fluxratio_star_model), log10(fluxratio_gas_model)')
    report=f_writetf_obs(alog10(apertures_star),alog10(apertures_gas),alog10(fluxratio_star),alog10(fluxratio_gas),err_star_log,err_gas_log,err_star_log_grey,err_gas_log_grey,galaxy,outputdir, $
                         '# log10(apertures_star), log10(apertures_gas), log10(fluxratio_star_obs), log10(fluxratio_gas_obs), ' $
                         +'log10(fluxratio_star_errorbar), log10(fluxratio_gas_errorbar), log10(fluxratio_star_errorgrey), log10(fluxratio_gas_errorgrey)')

    print,'     ==> writing probability distribution functions to output directory'
    report=f_writepdf(alog10(tgasarr),alog10(dtgas),alog10(probtgas),galaxy,outputdir,'tgas','# log10(tgas[Myr]), log10(dtgas[Myr]), log10(PDF[Myr^-1])')
    report=f_writepdf(alog10(toverarr),alog10(dtover),alog10(probtover),galaxy,outputdir,'tover','# log10(tover[Myr]), log10(dtover[Myr]), log10(PDF[Myr^-1])')
    report=f_writepdf(alog10(lambdaarr),alog10(dlambda),alog10(problambda),galaxy,outputdir,'lambda','# log10(lambda[pc]), log10(dlambda[pc]), log10(PDF[pc^-1])')

    print,'     ==> writing probability distribution arrays to temporary array directory'
    save,filename=arrdir+'probnorm.sav',probnorm
    save,filename=arrdir+'probtgastover.sav',probtgastover
    save,filename=arrdir+'probtgaslambda.sav',probtgaslambda
    save,filename=arrdir+'probtoverlambda.sav',probtoverlambda
    save,filename=arrdir+'probtgas.sav',probtgas
    save,filename=arrdir+'probtover.sav',probtover
    save,filename=arrdir+'problambda.sav',problambda
    save,filename=arrdir+'tgasarr.sav',tgasarr
    save,filename=arrdir+'toverarr.sav',toverarr
    save,filename=arrdir+'lambdaarr.sav',lambdaarr
    save,filename=arrdir+'dtgas.sav',dtgas
    save,filename=arrdir+'dtover.sav',dtover
    save,filename=arrdir+'dtlambda.sav',dlambda
    save,filename=arrdir+'tgasarr3d.sav',tgasarr3d
    save,filename=arrdir+'toverarr3d.sav',toverarr3d
    save,filename=arrdir+'lambdaarr3d.sav',lambdaarr3d
    save,filename=arrdir+'dtgasarr.sav',dtgasarr
    save,filename=arrdir+'dtoverarr.sav',dtoverarr
    save,filename=arrdir+'dtlambdaarr.sav',dlambdaarr

    return,[minchi2,tgas,tgas_errmin,tgas_errmax,tover,tover_errmin,tover_errmax,lambda,lambda_errmin,lambda_errmax,beta_star_best,beta_gas_best]

end
