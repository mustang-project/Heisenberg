;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DERIVE PHYSICAL QUANTITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function f_zeta,surfcontrast,ttot,tref,peak_prof ;Tracer peak concentration parameter
    COMMON numbers
    if peak_prof eq 0 then zeta=tiny
    if peak_prof eq 1 then zeta=sqrt(ttot/(surfcontrast*tref))
    if peak_prof eq 2 then zeta=1./sqrt(2.)*sqrt(ttot/(2.*alog(2.)*surfcontrast*tref))
    return,zeta
end

function f_rpeak,zeta,lambda ;Tracer peak radius
    rpeak=0.5*zeta*lambda
    return,rpeak
end

function f_esf,tgas,tdepl,fcl,fgmc ;Star formation efficiency per star formation event
    esf=min([fcl*tgas*1.d6/(fgmc*tdepl),1.]) ;cannot exceed unity
    return,esf
end

function f_mdotsf,tgas,lambda,surfgas,fgmc,esf ;Star formation rate per star formation event during the gas phase
    mdotsf=!pi*lambda^2.*fgmc*esf*surfgas/(4.*tgas*1.d6) ;in Msun/yr
    return,mdotsf
end

function f_mdotfb,tover,lambda,surfgas,fgmc,esf ;Mass outflow rate per star formation event during the overlap phase
    COMMON numbers
    mdotfb=max([!pi*lambda^2.*fgmc*(1.-esf)*surfgas/(4.*tover*1.d6),tiny]) ;in Msun/yr
    return,mdotfb
end

function f_vfb,tover,lambda ;Feedback-driven expansion velocity of ejecta
    COMMON astrconst
    vfb=0.5*(lambda*pc)/(tover*myr)/kms ;v in km/s
    return,vfb
end

function f_etainst,tgas,tover,esf ;Instantaneous mass loading factor
    COMMON numbers
    etainst=max([(1.-esf)*tgas/(esf*tover),tiny])
    return,etainst
end

function f_etaavg,esf ;Time-integrated mass loading factor
    COMMON numbers
    etaavg=max([(1.-esf)/esf,tiny])
    return,etaavg
end

function f_chie,tover,esf,vfb,psi ;Feedback energy efficiency
    COMMON numbers
    COMMON astrconst
    chie=max([(1.-esf)*(vfb*kms)^2./(2.*esf*(tover*myr)*psi),tiny])
    return,chie
end

function f_chip,tover,surfgas,esf,vfb,psi,phitrap,kappa0 ;Feedback momentum efficiency
    COMMON numbers
    COMMON physconst
    COMMON astrconst
    taueff=kappa0^2.*esf*(surfgas*msun/pc^2.)^2.*psi/(4.*sigmaboltz) ;J s^-1 m^-2 K^-4
    chip=max([(1.-esf)*(vfb*kms)*clight/((1.+taueff)*esf*(tover*myr)*psi),tiny])
    return,chip
end

function f_getdvar,var ;Note: input must be positive and real
    COMMON numbers
    nvar=n_elements(var)
    dvar=dblarr(nvar)
    for i=0,nvar-1 do begin
        if i eq 0 then dvar(i)=min([2.*var(i),var(i+1)-var(i)]) ;make sure that var-0.5*dvar cannot be negative
        if i eq nvar-1 then dvar(i)=var(i)-var(i-1)
        if i gt 0 and i lt nvar-1 then dvar(i)=.5*(var(i+1)-var(i-1))
        if dvar(i) eq 0. then dvar(i)=tiny
    endfor
    return,dvar
end

function f_createpdf,array,darray,cumpdf,nnew
    COMMON numbers
    norig=n_elements(array)
    minval=min(array)-.5*darray(0) ;assume linearly-symmetric bins as input
    maxval=max(array)+.5*darray(norig-1) ;assume linearly-symmetric bins as input
    pdf=dblarr(nnew)
    if minval le 0. then logrange=0 else logrange=1
    if minval ne maxval then begin
        if logrange then begin
            newarray=minval*(maxval/minval)^((dindgen(nnew)+.5)/nnew)
            left=newarray*(maxval/minval)^(-1./(2.*nnew)) ;left bin edge
            right=newarray*(maxval/minval)^(1./(2.*nnew)) ;right bin edge
        endif else begin
            newarray=minval+(maxval-minval)*((dindgen(nnew)+.5)/nnew)
            left=newarray-(maxval-minval)/(2.*nnew) ;left bin edge
            right=newarray+(maxval-minval)/(2.*nnew) ;right bin edge
        endelse
        newdarray=right-left
        for i=0,nnew-1 do begin
            if left(i) gt min(array) then cumpdfleft=interpol(cumpdf,array,left(i)) else begin ;left cumpdf value
                interpolmin=max([0,min(where(array gt min(array)))]) ;avoid interpolation over saturated maximum values (crashes interpol function) and just include number once 
                interpolarr=[0,findgen(n_elements(cumpdf)-interpolmin)+interpolmin] ;index array for interpolation
                cumpdfleft=max([0.,interpol(cumpdf(interpolarr),array(interpolarr),left(i))]) ;left cumpdf value
            endelse
            if right(i) lt max(array) then cumpdfright=interpol(cumpdf,array,right(i)) else begin ;right cumpdf value
                interpolmax=min([n_elements(cumpdf)-1,max(where(array lt max(array)))]) ;avoid interpolation over saturated maximum values (crashes interpol function) and just include number once 
                if interpolmax lt n_elements(cumpdf)-1 then interpolarr=[findgen(interpolmax+1),n_elements(cumpdf)-1] else interpolarr=findgen(interpolmax+1) ;index array for interpolation
                cumpdfright=min([interpol(cumpdf(interpolarr),array(interpolarr),right(i)),1.]) ;right cumpdf value
            endelse
            pdf(i)=(cumpdfright-cumpdfleft)/newdarray(i) ;pdf value in bin
        endfor
    endif else begin
        newarray=dblarr(nnew)
        newdarray=dblarr(nnew)+tiny
        for i=0,nnew-1 do newarray(i)=array(0)+(i-nnew)/2*tiny
        pdf(*)=(nnew*tiny)^(-1.)
    endelse
    if finite(total(pdf)) eq 0 then stop
    return,[[newarray],[newdarray],[pdf]]
end

function derivephys,surfsfr,surfsfr_err,surfgas,surfgas_err,area,tgas,tover,lambda,fcl,fgmc,tstariso,tstariso_rerrmin,tstariso_rerrmax, $
                    tstar_incl,surfcontrasts,surfcontrastg,psi,phitrap,kappa0,peak_prof,ntry,nphysmc,galaxy,outputdir,arrdir,figdir
        
    ;derived quantities for best-fitting model
    tdepl=surfgas/surfsfr
    if tstar_incl eq 0 then tstar=tstariso+tover else tstar=tstariso
    ttotal=tgas+tstar-tover
    zetastar=f_zeta(surfcontrasts,ttotal,tstar,peak_prof)
    zetagas=f_zeta(surfcontrastg,ttotal,tgas,peak_prof)
    rpeakstar=f_rpeak(zetastar,lambda)
    rpeakgas=f_rpeak(zetagas,lambda)
    esf=f_esf(tgas,tdepl,fcl,fgmc)
    mdotsf=f_mdotsf(tgas,lambda,surfgas,fgmc,esf)
    mdotfb=f_mdotfb(tover,lambda,surfgas,fgmc,esf)
    vfb=f_vfb(tover,lambda)
    etainst=f_etainst(tgas,tover,esf)
    etaavg=f_etaavg(esf)
    chie=f_chie(tover,esf,vfb,psi)
    chip=f_chip(tover,surfgas,esf,vfb,psi,phitrap,kappa0)
    
    restore,filename=arrdir+'probnorm.sav'
    restore,filename=arrdir+'probtgas.sav'
    restore,filename=arrdir+'probtover.sav'
    restore,filename=arrdir+'problambda.sav'
    restore,filename=arrdir+'tgasarr.sav'
    restore,filename=arrdir+'toverarr.sav'
    restore,filename=arrdir+'lambdaarr.sav'
    restore,filename=arrdir+'dtgas.sav'
    restore,filename=arrdir+'dtover.sav'
    restore,filename=arrdir+'dtlambda.sav'
    restore,filename=arrdir+'tgasarr3d.sav'
    restore,filename=arrdir+'toverarr3d.sav'
    restore,filename=arrdir+'lambdaarr3d.sav'
    restore,filename=arrdir+'dtgasarr.sav'
    restore,filename=arrdir+'dtoverarr.sav'
    restore,filename=arrdir+'dtlambdaarr.sav'
    
    probtgas3d=probnorm
    probtover3d=transpose(probnorm,[1,2,0])
    problambda3d=transpose(probnorm,[2,0,1])
    tgasarr3d=tgasarr3d
    toverarr3d=transpose(toverarr3d,[1,2,0])
    lambdaarr3d=transpose(lambdaarr3d,[2,0,1])
    
    ;get cumulative PDFs for use as generation functions
    probtgascum3d=total(probnorm*dtgasarr*dtoverarr*dlambdaarr,/cumulative)
    probtovercum3d=total(probtover3d*transpose(dtgasarr,[1,2,0])*transpose(dtoverarr,[1,2,0])*transpose(dlambdaarr,[1,2,0]),/cumulative)
    problambdacum3d=total(problambda3d*transpose(dtgasarr,[2,0,1])*transpose(dtoverarr,[2,0,1])*transpose(dlambdaarr,[2,0,1]),/cumulative)
    probtgascum=total(probtgas*dtgas,/cumulative)
    probtovercum=total(probtover*dtover,/cumulative)
    problambdacum=total(problambda*dlambda,/cumulative)
    
    ;create Monte Carlo arrays
    tgasmc=dblarr(nphysmc)
    tovermc=dblarr(nphysmc)
    lambdamc=dblarr(nphysmc)
    surfsfrmc=dblarr(nphysmc)
    surfgasmc=dblarr(nphysmc)
    tstarisomc=dblarr(nphysmc)
    
    tdeplmc=dblarr(nphysmc)+tdepl
    tstarmc=dblarr(nphysmc)+tstar
    ttotalmc=dblarr(nphysmc)+ttotal
    zetastarmc=dblarr(nphysmc)+zetastar
    zetagasmc=dblarr(nphysmc)+zetagas
    rpeakstarmc=dblarr(nphysmc)+rpeakstar
    rpeakgasmc=dblarr(nphysmc)+rpeakgas
    esfmc=dblarr(nphysmc)+esf
    mdotsfmc=dblarr(nphysmc)+mdotsf
    mdotfbmc=dblarr(nphysmc)+mdotfb
    vfbmc=dblarr(nphysmc)+vfb
    etainstmc=dblarr(nphysmc)+etainst
    etaavgmc=dblarr(nphysmc)+etaavg
    chiemc=dblarr(nphysmc)+chie
    chipmc=dblarr(nphysmc)+chip
    
    nrnd=10*nphysmc
    randomarr=randomu(systime(1),nrnd)
    gaussarr=randomn(systime(1),nrnd)
    irnd=long(0)
    igauss=long(0)
    
    for i=0L,nphysmc-1 do begin ;do Monte Carlo error propagation
        tstarisomc(i)=tstariso+gaussarr(igauss)*((gaussarr(igauss) lt 0.)*tstariso_rerrmin+(gaussarr(igauss) gt 0.)*tstariso_rerrmax)*tstariso
        xi=tstarisomc(i)/tstariso
        igauss+=1
        tgasmc(i)=10.^interpol(alog10(tgasarr3d),probtgascum3d,randomarr(irnd))*xi
        tovermc(i)=10.^interpol(alog10(toverarr3d),probtovercum3d,randomarr(irnd))*xi
        lambdamc(i)=10.^interpol(alog10(lambdaarr3d),problambdacum3d,randomarr(irnd))
        irnd+=1
        surfsfrmc(i)=surfsfr+gaussarr(igauss)*surfsfr_err
        igauss+=1
        surfgasmc(i)=surfgas+gaussarr(igauss)*surfgas_err
        igauss+=1
        
        tdeplmc(i)=surfgasmc(i)/surfsfrmc(i)
        if tstar_incl eq 0 then tstarmc(i)=tstariso*xi+tovermc(i) else tstarmc(i)=tstariso*xi
        ttotalmc(i)=tgasmc(i)+tstarmc(i)-tovermc(i)
        zetastarmc(i)=f_zeta(surfcontrasts,ttotalmc(i),tstarmc(i),peak_prof)
        zetagasmc(i)=f_zeta(surfcontrastg,ttotalmc(i),tgasmc(i),peak_prof)
        rpeakstarmc(i)=f_rpeak(zetastarmc(i),lambdamc(i))
        rpeakgasmc(i)=f_rpeak(zetagasmc(i),lambdamc(i))
        esfmc(i)=f_esf(tgasmc(i),tdeplmc(i),fcl,fgmc)
        mdotsfmc(i)=f_mdotsf(tgasmc(i),lambdamc(i),surfgasmc(i),fgmc,esfmc(i))
        mdotfbmc(i)=f_mdotfb(tovermc(i),lambdamc(i),surfgasmc(i),fgmc,esfmc(i))
        vfbmc(i)=f_vfb(tovermc(i),lambdamc(i))
        etainstmc(i)=f_etainst(tgasmc(i),tovermc(i),esfmc(i))
        etaavgmc(i)=f_etaavg(esfmc(i))
        chiemc(i)=f_chie(tovermc(i),esfmc(i),vfbmc(i),psi)
        chipmc(i)=f_chip(tovermc(i),surfgasmc(i),esfmc(i),vfbmc(i),psi,phitrap,kappa0)
        
        progress,'     ==> derived quantity error propagation progress',i,nphysmc-1
    endfor
    
    ;sort all arrays
    tgasmc=tgasmc(sort(tgasmc))
    tovermc=tovermc(sort(tovermc))
    lambdamc=lambdamc(sort(lambdamc))
    surfsfrmc=surfsfrmc(sort(surfsfrmc))
    surfgasmc=surfgasmc(sort(surfgasmc))
    tstarisomc=tstarisomc(sort(tstarisomc))
    
    tdeplmc=tdeplmc(sort(tdeplmc))
    dtdeplmc=f_getdvar(tdeplmc)
    tstarmc=tstarmc(sort(tstarmc))
    dtstarmc=f_getdvar(tstarmc)
    ttotalmc=ttotalmc(sort(ttotalmc))
    dttotalmc=f_getdvar(ttotalmc)
    zetastarmc=zetastarmc(sort(zetastarmc))
    dzetastarmc=f_getdvar(zetastarmc)
    zetagasmc=zetagasmc(sort(zetagasmc))
    dzetagasmc=f_getdvar(zetagasmc)
    rpeakstarmc=rpeakstarmc(sort(rpeakstarmc))
    drpeakstarmc=f_getdvar(rpeakstarmc)
    rpeakgasmc=rpeakgasmc(sort(rpeakgasmc))
    drpeakgasmc=f_getdvar(rpeakgasmc)
    esfmc=esfmc(sort(esfmc))
    desfmc=f_getdvar(esfmc)
    mdotsfmc=mdotsfmc(sort(mdotsfmc))
    dmdotsfmc=f_getdvar(mdotsfmc)
    mdotfbmc=mdotfbmc(sort(mdotfbmc))
    dmdotfbmc=f_getdvar(mdotfbmc)
    vfbmc=vfbmc(sort(vfbmc))
    dvfbmc=f_getdvar(vfbmc)
    etainstmc=etainstmc(sort(etainstmc))
    detainstmc=f_getdvar(etainstmc)
    etaavgmc=etaavgmc(sort(etaavgmc))
    detaavgmc=f_getdvar(etaavgmc)
    chiemc=chiemc(sort(chiemc))
    dchiemc=f_getdvar(chiemc)
    chipmc=chipmc(sort(chipmc))
    dchipmc=f_getdvar(chipmc)
    
    ;get error bars on best-fitting values
    cumdistr=findgen(nphysmc)/(nphysmc-1.)
    dummy=f_pdftovalues(tdeplmc,dtdeplmc,cumdistr)
    tdepl_errmin=dummy(1)
    tdepl_errmax=dummy(2)
    dummy=f_pdftovalues(tstarmc,dtstarmc,cumdistr)
    tstar_errmin=dummy(1)
    tstar_errmax=dummy(2)
    dummy=f_pdftovalues(ttotalmc,dttotalmc,cumdistr)
    ttotal_errmin=dummy(1)
    ttotal_errmax=dummy(2)
    dummy=f_pdftovalues(zetastarmc,dzetastarmc,cumdistr)
    zetastar_errmin=dummy(1)
    zetastar_errmax=dummy(2)
    dummy=f_pdftovalues(zetagasmc,dzetagasmc,cumdistr)
    zetagas_errmin=dummy(1)
    zetagas_errmax=dummy(2)
    dummy=f_pdftovalues(rpeakstarmc,drpeakstarmc,cumdistr)
    rpeakstar_errmin=dummy(1)
    rpeakstar_errmax=dummy(2)
    dummy=f_pdftovalues(rpeakgasmc,drpeakgasmc,cumdistr)
    rpeakgas_errmin=dummy(1)
    rpeakgas_errmax=dummy(2)
    dummy=f_pdftovalues(esfmc,desfmc,cumdistr)
    esf_errmin=dummy(1)
    esf_errmax=dummy(2)
    dummy=f_pdftovalues(mdotsfmc,dmdotsfmc,cumdistr)
    mdotsf_errmin=dummy(1)
    mdotsf_errmax=dummy(2)
    dummy=f_pdftovalues(mdotfbmc,dmdotfbmc,cumdistr)
    mdotfb_errmin=dummy(1)
    mdotfb_errmax=dummy(2)
    dummy=f_pdftovalues(vfbmc,dvfbmc,cumdistr)
    vfb_errmin=dummy(1)
    vfb_errmax=dummy(2)
    dummy=f_pdftovalues(etainstmc,detainstmc,cumdistr)
    etainst_errmin=dummy(1)
    etainst_errmax=dummy(2)
    dummy=f_pdftovalues(etaavgmc,detaavgmc,cumdistr)
    etaavg_errmin=dummy(1)
    etaavg_errmax=dummy(2)
    dummy=f_pdftovalues(chiemc,dchiemc,cumdistr)
    chie_errmin=dummy(1)
    chie_errmax=dummy(2)
    dummy=f_pdftovalues(chipmc,dchipmc,cumdistr)
    chip_errmin=dummy(1)
    chip_errmax=dummy(2)
            
    ;plot PDFs and write tables
    print,'     ==> plotting probability distribution functions and writing them to output directory'
    dummy=f_createpdf(tdeplmc,dtdeplmc,cumdistr,ntry)
    tdeplarr=dummy(*,0)
    dtdepl=dummy(*,1)
    probtdepl=dummy(*,2)
    report=f_plotdistr(tdeplarr,dtdepl,probtdepl,tdepl,galaxy,figdir,'tdepl','!8t!6!Ddepl!N','!6yr')
    report=f_writepdf(alog10(tdeplarr),alog10(dtdepl),alog10(probtdepl),galaxy,outputdir,'tdepl','# log10(tdepl[yr]), log10(dtdepl[yr]), log10(PDF[yr^-1])')
    dummy=f_createpdf(tstarmc,dtstarmc,cumdistr,ntry)
    tstararr=dummy(*,0)
    dtstar=dummy(*,1)
    probtstar=dummy(*,2)
    report=f_plotdistr(tstararr,dtstar,probtstar,tstar,galaxy,figdir,'tstar','!8t!6!Dstar!N','!6Myr')
    report=f_writepdf(alog10(tstararr),alog10(dtstar),alog10(probtstar),galaxy,outputdir,'tstar','# log10(tstar[Myr]), log10(dtstar[Myr]), log10(PDF[Myr^-1])')
    dummy=f_createpdf(ttotalmc,dttotalmc,cumdistr,ntry)
    ttotalarr=dummy(*,0)
    dttotal=dummy(*,1)
    probttotal=dummy(*,2)
    report=f_plotdistr(ttotalarr,dttotal,probttotal,ttotal,galaxy,figdir,'ttotal','!7s!6','!6Myr')
    report=f_writepdf(alog10(ttotalarr),alog10(dttotal),alog10(probttotal),galaxy,outputdir,'ttotal','# log10(ttotal[Myr]), log10(dttotal[Myr]), log10(PDF[Myr^-1])')
    dummy=f_createpdf(zetastarmc,dzetastarmc,cumdistr,ntry)
    zetastararr=dummy(*,0)
    dzetastar=dummy(*,1)
    probzetastar=dummy(*,2)
    report=f_plotdistr(zetastararr,dzetastar,probzetastar,zetastar,galaxy,figdir,'zetastar','!7f!6!Dstar!N','')
    report=f_writepdf(alog10(zetastararr),alog10(dzetastar),alog10(probzetastar),galaxy,outputdir,'zetastar','# log10(zetastar), log10(dzetastar), log10(PDF)')
    dummy=f_createpdf(zetagasmc,dzetagasmc,cumdistr,ntry)
    zetagasarr=dummy(*,0)
    dzetagas=dummy(*,1)
    probzetagas=dummy(*,2)
    report=f_plotdistr(zetagasarr,dzetagas,probzetagas,zetagas,galaxy,figdir,'zetagas','!7f!6!Dgas!N','')
    report=f_writepdf(alog10(zetagasarr),alog10(dzetagas),alog10(probzetagas),galaxy,outputdir,'zetagas','# log10(zetagas), log10(dzetagas), log10(PDF)')
    dummy=f_createpdf(rpeakstarmc,drpeakstarmc,cumdistr,ntry)
    rpeakstararr=dummy(*,0)
    drpeakstar=dummy(*,1)
    probrpeakstar=dummy(*,2)
    report=f_plotdistr(rpeakstararr,drpeakstar,probrpeakstar,rpeakstar,galaxy,figdir,'rpeakstar','!8r!6!Dstar!N','!6pc')
    report=f_writepdf(alog10(rpeakstararr),alog10(drpeakstar),alog10(probrpeakstar),galaxy,outputdir,'rpeakstar','# log10(rpeakstar[pc]), log10(drpeakstar[pc]), log10(PDF[pc^-1])')
    dummy=f_createpdf(rpeakgasmc,drpeakgasmc,cumdistr,ntry)
    rpeakgasarr=dummy(*,0)
    drpeakgas=dummy(*,1)
    probrpeakgas=dummy(*,2)
    report=f_plotdistr(rpeakgasarr,drpeakgas,probrpeakgas,rpeakgas,galaxy,figdir,'rpeakgas','!8r!6!Dgas!N','!6pc')
    report=f_writepdf(alog10(rpeakgasarr),alog10(drpeakgas),alog10(probrpeakgas),galaxy,outputdir,'rpeakgas','# log10(rpeakgas[pc]), log10(drpeakgas[pc]), log10(PDF[pc^-1])')
    dummy=f_createpdf(esfmc,desfmc,cumdistr,ntry)
    esfarr=dummy(*,0)
    desf=dummy(*,1)
    probesf=dummy(*,2)
    report=f_plotdistr(esfarr,desf,probesf,esf,galaxy,figdir,'esf','!7e!6!Dsf!N','')
    report=f_writepdf(alog10(esfarr),alog10(desf),alog10(probesf),galaxy,outputdir,'esf','# log10(esf), log10(desf), log10(PDF)')
    dummy=f_createpdf(mdotsfmc,dmdotsfmc,cumdistr,ntry)
    mdotsfarr=dummy(*,0)
    dmdotsf=dummy(*,1)
    probmdotsf=dummy(*,2)
    report=f_plotdistr(mdotsfarr,dmdotsf,probmdotsf,mdotsf,galaxy,figdir,'mdotsf','!6(d!8M!6/d!8t!6)!Dsf!N','!6M!D!9n!6!N yr!U-1!N')
    report=f_writepdf(alog10(mdotsfarr),alog10(dmdotsf),alog10(probmdotsf),galaxy,outputdir,'mdotsf','# log10(mdotsf[Msun yr^-1]), log10(dmdotsf[Msun yr^-1]), log10(PDF[Msun^-1 yr])')
    dummy=f_createpdf(mdotfbmc,dmdotfbmc,cumdistr,ntry)
    mdotfbarr=dummy(*,0)
    dmdotfb=dummy(*,1)
    probmdotfb=dummy(*,2)
    report=f_plotdistr(mdotfbarr,dmdotfb,probmdotfb,mdotfb,galaxy,figdir,'mdotfb','!6(d!8M!6/d!8t!6)!Dfb!N','!6M!D!9n!6!N yr!U-1!N')
    report=f_writepdf(alog10(mdotfbarr),alog10(dmdotfb),alog10(probmdotfb),galaxy,outputdir,'mdotfb','# log10(mdotfb[Msun yr^-1]), log10(dmdotfb[Msun yr^-1]), log10(PDF[Msun^-1 yr])')
    dummy=f_createpdf(vfbmc,dvfbmc,cumdistr,ntry)
    vfbarr=dummy(*,0)
    dvfb=dummy(*,1)
    probvfb=dummy(*,2)
    report=f_plotdistr(vfbarr,dvfb,probvfb,vfb,galaxy,figdir,'vfb','!8v!6!Dfb!N','!6km s!U-1!N')
    report=f_writepdf(alog10(vfbarr),alog10(dvfb),alog10(probvfb),galaxy,outputdir,'vfb','# log10(vfb[km s^-1]), log10(dvfb[km s^-1]), log10(PDF[km^-1 s])')
    dummy=f_createpdf(etainstmc,detainstmc,cumdistr,ntry)
    etainstarr=dummy(*,0)
    detainst=dummy(*,1)
    probetainst=dummy(*,2)
    report=f_plotdistr(etainstarr,detainst,probetainst,etainst,galaxy,figdir,'etainst','!7g!6!Dfb!N','')
    report=f_writepdf(alog10(etainstarr),alog10(detainst),alog10(probetainst),galaxy,outputdir,'etainst','# log10(etainst), log10(detainst), log10(PDF)')
    dummy=f_createpdf(etaavgmc,detaavgmc,cumdistr,ntry)
    etaavgarr=dummy(*,0)
    detaavg=dummy(*,1)
    probetaavg=dummy(*,2)
    report=f_plotdistr(etaavgarr,detaavg,probetaavg,etaavg,galaxy,figdir,'etaavg','!6<!7g!6!Dfb!N>','')
    report=f_writepdf(alog10(etaavgarr),alog10(detaavg),alog10(probetaavg),galaxy,outputdir,'etaavg','# log10(etaavg), log10(detaavg), log10(PDF)')
    dummy=f_createpdf(chiemc,dchiemc,cumdistr,ntry)
    chiearr=dummy(*,0)
    dchie=dummy(*,1)
    probchie=dummy(*,2)
    report=f_plotdistr(chiearr,dchie,probchie,chie,galaxy,figdir,'chie','!7v!6!Dfb,!8E!6!N','')
    report=f_writepdf(alog10(chiearr),alog10(dchie),alog10(probchie),galaxy,outputdir,'chie','# log10(chie), log10(dchie), log10(PDF)')
    dummy=f_createpdf(chipmc,dchipmc,cumdistr,ntry)
    chiparr=dummy(*,0)
    dchip=dummy(*,1)
    probchip=dummy(*,2)
    report=f_plotdistr(chiparr,dchip,probchip,chip,galaxy,figdir,'chip','!7v!6!Dfb,!8p!6!N','')
    report=f_writepdf(alog10(chiparr),alog10(dchip),alog10(probchip),galaxy,outputdir,'chip','# log10(chip), log10(dchip), log10(PDF)')        
    
    return,[tdepl,tdepl_errmin,tdepl_errmax, $
            tstar,tstar_errmin,tstar_errmax, $
            ttotal,ttotal_errmin,ttotal_errmax, $
            zetastar,zetastar_errmin,zetastar_errmax, $
            zetagas,zetagas_errmin,zetagas_errmax, $
            rpeakstar,rpeakstar_errmin,rpeakstar_errmax, $
            rpeakgas,rpeakgas_errmin,rpeakgas_errmax, $
            esf,esf_errmin,esf_errmax, $
            mdotsf,mdotsf_errmin,mdotsf_errmax, $
            mdotfb,mdotfb_errmin,mdotfb_errmax, $
            vfb,vfb_errmin,vfb_errmax, $
            etainst,etainst_errmin,etainst_errmax, $
            etaavg,etaavg_errmin,etaavg_errmax, $
            chie,chie_errmin,chie_errmax, $
            chip,chip_errmin,chip_errmax]
    
end







