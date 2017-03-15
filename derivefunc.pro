;Functions for derived quantities

function f_zeta,surfcontrast,ttot,tref,peak_prof ;Tracer peak concentration parameter
    COMMON numbers
    if peak_prof eq 0 then zeta=tiny
    if peak_prof eq 1 then zeta=sqrt(ttot/(surfcontrast*tref))
    if peak_prof eq 2 then zeta=sqrt(ttot/(2.*surfcontrast*tref))
    return,zeta
end

function f_rpeak,zeta,lambda ;Tracer peak radius
    rpeak=0.5*zeta*lambda
    return,rpeak
end

function f_esf,tgas,tdepl,fcl,fgmc ;Star formation efficiency per star formation event
    esf=min([fcl*tgas/(fgmc*tdepl*1.d3),1.]) ;cannot exceed unity
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

function f_chie,tover,esf,vfb,psie ;Feedback energy efficiency
    COMMON numbers
    COMMON astrconst
    chie=max([(1.-esf)*(vfb*kms)^2./(2.*esf*(tover*myr)*psie),tiny])
    return,chie
end

function f_chip,tover,esf,vfb,psip ;Feedback momentum efficiency
    COMMON numbers
    COMMON physconst
    COMMON astrconst
    chip=max([(1.-esf)*(vfb*kms)/(esf*(tover*myr)*psip),tiny])
    return,chip
end