! FORTRAN MODULE FOR APPLYING THE UNCERTAINY PRINCIPLE OF 
! KRUIJSSEN & LONGMORE (2014), MNRAS ACCEPTED, ARXIV:1401.4459
! Copyright (C) 2013-2014 J. M. Diederik Kruijssen & Steven N. Longmore
!
! NAMES:
!       principle_mod (module)
!       F_DELTAX (function)
!       F_SCATTER (function)
!       F_RATIOBIAS_1 (function)
!       F_RATIOBIAS_2 (function)
!
! PURPOSE:
!       f_deltax - Predict spatial scale below which scaling relation fails 
!       f_scatter - Predict scatter on log(tracer ratio) (e.g. log(t_depl))
!       f_ratiobias_1 - Predict the bias of the tracer ratio (e.g. the gas
!                       depletion time-scale) when focussing on phase 1 (e.g.
!                       gas) tracer peaks
!       f_ratiobias_2 - Predict the bias of the tracer ratio (e.g. the gas
!                       depletion time-scale) when focussing on phase 2 (e.g.
!                       SF) tracer peaks
!
! TERMS OF USE:
!       If you use these routines while preparing a paper, please cite the paper
!       in which this model was presented:
!       Kruijssen J. M. D., Longmore S. N. (2014), MNRAS accepted, arXiV:1401.4459
!       (this file will be updated with the final MNRAS reference later)
!
!       If it is the first time you apply the uncertainty principle, the authors
!       would appreciate being involved in the project to provide assistance
!       where necessary. However, this is not required and you are free to use
!       these routines as you see fit.
!
! CALLING SEQUENCE:
!       f_deltax(t_1,t_2,t_over,lambda,sfr_min,surfsfr,veldis)
!       f_scatter(t_1,t_2,t_over,l_ap,lambda,beta_1,beta_2,sigma_1evo1,
!                 sigma_2evo1,sigma_mf1,sigma_obs)
!       f_ratiobias_1(t_1,t_2,t_over,l_ap,lambda,beta_2)
!       f_ratiobias_2(t_1,t_2,t_over,l_ap,lambda,beta_1)
!
! INPUT PARAMETERS USED:
!       NOTE: ALL INPUT SHOULD BE IN CONSISTENT UNITS!
!             E.G. IF PARSECS ARE USED TO EXPRESS SIZES, DO NOT EXPRESS THE SFR
!             DENSITY IN UNITS OF MSUN YR-1 KPC-2 BUT IN UNITS OF MSUN YR-1 PC-2
!       t_1 - duration of phase 1
!       t_2 - duration of phase 2
!       t_over - overlap between phases 1 and 2
!       lambda - mean separation between independent regions
!       sfr_min - minimum SFR for SF tracer to be well-sampled from IMF
!                 (specific for SF relations, set to 0 to omit IMF effect)
!       surfsfr - SFR surface density
!                 (specific for SF relations, do not set to 0)
!       veldis - velocity dispersion
!                 (set to 0 to omit drift effect)
!       l_ap - aperture size (diameter)
!       beta_1 - overlap-to-isolation flux ratio of phase 1 tracer
!       beta_2 - overlap-to-isolation flux ratio of phase 2 tracer
!       sigma_1evo1 - scatter due to luminosity evolution of phase 1 tracer for
!                     a single region
!       sigma_2evo1 - scatter due to luminosity evolution of phase 2 tracer for
!                     a single region
!       sigma_mf1 - scatter due to mass spectrum for a single region
!       sigma_obs - scatter due to intrinsic observational errors
!
! OUTPUT:
!       f_deltax - array containing 4 elements:
!                  deltax_max (scale below which scaling relation fails)
!                  deltax_samp (failure scale due to Poisson statistics)
!                  deltax_imf (failure scale due to IMF sampling)
!                  deltax_drift (failure scale due to spatial drift)
!       f_scatter - array containing 5 elements:
!                   scatter_tot (total scatter of log(tracer ratio))
!                   scatter_samp (due to Poisson statistics)
!                   scatter_evo (due to luminosity evolution)
!                   scatter_mf (due to mass spectrum)
!                   scatter_obs (due to intrinsic observational errors)
!       f_ratiobias_1 - single number indicating the bias of the tracer ratio 
!                       (e.g. the gas depletion time-scale) when focussing on
!                       phase 1 tracer peaks
!       f_ratiobias_2 - single number indicating the bias of the tracer ratio 
!                       (e.g. the gas depletion time-scale) when focussing on
!                       phase 2 tracer peaks
!
! EXAMPLE FOR ALL PROGRAMS ARE GIVEN IN THE FILE test_principle.F90:
!       Calculate the bias of the tracer ratio (e.g. the gas depletion
!       time-scale) when focussing on phase 1 (e.g. gas) tracer peaks
!
!       PROGRAM test_principle
!           use principle_mod
!           REAL t_1,t_2,t_over,l_ap,lambda,beta_2,ratiobias_1
!           t_1=10.e6
!           t_2=5.e6
!           t_over=3.e6
!           l_ap=50.
!           lambda=130.
!           beta_2=0.5
!           ratiobias_1=f_ratiobias_1(t_1,t_2,t_over,l_ap,lambda,beta_2)
!           PRINT*,ratiobias_1
!       END PROGRAM
!
!       > g95 f_principle.f90 test_principle.f90 -o test_principle
!       > ./test_principle
!          2.03548646
!
! MODULE STRUCTURE: 
!       The functions are listed below in the order of their appearance in the
!       documentation.
!
! REVISION HISTORY:
!       Written by Diederik Kruijssen, November 2013
!

MODULE principle_mod
    CONTAINS

FUNCTION f_deltax(t_1,t_2,t_over,lambda,sfr_min,surfsfr,veldis) result(deltax) !size-scale below which galactic SF relations fail
    REAL, INTENT(in) :: t_1,t_2,t_over,lambda,sfr_min,surfsfr,veldis
    REAL :: pi,tau,tshort,deltax_samp,deltax_imf,deltax_drift,deltax_max
    REAL :: deltax(1:4)
    pi=3.14159265358979323846 !pi
    tau=t_1+t_2-t_over !total duration of process
    tshort=MIN(t_1,t_2) !duration of shortest phase
    deltax_samp=sqrt(tau/tshort)*lambda !size-scale below which scaling relations fail due to inadequate Poisson sampling of independent regions
    deltax_imf=sqrt(4./pi*sfr_min/surfsfr) !size-scale below which galactic SF relations fail due to inadequate sampling of massive stars from the IMF
    deltax_drift=0.5*veldis*tau !size-scale below which scaling relations fail due to phase 1 and 2 tracers undergoing spatial drift
    deltax_max=MAX(deltax_samp,deltax_imf,deltax_drift) !maximum of the previous three
    deltax=[deltax_max,deltax_samp,deltax_imf,deltax_drift] !array containing deltax_max, deltax_samp, deltax_imf, and deltax_drift
END FUNCTION

FUNCTION f_scatter(t_1,t_2,t_over,l_ap,lambda,beta_1,beta_2,sigma_1evo1,sigma_2evo1,sigma_mf1,sigma_obs) result(scatter) !scatter in log(tracer ratio)
    REAL, INTENT(in) :: t_1,t_2,t_over,l_ap,lambda,beta_1,beta_2,sigma_1evo1,sigma_2evo1,sigma_mf1,sigma_obs
    REAL :: tau,t_1iso,t_2iso,n_1iso,n_2iso,n_over
    REAL :: exp_1o,exp_2o,denominator_1,denominator_2
    REAL :: dlntdn1,dlntdn2,dlntdno
    REAL :: p_1,p_2,p_3,p_tot,n_1,n_2,n_tot
    REAL :: scatter_samp,scatter_evo,scatter_mf,scatter_obs,scatter_tot
    REAL :: scatter(1:5)
    tau=t_1+t_2-t_over !total duration of process
    t_1iso=t_1-t_over !duration of isolated phase 1
    t_2iso=t_2-t_over !duration of isolated phase 2
    n_1iso=t_1iso/tau*(l_ap/lambda)**2. !expected number of peaks in isolated phase 1
    n_2iso=t_2iso/tau*(l_ap/lambda)**2. !expected number of peaks in isolated phase 2
    n_over=t_over/tau*(l_ap/lambda)**2. !expected number of peaks in the overlap phase
    exp_1o=EXP(-n_1iso-n_over) !often-used exponential depending on n_1iso and n_over
    exp_2o=EXP(-n_2iso-n_over) !often-used exponential depending on n_2iso and n_over
    denominator_1=beta_1*n_over+n_1iso*(1.-exp_2o) !often-used term in the derivatives of ln(tracer ratio) mainly containing phase 1-related variables
    denominator_2=beta_2*n_over+n_2iso*(1.-exp_1o) !often-used term in the derivatives of ln(tracer ratio) mainly containing phase 2-related variables
    dlntdn1=(1.-exp_2o)/denominator_1-n_2iso*exp_1o/denominator_2 !partial derivative of ln(tracer ratio) with respect to n_1iso
    dlntdn2=n_1iso*exp_2o/denominator_1-(1.-exp_1o)/denominator_2 !partial derivative of ln(tracer ratio) with respect to n_2iso
    dlntdno=(beta_1+n_1iso*exp_2o)/denominator_1-(beta_2+n_2iso*exp_1o)/denominator_2 !partial derivative of ln(tracer ratio) with respect to n_over
    scatter_samp=1./ALOG(10.)*SQRT(dlntdn1**2.*n_1iso+dlntdn2**2.*n_2iso+dlntdno**2.*n_over) !scatter due to inadequate Poisson sampling of independent regions
    p_1=EXP(-n_1iso-n_2iso-n_over) !probability of disallowed combination #1 (no isolated phase 1, no isolated phase 2, no overlap)
    p_2=EXP(-n_2iso-n_over)*(1.-EXP(-n_1iso)) !probability of disallowed combination #2 (no isolated phase 2, no overlap)
    p_3=EXP(-n_1iso-n_over)*(1.-EXP(-n_2iso)) !probability of disallowed combination #3 (no isolated phase 1, no overlap)
    p_tot=1.-p_1-p_2-p_3 !total probability of allowed combinations of n_1iso, n_2iso and n_over
    n_1=(n_1iso*(1.-exp_2o)+n_over)/p_tot !number of phase 1 regions in the aperture (including overlap), accounting for the selection bias of requiring both tracers to be present
    n_2=(n_2iso*(1.-exp_1o)+n_over)/p_tot !number of phase 2 regions in the aperture (including overlap), accounting for the selection bias of requiring both tracers to be present
    n_tot=(n_1iso*(1.-exp_2o)+n_2iso*(1.-exp_1o)+n_over)/p_tot !total number of regions in the aperture, accounting for the selection bias of requiring both tracers to be present
    scatter_evo=SQRT(sigma_1evo1**2./n_1+sigma_2evo1**2./n_2) !scatter due to the luminosity evolution of independent regions
    scatter_mf=sigma_mf1/SQRT(n_tot) !scatter due to the mass spectrum of independent regions
    scatter_obs=sigma_obs !scatter due to intrinsic observational errors
    scatter_tot=SQRT(scatter_samp**2.+scatter_evo**2.+scatter_mf**2.+scatter_obs**2.) !total scatter
    scatter=[scatter_tot,scatter_samp,scatter_evo,scatter_mf,scatter_obs] !array containing total scatter and its four constituent scatter terms
END FUNCTION

FUNCTION f_ratiobias_1(t_1,t_2,t_over,l_ap,lambda,beta_2) result(ratiobias_1) !bias of the tracer ratio (e.g. t_depl) when focussing apertures on phase 1 tracer peaks
!NOTE: the asymptote for l_ap=>0 is (t_2/t_over-1.)/beta_2+1.
    REAL, INTENT(in) :: t_1,t_2,t_over,l_ap,lambda,beta_2
    REAL :: tau,focus1,focus2,random
    REAL :: ratiobias_1
    tau=t_1+t_2-t_over !total duration of process
    focus1=1. !contribution of central peak to phase 1 tracer flux
    focus2=beta_2*t_over/t_2/(1.+(beta_2-1.)*t_over/t_2) !contribution of central peak to phase 2 tracer flux
    random=(t_1/tau)*(l_ap/lambda)**2. !random regions covered by aperture
    ratiobias_1=(focus1+random)/(focus2+random) !bias of the tracer ratio (e.g. t_depl) when focussing apertures on phase 1 tracer peaks
END FUNCTION

FUNCTION f_ratiobias_2(t_1,t_2,t_over,l_ap,lambda,beta_1) result(ratiobias_2) !bias of the tracer ratio (e.g. t_depl) when focussing apertures on phase 2 tracer peaks
!NOTE: the asymptote for l_ap=>0 is [(t_1/t_over-1.)/beta_1+1.]**(-1.)
    REAL, INTENT(in) :: t_1,t_2,t_over,l_ap,lambda,beta_1
    REAL :: tau,focus1,focus2,random
    REAL :: ratiobias_2
    tau=t_1+t_2-t_over !total duration of process
    focus1=beta_1*t_over/t_1/(1.+(beta_1-1.)*t_over/t_1) !contribution of central peak to phase 1 tracer flux
    focus2=1. !contribution of central peak to phase 2 tracer flux
    random=(t_2/tau)*(l_ap/lambda)**2. !random regions covered by aperture
    ratiobias_2=(focus1+random)/(focus2+random) !bias of the tracer ratio (e.g. t_depl) when focussing apertures on phase 2 tracer peaks
END FUNCTION

END MODULE


