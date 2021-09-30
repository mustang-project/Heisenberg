! TEST PROGRAM - READ THE HEADER OF F_PRINCIPLE.F90 FOR DETAILS

PROGRAM test_principle
    use principle_mod
    REAL t_1,t_2,t_over,lambda,sfr_min,surfsfr,veldis
    REAL l_ap,beta_1,beta_2,sigma_1evo1,sigma_2evo1,sigma_mf1,sigma_obs
    REAL deltax(1:4),scatter(1:5),ratiobias_1,ratiobias_2
    t_1=10.e6 !yr
    t_2=5.e6 !yr
    t_over=3.e6 !yr
    lambda=130. !pc
    sfr_min=1.e-3 !Msun/yr
    surfsfr=0.01e-6 !Msun/yr/pc**2
    veldis=1.e-5 !pc/yr
    l_ap=50. !pc
    beta_1=0.5
    beta_2=0.5
    sigma_1evo1=0.3
    sigma_2evo1=0.3
    sigma_mf1=0.8
    sigma_obs=0.15
    deltax=f_deltax(t_1,t_2,t_over,lambda,sfr_min,surfsfr,veldis)
    scatter=f_scatter(t_1,t_2,t_over,l_ap,lambda,beta_1,beta_2,sigma_1evo1,sigma_2evo1,sigma_mf1,sigma_obs)
    ratiobias_1=f_ratiobias_1(t_1,t_2,t_over,l_ap,lambda,beta_2)
    ratiobias_2=f_ratiobias_1(t_1,t_2,t_over,l_ap,lambda,beta_1)
    PRINT*,deltax(1:4)
    PRINT*,scatter(1:5)
    PRINT*,ratiobias_1
    PRINT*,ratiobias_2
END PROGRAM