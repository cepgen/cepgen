        subroutine F2_fit_luxlike(xbj,q2,F2,FL)
c       input: x,q2
c       output: F2,F1
        implicit real*8 (a-h,o-z)
c       -----------------------------        
        amp = 0.9382727
        am_pi = 0.135
        alpha_em = 1.d0/137.d0
c       -----------------------------
        w_thr = amp+am_pi
        q2_cut = 9.d0
        W2_hi = 4.d0
        W2_lo = 3.d0

        w2 = amp**2 + q2*(1.d0-xbj)/xbj
        w = dsqrt(w2)

        omega = (w2-w2_lo)/(w2_hi-w2_lo)
        rho = 2.d0*omega**2 - omega**4 


        if(q2.gt.q2_cut.and.w2.gt.w2_hi) then
                        call F2_perturbative(xbj,q2,F2p,FLp)
                        F2 = F2p
                        FL = FLp        
        elseif(q2.gt.q2_cut.and.w2.lt.w2_hi) then
                        call F2_cont(xbj,q2,F2c,FLc)
                        F2 = F2c
                        FL = FLc
        elseif(q2.lt.q2_cut.and.w2.lt.w2_lo) then
                        call CepGen_F2_ChristyBosted(xbj,Q2,F2r,FLr)
                        F2 = F2r
                        FL = FLr
        elseif(q2.lt.q2_cut.and.w2.gt.w2_lo.and.w2.lt.w2_hi) then
                        call CepGen_F2_ChristyBosted(xbj,Q2,F2r,FLr)
                        call F2_cont(xbj,q2,F2c,FLc)
                        F2 = (1.d0-rho)*F2r + rho*F2c
                        FL=  (1.d0-rho)*FLr + rho*FLc        
        elseif(q2.lt.q2_cut.and.w2.gt.w2_hi) then
                        call F2_cont(xbj,q2,F2c,FLc)
                        F2 = F2c
                        FL = FLc
        endif
               
        return
        end
c       ---------------------------------------------------------
        subroutine F2_CONT(xbj,q2,F2c,FLc)
        implicit real*8 (a-h,o-z)
        amp = 0.938d0

        call CepGen_F2_GD11P(xbj,Q2,F2)
        F2c = F2
        R = R_1998(xbj,Q2)
        tau = 4.d0*xbj**2*amp**2/q2
        FLc = F2c*(1+tau)*R/(1.d0+R)

        return
        end        
c       ---------------------------------------------------------       
        subroutine F2_perturbative(xbj,q2,F2pert,FLpert)
      IMPLICIT NONE
      INTEGER iord,ipn,ieigen,neigen,ix,iq
      PARAMETER(neigen=20)
      COMMON/iordCommon/iord
      DOUBLE PRECISION ALPHAS,x,q,GetOnePDF,xg,
     &     f2,f2c,f2b,fl,flc,flb,
     &     f2_p,f2c_p,f2b_p,fl_p,flc_p,flb_p,
     &     f2_m,f2c_m,f2b_m,fl_m,flc_m,flb_m,
     &     f2_max,f2c_max,f2b_max,fl_max,flc_max,flb_max,
     &     f2_min,f2c_min,f2b_min,fl_min,flc_min,flb_min
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION factor_HT,xbj,q2,F2pert,FLpert
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      INTEGER ISET
      CHARACTER prefix*50,cl*4,pdfset*50
      COMMON/mstwfiles/iset,prefix,cl

      CALL WATE96               ! needed for Gaussian integration

      IORD = 2                  ! 0 = LO, 1 = NLO, 2 = NNLO
      IPN = 1                   ! 1 = proton, 2 = neutron

      IF (IORD.EQ.0) THEN
         pdfset = "mstw2008lo"
      ELSE IF (IORD.EQ.1) THEN
         pdfset = "mstw2008nlo"
      ELSE IF (IORD.EQ.2) THEN
         pdfset = "mstw2008nnlo"
      END IF
      prefix = "Grids/"//pdfset
      cl = '68cl'
C      cl = '90cl'
      iset = 0                  ! central fit
C--   Heavy quark masses and alphaS values
C--   are stored in a COMMON block.
      xg = GetOnePDF(prefix,iset,0.1d0,10d0,0) ! dummy call
c      WRITE(6,*) "mCharm = ",mCharm,", mBottom = ",mBottom
c      WRITE(6,*) "alphaS(Q0) = ",alphaSQ0,", alphaS(MZ) = ",
c     &     alphaSMZ,", alphaSorder = ",alphaSorder,
c     &     ", alphaSnfmax = ",alphaSnfmax

C--   Call the initialisation routine with alpha_S(Q_0).
      CALL INITALPHAS(alphaSorder,1.D0,1.D0,alphaSQ0,
     &     mCharm,mBottom,1.D10)
C--   Check calculated value of alpha_S(M_Z) matches stored value.
cc      WRITE(6,'(" alphaS(MZ) = ",F7.5," = ",F7.5)')
cc     &     ALPHAS(91.1876D0),alphaSMZ


c       in the Lux-paper a "higher-twist" correction is applied to F2

        factor_HT = 1.d0 + 5.5d0/q2
        q = dsqrt(q2)
        x = xbj
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        F2pert = F2*factor_HT
        FLpert = FL
        return
        end
c       --------------------------------------------------------                
        subroutine F2_ALLM(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        real*8 q2
c       From  Abramowicz, Levy et al.
c       --------------------------
        data cp1,cp2,cp3/0.28076d0,0.22291d0,2.1979d0/
        data ap1,ap2,ap3/-0.0808d0,-0.44812d0,1.1709d0/
        data bp1,bp2,bp3/0.36292d0,1.8917d0,1.8439d0/
c       --------------------------
        data cr1,cr2,cr3/0.80107d0,0.97307d0,3.4924d0/
        data ar1,ar2,ar3/0.58400d0,0.37888d0,2.6063d0/
        data br1,br2,br3/0.01147d0,3.7582d0,0.49338d0/
c       --------------------------
        data am02,amp2,amr2/0.31985d0,49.457d0,0.15052d0/
        data q02,alam2/0.52544d0,0.06526d0/
c       --------------------------
        factor = Q2/(Q2+am02)
        W2_eff = Q2*(1.d0-x)/x
        xp = (Q2+amp2)/(Q2+W2_eff+amp2)
        xr = (Q2+amr2)/(Q2+W2_eff+amr2)
c       ---------------------------
        xlog1 = dlog((Q2+q02/alam2))
        xlog2 = dlog(q02/alam2)
        t = dlog(xlog1/xlog2)
c       ----------------------------
        cpom = cp1 + (cp1-cp2)*(1.d0/(1.d0 + t**cp3) - 1.d0)
        apom = ap1 + (ap1-ap2)*(1.d0/(1.d0 + t**ap3) - 1.d0)
        creg = cr1 + cr2*t**cr3        
        areg = ar1 + ar2*t**ar3        
c        bpom = bp1**2 + bp2**2 * t**bp3
c        breg = br1**2 + br2**2 * t**br3
        bpom = bp1 + bp2 * t**bp3
        breg = br1 + br2 * t**br3

c       -----------------------------
c        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**(1.+bpom)
c        F2_Reg = factor*creg*xr**areg *(1.d0-x)**(1.+breg)

        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**bpom
        F2_Reg = factor*creg*xr**areg *(1.d0-x)**breg

cc        F2 = F2_Pom
        F2 = F2_Pom + F2_Reg

        return
        end

c       -------------------------------------------------

        subroutine F2_GD11P(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        real*8 q2
c       Refit of ALLM by the HERMES collaboration
c       --------------------------
        data cp1,cp2,cp3/0.3638d0,0.1211d0,1.166d0/
        data ap1,ap2,ap3/-0.11895,-0.4783d0,1.353d0/
        data bp1,bp2,bp3/1.0833d0,2.656d0,1.771d0/
c       --------------------------
        data cr1,cr2,cr3/1.3633d0,2.256d0,2.209d0/
        data ar1,ar2,ar3/0.3425d0,1.0603d0,0.5164d0/
        data br1,br2,br3/-10.408d0,14.857d0,0.07739d0/
c       --------------------------
        data am02,amp2,amr2/0.5063d0,34.75d0,0.03190/
        data q02,alam2/1.374d0,0.06527d0/
c       --------------------------
        factor = Q2/(Q2+am02)
        W2_eff = Q2*(1.d0-x)/x
        xp = (Q2+amp2)/(Q2+W2_eff+amp2)
        xr = (Q2+amr2)/(Q2+W2_eff+amr2)
c       ---------------------------
        xlog1 = dlog((Q2+q02)/alam2)
        xlog2 = dlog(q02/alam2)
        t = dlog(xlog1/xlog2)
c       ----------------------------
        cpom = cp1 + (cp1-cp2)*(1.d0/(1.d0 + t**cp3) - 1.d0)
        apom = ap1 + (ap1-ap2)*(1.d0/(1.d0 + t**ap3) - 1.d0)
        creg = cr1 + cr2*t**cr3        
        areg = ar1 + ar2*t**ar3        
c        bpom = bp1**2 + bp2**2 * t**bp3
c        breg = br1**2 + br2**2 * t**br3
        bpom = bp1 + bp2 * t**bp3
        breg = br1 + br2 * t**br3

c       -----------------------------
c        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**(1.+bpom)
c        F2_Reg = factor*creg*xr**areg *(1.d0-x)**(1.+breg)

        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**bpom
        F2_Reg = factor*creg*xr**areg *(1.d0-x)**breg

cc        F2 = F2_Pom
        F2 = F2_Pom + F2_Reg

        return
        end
c       -------------------------------------------------
        double precision function R_1998(x,q2)
        implicit real*8 (a-h,o-z)

c       from hep-ex/9808028v1, K.Abe et al.(SLAC)

        data a1,a2,a3/0.0485d0,0.5470d0,2.0621d0/
        data a4,a5,a6/-0.3804d0,0.5090d0,-0.0285d0/
        data b1,b2,b3/0.0481d0,0.6114d0,-0.3509d0/
        data b4,b5,b6/-0.4611d0,0.7172d0,-0.0317d0/ 
        data c1,c2,c3/0.0577d0,0.4644d0,1.8288d0/
        data c4,c5,c6/12.3708d0,-43.1043d0,41.7415d0/

        q2_b = 0.34
        u = q2/q2_b
        xl = dlog(q2/0.04d0)

        pa = (1.d0 + a4*x +a5*x**2)*x**a6
        pb = (1.d0 + b4*x +b5*x**2)*x**b6
        tt = theta(x,q2)
        q2_thr = c4*x +c5*x**2+c6*x**3

        ra = a1/xl*tt+a2/(q2**4.d0+a3**4.d0)**.25d0*pa
        rb = b1/xl*tt+(b2/q2+b3/(q2**2+.3d0**2))*pb
        rc = c1/xl*tt+c2*((q2-q2_thr)**2 +c3**2)**(-.5d0)

       
        R =(ra+rb+rc)/3.d0

        if(q2.gt.q2_b) then
                R_1998 = R
        else
                R_1998 = R*(3.d0*u - u**3)/2.d0
        endif

        return
        end

c       -------------------------------------------------
        double precision function theta(x,q2)
c       
c       function needed for the 1998 parametrisation of R=sigma_L/sigma_T
c
        implicit real*8(a-h,o-z)
        theta = 1.d0 + 12.d0*q2/(q2+1.d0)
     >                  *(0.125d0**2/(.125d0**2 +x**2))
        return
        end
c       --------------------------------------------------
