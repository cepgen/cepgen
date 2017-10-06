        subroutine F2_fit_luxlike(xbj,q2,F2,FL)
c       input: x,q2
c       output: F2,F1
        implicit real*8 (a-h,o-z)
        common/luxlike_params/amp,am_pi,alpha_em,
     &     q2_cut,w2_lo,w2_hi,
     &     ires_model,icont_model
c       -----------------------------
        w_thr = amp+am_pi
        w2 = amp**2 + q2*(1.d0-xbj)/xbj
        w = dsqrt(w2)
c       -----------------------------
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
           if(ires_model.eq.1) then
              call CepGen_F2_ChristyBosted(xbj,Q2,F2r,FLr)
           elseif(ires_model.eq.2) then
              call CepGen_F2_FioreBrasse(xbj,Q2,F2r,FLr)
           endif
           F2 = F2r
           FL = FLr
        elseif(q2.lt.q2_cut.and.w2.gt.w2_lo.and.w2.lt.w2_hi) then
           if(ires_model.eq.1) then
              call CepGen_F2_ChristyBosted(xbj,Q2,F2r,FLr)
           elseif(ires_model.eq.2) then
              call CepGen_F2_FioreBrasse(xbj,Q2,F2r,FLr)
           endif
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
        common/luxlike_params/amp,am_pi,alpha_em,
     &     q2_cut,w2_lo,w2_hi,
     &     ires_model,icont_model
        if(icont_model.eq.1) then
           call CepGen_F2_GD11P(xbj,Q2,F2)
        elseif(icont_model.eq.2) then
           call CepGen_F2_ALLM91(xbj,Q2,F2)
        elseif(icont_model.eq.3) then
           call CepGen_F2_ALLM97(xbj,Q2,F2)
        endif
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
