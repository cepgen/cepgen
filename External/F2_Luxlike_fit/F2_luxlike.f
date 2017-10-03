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
                        call F2_res(xbj,q2,F2r,FLr)
                        F2 = F2r
                        FL = FLr
           elseif(q2.lt.q2_cut.and.w2.gt.w2_lo.and.w2.lt.w2_hi) then
                        call F2_res(xbj,q2,F2r,FLr)
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
c       -----------------------------------------------------
        subroutine F2_res(xbj,q2,F2,FL)
        implicit real*8 (a-h,o-z)
c       -----------------------------
c       Resonance contribution at low W2.
c       This one uses the Christy-Bosted fit, one could also use Fiore et al
c       here. Maybe add a choice later??
c       -----------------------------        
        amp = 0.9382727
        am_pi = 0.135
        w_thr = amp+am_pi
        alpha_em = 1.d0/137.d0
        pi = dacos(-1.d0)
c       -----------------------------
        w2 = amp**2 + q2*(1.d0-xbj)/xbj
        w = dsqrt(w2)
        if(w.lt.w_thr) then
                F2= 0.d0
                FL=0.d0
                return
        endif
c       -----------------------------
c       modification of Christy-Bosted at large q2 as described in
c       the LUXqed paper
c       -----------------------------
        q2_1 = 30.d0
        q2_0 = 8.d0
        delq2 = q2 - q2_0
        qq = q2_1 - q2_0
        factor_mod = q2_1/(q2_1 + delq2)
        q2_mod = q2_0 + delq2/(1.d0+delq2/qq)
c       ------------------------------

        if(q2.lt.q2_0) then
                tau = 4.d0*xbj**2*amp**2/q2
                prefac = 1.d0/(4.d0*pi**2*alpha_em)
     >                          *q2*(1.d0-xbj)/(1+tau)
                call christy507(w2,q2,f1,R,sigt,sigl)
                F2 = prefac*(sigt+sigl)/0.3894e3
                FL = F2*(1+tau)*R/(1.d0+R)
        else
                w2 = amp**2 + q2_mod*(1.d0-xbj)/xbj
                tau = 4.d0*xbj**2*amp**2/q2_mod
                prefac = 1.d0/(4.d0*pi**2*alpha_em)
     >                          *q2_mod*(1.d0-xbj)/(1+tau)
                call christy507(w2,q2_mod,f1,R,sigt,sigl)
                F2 = prefac*(sigt+sigl)/0.3894e3*factor_mod
                FL = F2*(1+tau)*R/(1.d0+R)*factor_mod
        endif

        return
        end
c       ---------------------------------------------------------
        subroutine F2_CONT(xbj,q2,F2c,FLc)
        implicit real*8 (a-h,o-z)
        amp = 0.938d0

        call F2_GD11P(xbj,Q2,F2)
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
        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**(1.+bpom)
        F2_Reg = factor*creg*xr**areg *(1.d0-x)**(1.+breg)

c        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**bpom
c        F2_Reg = factor*creg*xr**areg *(1.d0-x)**breg

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
        data ap1,ap2,ap3/-0.11859,-0.4783d0,1.353d0/
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
        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**(1.+bpom)
        F2_Reg = factor*creg*xr**areg *(1.d0-x)**(1.+breg)

c        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**bpom
c        F2_Reg = factor*creg*xr**areg *(1.d0-x)**breg

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
c cannibalized from:
c Fit to proton F1, R, sigma_T, and Sigma_L from

      SUBROUTINE christy507(W2,Q2,F1,R,sigt,sigl)
c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision 
c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
c    (arXiv:0712.3731). To be submitted to Phys. Rev. C.

      IMPLICIT NONE

      double precision w2,q2,xval1(50),xvall(50),xval(100)
      double precision mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end

 
      SUBROUTINE RESMOD507(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      DOUBLE PRECISION W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig
      double precision xval(50),mass(7),width(7)
      DOUBLE PRECISION height(7),sig_del,sig_21,sig_22,sig_31,sig_32
      double precision rescoef(6,4)
      DOUBLE PRECISION nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      DOUBLE PRECISION mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7)
      double precision ppicm,ppi2cm
      DOUBLE PRECISION petacm,ppicmr(7),ppi2cmr(7),petacmr(7)
      double precision epicmr(7),epi2cmr(7)
      DOUBLE PRECISION eetacmr(7),epicm,epi2cm,eetacm,br(7,3)
      double precision spin(7),ang(7)
      DOUBLE PRECISION pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      DOUBLE PRECISION sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      double precision sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      double precision sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       


c 1000  format(8f12.5)

      RETURN 
      END 



