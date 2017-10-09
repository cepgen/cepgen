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

