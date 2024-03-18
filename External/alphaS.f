C----------------------------------------------------------------------
C--   Stand-alone code for alpha_s cannibalised (with permission)
C--   from Andreas Vogt QCD-PEGASUS package (hep-ph/0408244).
C--   The running coupling alpha_s is obtained at N^mLO (m = 0,1,2,3)
C--   by solving the renormalisation group equation in the MSbar scheme
C--   by a fourth-order Runge-Kutta integration.  Transitions from
C--   n_f to n_f+1 flavours are made when the factorisation scale
C--   mu_f equals the pole masses m_h (h = c,b,t).  At exactly
C--   the thresholds m_{c,b,t}, the number of flavours n_f = {3,4,5}.
C--   The top quark mass should be set to be very large to evolve with
C--   a maximum of five flavours.  The factorisation scale mu_f may be
C--   a constant multiple of the renormalisation scale mu_r.  The input
C--   factorisation scale mu_(f,0) should be less than or equal to
C--   the charm quark mass.  However, if it is greater than the
C--   charm quark mass, the value of alpha_s at mu_(f,0) = 1 GeV will
C--   be found using a root-finding algorithm.
C--
C--   Example of usage.
C--   First call the initialisation routine (only needed once):
C--
C--    IORD = 2                  ! perturbative order (N^mLO,m=0,1,2,3)
C--    FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2
C--    MUR = 1.D0                ! input mu_r in GeV
C--    ASMUR = 0.5D0             ! input value of alpha_s at mu_r
C--    MC = 1.4D0                ! charm quark mass
C--    MB = 4.75D0               ! bottom quark mass
C--    MT = 1.D10                ! top quark mass
C--    CALL INITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
C--
C--   Then get alpha_s at a renormalisation scale mu_r with:
C--
C--    MUR = 100.D0              ! renormalisation scale in GeV
C--    ALFAS = ALPHAS(MUR)
C--
C----------------------------------------------------------------------
C--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>
C----------------------------------------------------------------------

      SUBROUTINE INITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
C--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
C--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
C--   MUR = input renormalisation scale (in GeV) for alpha_s.
C--   ASMUR = input value of alpha_s at the renormalisation scale MUR.
C--   MC,MB,MT = heavy quark masses in GeV.
      IMPLICIT NONE
      INTEGER IORD,IORDc,MAXF,MODE
      DOUBLE PRECISION FR2,MUR,ASMUR,MC,MB,MT,EPS,A,B,DZEROX,
     &     R0c,FR2c,MURc,ASMURc,MCc,MBc,MTc,FINDALPHASR0,R0,ASI
      COMMON / DZEROXcommon / FR2c,MURc,ASMURc,MCc,MBc,MTc,R0c,IORDc
      PARAMETER(EPS=1.D-10,MAXF=10000,MODE=1)
      EXTERNAL FINDALPHASR0

      IF (MUR*sqrt(FR2).LE.MC) THEN ! Check that MUF <= MC.
         R0 = MUR
         ASI = ASMUR
      ELSE                      ! Solve for alpha_s at R0 = 1 GeV.
C--   Copy variables to common block.
         R0c = 1.D0/sqrt(FR2)
         IORDc = IORD
         FR2c = FR2
         MURc = MUR
         ASMURc = ASMUR
         MCc = MC
         MBc = MB
         MTc = MT
C--   Now get alpha_s(R0) corresponding to alpha_s(MUR).
         A = 0.02D0              ! lower bound for alpha_s(R0)
         B = 2.00D0              ! upper bound for alpha_s(R0)
         R0 = R0c
         ASI = DZEROX(A,B,EPS,MAXF,FINDALPHASR0,MODE)
      END IF

      CALL INITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)

      RETURN
      END

C----------------------------------------------------------------------

C--   Find the zero of this function using DZEROX.
      DOUBLE PRECISION FUNCTION FINDALPHASR0(ASI)
      IMPLICIT NONE
      INTEGER IORD
      DOUBLE PRECISION FR2, R0, ASI, MC, MB, MT, MUR, ASMUR, ALPHAS
      COMMON / DZEROXcommon / FR2, MUR, ASMUR, MC, MB, MT, R0, IORD

      CALL INITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)
      FINDALPHASR0 = ALPHAS(MUR) - ASMUR ! solve equal to zero

      RETURN
      END

C----------------------------------------------------------------------

      SUBROUTINE INITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)
C--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
C--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
C--   R0 = input renormalisation scale (in GeV) for alphas_s.
C--   ASI = input value of alpha_s at the renormalisation scale R0.
C--   MC,MB,MT = heavy quark masses in GeV.
C--   Must have R0*sqrt(FR2) <= MC to call this subroutine.
      IMPLICIT NONE
      INTEGER IORD,NAORD,NASTPS,IVFNS,NFF
      DOUBLE PRECISION FR2,R0,ASI,MC,MB,MT,LOGFR,R20,
     &     PI,ZETA,CF,CA,TR,AS0,M20,MC2,MB2,MT2
      PARAMETER(PI = 3.1415 92653 58979 D0)

      COMMON / RZETA  / ZETA(6)
      COMMON / COLOUR / CF, CA, TR
      COMMON / ASINP  / AS0, M20
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / VARFLV / IVFNS
      COMMON / NFFIX  / NFF
      COMMON / FRRAT  / LOGFR

c
c ..QCD colour factors
c
      CA = 3.D0
      CF = 4./3.D0
      TR = 0.5 D0
c
c ..The lowest integer values of the Zeta function
c
      ZETA(1) = 0.5772 15664 90153 D0
      ZETA(2) = 1.64493 40668 48226 D0
      ZETA(3) = 1.20205 69031 59594 D0
      ZETA(4) = 1.08232 32337 11138 D0
      ZETA(5) = 1.03692 77551 43370 D0
      ZETA(6) = 1.01734 30619 84449 D0

      IVFNS = 1                 ! variable flavour-number scheme (VFNS)
c      IVFNS = 0                 ! fixed flavour-number scheme (FFNS)
      NFF = 4                   ! number of flavours for FFNS
      NAORD = IORD              ! perturbative order of alpha_s
      NASTPS = 20               ! num. steps in Runge-Kutta integration
      R20 = R0**2               ! input renormalisation scale
      MC2 = MC**2               ! mu_f^2 for charm threshold
      MB2 = MB**2               ! mu_f^2 for bottom threshold
      MT2 = MT**2               ! mu_f^2 for top threshold
      LOGFR = LOG(FR2)          ! log of ratio of mu_f^2 to mu_r^2
      M20 = R20 * FR2           ! input factorisation scale

c
c ..Stop some nonsense
c
      IF ( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
c
      IF ( NAORD .GT. 3 ) THEN
         WRITE (6,*) 'Specified order in a_s too high. STOP'
         STOP
      END IF
c
      IF ( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001D0) ) THEN
         WRITE (6,*) 'Too low mu_r for VFNS evolution. STOP'
         STOP
      END IF
c
      IF ( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
         WRITE (6,*) 'Too high mu_0 for VFNS evolution. STOP'
         STOP
      END IF
c
      IF ( (ASI .GT. 2.D0) .OR. (ASI .LT. 2.D-2) ) THEN
         WRITE (6,*) 'alpha_s out of range. STOP'
         STOP
      END IF
c
      IF ( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
         WRITE (6,*) 'Wrong charm-bottom mass hierarchy. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
         WRITE (6,*) 'Wrong bottom-top mass hierarchy. STOP'
         STOP
      END IF
c

C--   Store the beta function coefficients in a COMMON block.
      CALL BETAFCT

C--   Store a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
      AS0 = ASI / (4.D0* PI)

C--   Store alpha_s at the heavy flavour thresholds in a COMMON block.
       IF (IVFNS .NE. 0) THEN
          CALL EVNFTHR (MC2, MB2, MT2)
       END IF

      RETURN
      END

C----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ALPHAS(MUR)
      IMPLICIT NONE
      INTEGER NFF,IVFNS,NF
      DOUBLE PRECISION PI,LOGFR,AS0,M20,ASC,M2C,ASB,M2B,AST,M2T,M2,MUR,
     &     R2,ASI,ASF,R20,R2T,R2B,R2C,AS
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
c
c ..Input common blocks
c
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

       R2 = MUR**2
       M2 = R2 * EXP(+LOGFR)
       IF (IVFNS .EQ. 0) THEN
c
c   Fixed number of flavours
c
          NF  = NFF
          R20 = M20 * R2/M2
          ASI = AS0
          ASF = AS (R2, R20, AS0, NF)
c
      ELSE
c
c ..Variable number of flavours
c
          IF (M2 .GT. M2T) THEN
             NF = 6
             R2T = M2T * R2/M2
             ASI = AST
             ASF = AS (R2, R2T, AST, NF)
c
          ELSE IF (M2 .GT. M2B) THEN
             NF = 5
             R2B = M2B * R2/M2
             ASI = ASB
             ASF = AS (R2, R2B, ASB, NF)
c
          ELSE IF (M2 .GT. M2C) THEN
             NF = 4
             R2C = M2C * R2/M2
             ASI = ASC
             ASF = AS (R2, R2C, ASC, NF)
c
          ELSE
             NF = 3
             R20 = M20 * R2/M2
             ASI = AS0
             ASF = AS (R2, R20, AS0, NF)
c
          END IF
c
       END IF
c
c ..Final value of alpha_s
c
       ALPHAS = 4.D0*PI*ASF

       RETURN
       END
c
c =================================================================av==


c =====================================================================
c
c ..The threshold matching of the QCD coupling in the MS(bar) scheme,
c    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours
c    up to order a_s^4 (NNNLO).
c
c ..The value  ASNF  of a_s for NF flavours at the matching scale, the
c    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
c    of the heavy quark -- and  NF  are passed as arguments to the
c    function  ASNF1.  The order of the expansion  NAORD  (defined as
c    the 'n' in N^nLO) is provided by the common-block  ASPAR.
c
c ..The matching coefficients are inverted from Chetyrkin, Kniehl and
c    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
c    factors have been hard-wired in these results. The lowest integer
c    values of the Zeta function are given by the common-block  RZETA.
c
c =====================================================================
c

      DOUBLE PRECISION FUNCTION ASNF1 (ASNF, LOGRH, NF)
c
      IMPLICIT NONE
      INTEGER NF, NAORD, NASTPS, PRVCLL, K1, K2
      DOUBLE PRECISION ASNF,LOGRH,ZETA,CMC,CMCI30,CMCF30,CMCF31,
     &     CMCI31,ASP,LRHP

      DIMENSION CMC(3,0:3)
c
c ---------------------------------------------------------------------
c
c ..Input common-blocks
c
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / RZETA  / ZETA(6)
c
c ..Variables to be saved for the next call
c
      SAVE CMC, CMCI30, CMCF30, CMCF31, CMCI31, PRVCLL
c
c ---------------------------------------------------------------------
c
c ..The coupling-constant matching coefficients (CMCs) up to NNNLO
c   (calculated and saved in the first call of this routine)
c
       IF (PRVCLL .NE. 1) THEN
c
         CMC(1,0) =  0.D0
         CMC(1,1) =  2./3.D0
c
         CMC(2,0) = 14./3.D0
         CMC(2,1) = 38./3.D0
         CMC(2,2) =  4./9.D0
c
         CMCI30 = + 80507./432.D0 * ZETA(3) + 58933./1944.D0
     1            + 128./3.D0 * ZETA(2) * (1.+ DLOG(2.D0)/3.D0)
         CMCF30 = - 64./9.D0 * (ZETA(2) + 2479./3456.D0)
         CMCI31 =   8941./27.D0
         CMCF31 = - 409./27.D0
         CMC(3,2) = 511./9.D0
         CMC(3,3) = 8./27.D0
c
         PRVCLL = 1
c
       END IF
c
c ---------------------------------------------------------------------
c
c ..The N_f dependent CMCs, and the alpha_s matching at order NAORD
c
       CMC(3,0) = CMCI30 + NF * CMCF30
       CMC(3,1) = CMCI31 + NF * CMCF31
c
       ASNF1 = ASNF
       IF (NAORD .EQ. 0) GO TO 1
       ASP   = ASNF
c
       DO 11 K1 = 1, NAORD
         ASP = ASP * ASNF
         LRHP = 1.D0
c
       DO 12 K2 = 0, K1
         ASNF1 = ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
c
  12   CONTINUE
  11   CONTINUE
c
c ---------------------------------------------------------------------
c
   1   RETURN
       END

c
c =================================================================av==
c
c ..The subroutine  EVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
c    from a three-flavour initial scale to the four- to six-flavour
c    thresholds (identified with the squares of the corresponding quark
c    masses).  The results are written to the common-block  ASFTHR.
c
c ..The input scale  M20 = mu_(f,0)^2  and the corresponding value
c    AS0  of a_s  are provided by  ASINP.  The fixed scale logarithm
c    LOGFR = ln (mu_f^2/mu_r^2) is specified in  FRRAT.  The alpha_s
c    matching is done by the function ASNF1.
c
c =====================================================================
c
c
       SUBROUTINE EVNFTHR (MC2, MB2, MT2)
c
       IMPLICIT NONE
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C,
     1                  R2B, R2T, AS, ASNF1, AS0, ASC, ASB, AST,
     2                  ASC3, ASB4, AST5, LOGFR, SC, SB, ST
c
c ---------------------------------------------------------------------
c
c ..Input common blocks
c
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
c
c ..Output common blocks
c
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

c ---------------------------------------------------------------------
c
c ..Coupling constants at and evolution distances to/between thresholds
c
       R20 = M20 * EXP(-LOGFR)
c
c ..Charm
c
       M2C  = MC2
       R2C  = M2C * R20/M20
       ASC3 = AS (R2C, R20, AS0, 3)
       SC   = LOG (AS0 / ASC3)
       ASC  = ASNF1 (ASC3, -LOGFR, 3)
c
c ..Bottom
c
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = AS (R2B, R2C, ASC, 4)
       SB   = LOG (ASC / ASB4)
       ASB  = ASNF1 (ASB4, -LOGFR, 4)
c
c ..Top
c
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = AS (R2T, R2B, ASB, 5)
       ST   = LOG (ASB / AST5)
       AST  = ASNF1 (AST5, -LOGFR, 5)

       RETURN
       END

c
c =================================================================av==
c
c ..The running coupling of QCD,
c
c         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
c
c    obtained by integrating the evolution equation for a fixed number
c    of massless flavours  NF.  Except at leading order (LO),  AS  is
c    obtained using a fourth-order Runge-Kutta integration.
c
c ..The initial and final scales  R20  and  R2,  the value  AS0  at
c    R20, and  NF  are passed as function arguments.  The coefficients
c    of the beta function up to  a_s^5 (N^3LO)  are provided by the
c    common-block  BETACOM.  The order of the expansion  NAORD (defined
c    as the 'n' in N^nLO) and the number of steps  NASTPS  for the
c    integration beyond LO are given by the common-block  ASPAR.
c
c =====================================================================
c
c
      DOUBLE PRECISION FUNCTION AS (R2, R20, AS0, NF)
c
      IMPLICIT NONE
      INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
      DOUBLE PRECISION R2, R20, AS0, SXTH, BETA0, BETA1, BETA2, BETA3,
     &     FBETA1,FBETA2,FBETA3,A,LRRAT,DLR,XK0,XK1,XK2,XK3
      PARAMETER (NFMIN = 3, NFMAX = 6)
      PARAMETER ( SXTH = 0.16666 66666 66666 D0 )
c
c ---------------------------------------------------------------------
c
c ..Input common-blocks
c
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     ,                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
c
c ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
c
       FBETA1(A) = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * (BETA2(NF) + A * BETA3(NF)) ) )
c
c ---------------------------------------------------------------------
c
c ..Initial value, evolution distance and step size
c
       AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
c
c ..Solution of the evolution equation depending on  NAORD
c   (fourth-order Runge-Kutta beyond the leading order)
c
       IF (NAORD .EQ. 0) THEN
c
         AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
c
       ELSE IF (NAORD .EQ. 1) THEN
c
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS)
         XK1 = DLR * FBETA1 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA1 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA1 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
c
       ELSE IF (NAORD .EQ. 2) THEN
c
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
c
       ELSE IF (NAORD .EQ. 3) THEN
c
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
c
c ---------------------------------------------------------------------
c
       RETURN
       END

c
c =================================================================av==
c
c ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the
c    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by
c
c
c        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ...
c
c    with  a_s = alpha_s/(4*pi).
c
c ..The MSbar coefficients are written to the common-block BETACOM for
c
c   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
c
c ..The factors CF, CA and TF  are taken from the common-block  COLOUR.
c    Beyond NLO the QCD colour factors are hard-wired in this routine,
c    and the numerical coefficients are truncated to six digits.
c
c =====================================================================
c
c
       SUBROUTINE BETAFCT
c
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
c
c ---------------------------------------------------------------------
c
c ..Input common-block
c
       COMMON / COLOUR / CF, CA, TR
c
c ..Output common-block
c
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

c
c ---------------------------------------------------------------------
c
c ..The full LO and NLO coefficients
c
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
c
c ..Flavour-number loop and output to the array
c
       DO 1 NF = NFMIN, NFMAX
c
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
c
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2
     1             + 1.49931 * NF**3
c
c ---------------------------------------------------------------------
c
  1    CONTINUE
c
       RETURN
       END
c
c =================================================================av==


C--   G.W. DZEROX taken from CERNLIB to find the zero of a function.
      DOUBLE PRECISION FUNCTION DZEROX(A0,B0,EPS,MAXF,F,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Based on
C
C        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
C        Guaranteed Convergence for Finding a Zero of a Function,
C        ACM Trans. Math. Software 1 (1975) 330-345.
C
C        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
      CHARACTER*80 ERRTXT
      LOGICAL LMT
      DIMENSION IM1(2),IM2(2),LMT(2)
      PARAMETER (Z1 = 1, HALF = Z1/2)
      DATA IM1 /2,3/, IM2 /-1,3/
      D=0.D0
      FD=0.D0
      DZEROX = 0.D0             ! G.W. to prevent compiler warning
      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
       C=0
       WRITE(ERRTXT,101) MODE
       WRITE(6,*) ERRTXT
       GO TO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=0
       WRITE(ERRTXT,102) A0,B0
       WRITE(6,*) ERRTXT
       GO TO 99
      ENDIF
      ATL=ABS(EPS)
      B=A0
      A=B0
      LMT(2)=.TRUE.
      MF=2
    1 C=A
      FC=FA
    2 IE=0
    3 IF(ABS(FC) .LT. ABS(FB)) THEN
       IF(C .NE. A) THEN
        D=A
        FD=FA
       END IF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      END IF
      TOL=ATL*(1+ABS(C))
      H=HALF*(C+B)
      HB=H-B
      IF(ABS(HB) .GT. TOL) THEN
       IF(IE .GT. IM1(MODE)) THEN
        W=HB
       ELSE
        TOL=TOL*SIGN(Z1,HB)
        P=(B-A)*FB
        LMT(1)=IE .LE. 1
        IF(LMT(MODE)) THEN
         Q=FA-FB
         LMT(2)=.FALSE.
        ELSE
         FDB=(FD-FB)/(D-B)
         FDA=(FD-FA)/(D-A)
         P=FDA*P
         Q=FDB*FA-FDA*FB
        END IF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        END IF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        END IF
       END IF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
        WRITE(6,*) "Error in DZEROX: TOO MANY FUNCTION CALLS"
        GO TO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GO TO 1
       IF(W .EQ. HB) GO TO 2
       IE=IE+1
       GO TO 3
      END IF
      DZEROX=C
   99 CONTINUE
      RETURN
  101 FORMAT('Error in DZEROX: MODE = ',I3,' ILLEGAL')
  102 FORMAT('Error in DZEROX: F(A) AND F(B) HAVE THE SAME SIGN, A = ',
     1     1P,D15.8,', B = ',D15.8)
      END

C ---------------------------------------------------------------------
