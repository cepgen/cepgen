      subroutine INCqqbar(aintegrand,mode,q1t,q2t,phiq1t,phiq2t,
     +     y1,y2,ptdiff,phiptdiff,amx,amy)
      implicit none
c     =================================================================
c     Vegas subroutine:
c
c     different distributions
c
c     for the reaction
c     in p p --> W- W+ reaction
c     in the k_t-factorization approach
c     =================================================================

      integer mode
      integer iread

      double precision s,s12,am_l,am_l2,am_p,pdif
      double precision dely_min,dely_max,pi,alpha_em,units
      double precision ak10,ak1x,ak1y,ak1z,dely,y1,y2,pt1,pt2
      double precision alpha1,alpha2,amt1,amt2,ak20,ak2x,ak2y,ak2z
      double precision sig_ptpt
cc      double precision pt_cutoff
      double precision beta1,beta2,x1,x2,xi_x1,xi_x2
      double precision z1p,z1m,z2p,z2m,qcaptx,qcapty,q1tx,q1ty
      double precision q2tx,q2ty,q10,q1z,q20,q2z,q1x,q1y,q2x,q2y
      double precision akapt12,akapt22
      double precision pt1x,pt1y,pt2x,pt2y
      double precision pcaptx,pcapty
      double precision p10,p1x,p1y,p1z
      double precision p20,p2x,p2y,p2z
      double precision p12,p22
      double precision q1t2,q2t2,q1t,q2t,q12,q22
      double precision ptsum,eta1,eta2,shat,aM342,that1,that2
      double precision stild
      double precision uhat1,uhat2,aMll
      double precision beta
      double precision uhat,that
      double precision term1,term2,term3,term4,term5,term6,term7
      double precision term8,term9,term10,auxil_gamgam,g_em,amat2
      double precision ak1_x,ak1_y,ak2_x,ak2_y,amx,amy
      double precision eps12,eps22,phi10,phi11_x,phi11_y,phi102,phi112
      double precision xcostheta1,xsintheta1,xcostheta2,xsintheta2
      double precision p1mod,p1_x,p1_y,p1_z,p2mod,p2_x,p2_y,p2_z
      double precision phi20,phi21_x,phi21_y,phi202,phi212
      double precision aux2_1,aux2_2,f1,f2
      double precision f_ela,f_ine,f_ine_fit,f_ine_old
      double precision f_ine_Budnev,f_ela_Budnev,f_ela_nucleus
      double precision delta_x1
      double precision Phi11_dot_e,Phi11_cross_e
      double precision Phi21_dot_e,Phi21_cross_e
      double precision aintegral,aintegrand
      double precision ww
      integer icontri,imode,ilepton,imethod,imat1,imat2,idif
      integer ishat,idely

      double precision px_plus,px_minus,px_0,px_x,px_y,px_z
      double precision py_plus,py_minus,py_0,py_x,py_y,py_z
      double precision wgt
      double precision xi_q1t, xi_q2t
      double precision ptdiff, r1,r2
      double precision ptdiffx, ptdiffy
      double precision ptsumx, ptsumy
      double precision invm,invm2,s1_eff,s2_eff
      double precision phiq1t, phiq2t
      double precision sudakov_2,ratio2
      double precision sudakov_1,ratio1,xx1,xx2
      double precision phiptdiff
      double precision t1abs,t2abs,xb1,xb2
      double precision t_max
      double precision a_nuc,inp_A,am_nuc
      double precision p1_plus, p2_minus
      double precision phidiff,phisum
      double precision amat2_0,amat2_1,amat2_interf,amat2_2
      double precision ampli_pp,ampli_mm,ampli_pm,ampli_mp
      double precision fmap_mstw_F2
      double precision fmap_mstw_FL
cc      double precision F2_inter_MSTW
cc      double precision fmap_mstw_F2
      integer iterm11,iterm22,iterm12,itermtt
      integer lam1,lam2,lam3,lam4
c      integer amat2_gamgam_lW_aux, amat2_gamgam_lW
c     =================================================================
c     common blocks (for VegasLL code)
c     =================================================================

      double precision inp1,inp2
      double precision q1t_min,q1t_max,q2t_min,q2t_max
      double precision ptdiff_min,ptdiff_max
      double precision pt_min,pt_max
      double precision y_min,y_max
      double precision amx_min,amx_max,amy_min,amy_max
      double precision ay1,ay2

      common/cuts/inp1,inp2,
     +     ptdiff_min,ptdiff_max,
     +     y_min,y_max,
     +     pt_min,pt_max,
     +     q1t_min,q1t_max,q2t_min,q2t_max,
     +     amx_min,amx_max,amy_min,amy_max

      common/kinematics/pt1x,pt1y,ay1,pt2x,pt2y,ay2,am_l,
     +     ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z,
     +     px_0,px_x,px_y,px_z,py_0,py_x,py_y,py_z

ccc      common/transfer1/fmap_mstw_F2(120,120),
ccc     2                fmap_mstw_FL(120,120)

c     =================================================================

      ay1 = y1
      ay2 = y2

c ===================================================================
c       INPUT:
c       inp1 = proton energy in lab frame
c       inp2 = nucleus energy **per nucleon** in LAB frame
c       Collision is along z-axis
c====================================================================

c       go to energy of Nucleus = A*inp in order that generator puts out
c       proper momenta in LAB frame
c
        a_nuc = 208
        inp_A = inp2*a_nuc


c     =================================================================
c     four-momenta for incoming beams in LAB !!!!
c     =================================================================

        r1 = dsqrt(1.d0-am_p**2/inp1**2)
        r2 = dsqrt(1.d0-(a_nuc*am_p)**2/inp_A**2)


      ak10 = inp1
      ak1x = 0.d0
      ak1y = 0.d0
      ak1z = inp1*r1

      ak20 = inp_A
      ak2x = 0.d0
      ak2y = 0.d0
      ak2z = -inp_A*r2


        s = 4.*inp1*inp_A*(1.d0 + r1*r2)/2.d0
     >                       + am_p**2 + (a_nuc*am_p)**2

c        s = 4.*inp1*inp_A
        s12 = dsqrt(s)

        p1_plus = (ak10+ak1z)/dsqrt(2.d0)
        p2_minus = (ak20-ak2z)/dsqrt(2.d0)


c     =================================================================
c     contribution included
c         icontri = 1: elastic-elastic
c         icontri = 2: elastic-inelastic
c         icontri = 3: inelastic-elastic
c         icontri = 4: inelastic-inelastic
c     =================================================================

      icontri = mode
cc      icontri = 4
cc      icontri = 2
c     =================================================================
c     choice of F2 structure function
c         imode = 1: Szczurek-Uleshchenko
c         imode = 2: Fiore et al.(parametrization of JLAB data)
c         imode = 3: Suri-Yennie
c         imode = 4: ALLM
c         imode = 5: LUX

c     =================================================================

      imode = 5 ! 5

c     =================================================================
c     choice of the quark-antiquark pair
c         ilepton = 1: electrons
c         ilepton = 2: muons
c         ilepton = 3: tauons
c     =================================================================

      ilepton = 2

      if(ilepton.eq.1) then
         am_l = 0.000510998928
      elseif(ilepton.eq.2) then
         am_l = 0.1056583715
      elseif(ilepton.eq.3) then
         am_l = 1.77682
      endif

      am_l2 = am_l*am_l
c     terms in the matrix element
c
      iterm11 = 1         ! LL
      iterm22 = 1         ! TT
      iterm12 = 1         ! LT
      itermtt = 1         ! TT'

c     =================================================================
c     How matrix element is calculated
c         imethod = 0: on-shell formula
c         imethod = 1: off-shell formula
c     =================================================================
        imethod = 1
c     =================================================================
c     two terms in the Wolfgang's formula for
c     off-shell gamma gamma --> l^+ l^-
c     =================================================================
      imat1 = 1
      imat2 = 1
c     =================================================================
c     mass of the proton
c     =================================================================
      am_p = 0.93827203d0
c     =================================================================
c     way of calculating shat
c         ishat = 1: s x1 x2
c         ishat = 2: exact formula
c     =================================================================
      ishat = 2
c     =================================================================
c     extra cuts on the p1t(l) and p2t(l) plane
c         idif = 0: no extra cut
c         idif = 1: extra cut
c     =================================================================
      idif = 0          ! 0 is a standard
      pdif = 2.5d0
c     =================================================================
c     the distance in rapidity between W^+ and W^-
c     =================================================================
      idely = 0         ! 0 or 1
      dely_min = 4.0d0
      dely_max = 5.0d0
c     =================================================================
c     fundamental constants
c     =================================================================
      pi = 4.d0*datan(1.d0)
      alpha_em = 1.d0/137.035d0
c     =================================================================
c     conversion factor
c     1/GeV^2 --> nb
c     =================================================================
      units = 10.d0*(197.3271d0)**2
      stild = s/2.d0*(1+dsqrt(1.d0-(4*am_p**4)/s**2))


c     =================================================================
c     Outgoing proton final state's mass
c     =================================================================
      if((icontri.eq.1).or.(icontri.eq.2)) amx = am_p
      if((icontri.eq.1).or.(icontri.eq.3)) amy = am_p

      q1tx = q1t*cos(phiq1t)
      q1ty = q1t*sin(phiq1t)

      q2tx = q2t*cos(phiq2t)
      q2ty = q2t*sin(phiq2t)

      ptsumx = q1tx+q2tx
      ptsumy = q1ty+q2ty

      ptsum = sqrt(ptsumx**2+ptsumy**2)

      ptdiffx = ptdiff*cos(phiptdiff)
      ptdiffy = ptdiff*sin(phiptdiff)

      pt1x = 0.5*(ptsumx+ptdiffx)
      pt1y = 0.5*(ptsumy+ptdiffy)

      pt2x = 0.5*(ptsumx-ptdiffx)
      pt2y = 0.5*(ptsumy-ptdiffy)

      pt1 = sqrt(pt1x**2+pt1y**2)
      pt2 = sqrt(pt2x**2+pt2y**2)

      if(pt1.lt.pt_min.or.pt2.lt.pt_min) then
        aintegrand = 0.d0
        goto 100
      endif

      amt1 = dsqrt(pt1**2+am_l2)
      amt2 = dsqrt(pt2**2+am_l2)

        invm2 = amt1**2 + amt2**2 + 2.d0*amt1*amt2*dcosh(y1-y2)
     >      -ptsum**2

        invm = dsqrt(invm2)

c        if(invm.lt.1500.0d0) then
c        aintegrand = 0.d0
c        goto 100
c        endif



c	beta = dsqrt(1.d0 - am_l*am_l/invm2)


c     =================================================================
      if(idif.eq.1) then
      if(abs(pt1-pt2).gt.pdif) then
      sig_ptpt = 0.0                                                                         !?????????????????
      goto 100    ! ---->
      endif
      endif
c     =================================================================

      pcaptx = pt1x + pt2x
      pcapty = pt1y + pt2y

      dely = dabs(y1-y2)
c     =================================================================
c     a window in rapidity distance
c     =================================================================
      if(idely.eq.1) then
      if(dely.lt.dely_min.or.dely.gt.dely_max) goto 100
      endif

c     =================================================================
c     auxiliary quantities
c     =================================================================

      alpha1 = amt1/(dsqrt(2.d0)*p1_plus)*dexp( y1)
      alpha2 = amt2/(dsqrt(2.d0)*p1_plus)*dexp( y2)
      beta1  = amt1/(dsqrt(2.d0)*p2_minus)*dexp(-y1)
      beta2  = amt2/(dsqrt(2.d0)*p2_minus)*dexp(-y2)

       q1t2 = q1tx**2 + q1ty**2
      q2t2 = q2tx**2 + q2ty**2

      delta_x1 = (amx**2 + q2t2)/((1.d0-x2)*s)

      x1 = alpha1 + alpha2
      x2 = beta1  + beta2

      xi_x1 = dlog10(x1)
      xi_x2 = dlog10(x2)

      z1p = alpha1/x1
      z1m = alpha2/x1
      z2p = beta1/x2
      z2m = beta2/x2


c     -----------------------------------------------------------------
      if(x1.gt.1.0.or.x2.gt.1.0) then
c      if(x1.gt.0.1.or.x2.gt.0.1) then
        aintegrand=0.d0
        goto 100
        endif
c     -----------------------------------------------------------------


        s1_eff = x1*s - q1t**2
        s2_eff = x2*s - q2t**2

c-------------------------------------------------------------------
c     Additional conditions for energy-momentum conservation
c     -----------------------------------------------------------------
      if(((icontri.eq.2).or.(icontri.eq.4))
     1       .and.(dsqrt(s1_eff).le.(amy+invm))) then
        aintegrand=0.d0
        goto 100
        endif
      if(((icontri.eq.3).or.(icontri.eq.4))
     1       .and.(dsqrt(s2_eff).le.(amx+invm))) then
        aintegrand=0.d0
        goto 100
        endif

c        pt_cutoff = 20.d0
c        if((pt1.gt.pt_cutoff).or.(pt2.gt.pt_cutoff)) then
c                aintegrand = 0.d0
c        goto 100
c        endif


c       if (ptsum.gt.50.d0) then
c           aintegrand=0.d0
c           goto 100
c       endif

c     -----------------------------------------------------------------

c      if(icontri.eq.2) then
c      sudakov_2 = (amy**2 - am_p**2 + q2t2 + x2*am_p**2)
c     >      /((1.d0-x2)*s)
c      sudakov_1 = (q1t2 + x1*am_p**2)/((1.d0-x1)*s)
c      ratio1 = sudakov_1 / x1
c      ratio2 = sudakov_2 / x2
c      if(ratio2.gt.0.1) then
c        aintegrand = 0.d0
c        goto 100
c      endif
c        endif

      qcaptx = pcaptx
      qcapty = pcapty

c     =================================================================
c     four-momenta of the outgoing protons (or remnants)
c     =================================================================

       px_plus = (1.d0-x1) * p1_plus
      px_minus = (amx**2 + q1tx**2 + q1ty**2)/2.d0/px_plus

      px_0 = (px_plus + px_minus)/dsqrt(2.d0)
      px_z = (px_plus - px_minus)/dsqrt(2.d0)
      px_x = - q1tx
      px_y = - q1ty


      am_nuc = a_nuc*amy
      py_minus = (1.d0-x2) * p2_minus
      py_plus =  (am_nuc**2 + q2tx**2 + q2ty**2)/2.d0/py_minus

      py_0 = (py_plus + py_minus)/dsqrt(2.d0)
      py_z = (py_plus - py_minus)/dsqrt(2.d0)
      py_x = - q2tx
      py_y = - q2ty

      q1t = dsqrt(q1t2)
      q2t = dsqrt(q2t2)

      xi_q1t = dlog10(q1t)                   ! new
      xi_q2t = dlog10(q2t)                   ! new



c     =================================================================
c     four-momenta of the outgoing l^+ and l^-
c     =================================================================

      p10 = alpha1*ak10 + beta1*ak20
      p1x = pt1x
      p1y = pt1y
      p1z = alpha1*ak1z + beta1*ak2z

      p20 = alpha2*ak10 + beta2*ak20
      p2x = pt2x
      p2y = pt2y
      p2z = alpha2*ak1z + beta2*ak2z

      p12 = p10**2-p1x**2-p1y**2-p1z**2
      p22 = p20**2-p2x**2-p2y**2-p2z**2

cc      invm = dsqrt((p10+p20)**2-(p1x+p2x)**2-(p1y+p2y)**2-(p1z+p2z)**2)
c     ptsum = dsqrt((p1x+p2x)**2 + (p1y+p2y)**2)
c     =================================================================
c     pseudorapidities of l^+ and l^-
c     =================================================================
c
      eta1 = 0.5d0*dlog((dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2) +
     2       amt1*dsinh(y1))/(dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2)
     3       - amt1*dsinh(y1)))

      eta2 = 0.5d0*dlog((dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2) +
     2       amt2*dsinh(y2))/(dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2)
     3       - amt2*dsinh(y2)))

c
c     matrix element squared
c     averaged over initial spin polarizations
c     and summed over final spin polarizations
c     (--> see Wolfgang's notes
c     =================================================================
c     four-momenta squared of virtual photons
c     =================================================================
      q12 = q10**2-q1tx**2-q1ty**2-q1z**2
      q22 = q20**2-q2tx**2-q2ty**2-q2z**2

ctest      write(*,*) q12,q22
c     =================================================================
c     Mendelstam variables
c     =================================================================
      if(ishat.eq.1) then
      shat = s*x1*x2
      elseif(ishat.eq.2) then
      shat = (q10+q20)**2-(q1tx+q2tx)**2-(q1ty+q2ty)**2-(q1z+q2z)**2
      endif

      aM342 = shat

      that1 = (q10-p10)**2 -(q1tx-p1x)**2 -(q1ty-p1y)**2 -(q1z-p1z)**2
      uhat1 = (q10-p20)**2 -(q1tx-p2x)**2 -(q1ty-p2y)**2 -(q1z-p2z)**2
      that2 = (q20-p20)**2 -(q2tx-p2x)**2 -(q2ty-p2y)**2 -(q2z-p2z)**2
      uhat2 = (q20-p10)**2 -(q2tx-p1x)**2 -(q2ty-p1y)**2 -(q2z-p1z)**2

      aMll = dsqrt(aM342)

      that = (that1+that2)/2.d0
      uhat = (uhat1+uhat2)/2.d0
c     =================================================================
c     polar angles of l^+ and l^-
c     =================================================================

      p1mod = dsqrt(p1_x**2+p1_y**2+p1_z**2)
      p2mod = dsqrt(p2_x**2+p2_y**2+p2_z**2)

      xcostheta1 = p1_z / p1mod
      xsintheta1 = dsqrt(p1_x**2+p1_y**2) / p1mod

      xcostheta2 = p2_z / p2mod
      xsintheta2 = dsqrt(p2_x**2+p2_y**2) / p2mod

c     =================================================================
c     matrix elements
c     =================================================================
      if(imethod.eq.0) then
c     =================================================================
c     on-shell formula for M^2
c     =================================================================
      term1 = 6.d0*am_l**8
      term2 = -3.d0*am_l**4*that**2
      term3 = -14.d0*am_l**4*that*uhat
      term4 = -3.d0*am_l**4*uhat**2
      term5 = am_l**2*that**3
      term6 = 7.d0*am_l**2*that**2*uhat
      term7 = 7.d0*am_l**2*that*uhat**2
      term8 = am_l**2*uhat**3
      term9  = -that**3*uhat
      term10 = -that*uhat**3

      auxil_gamgam = -2.d0*(  term1+term2+term3+term4+term5
     2                    +term6+term7+term8+term9+term10 )
     3             / ( (am_l2-that)**2 * (am_l2-uhat)**2 )

      g_em = dsqrt(4.d0*pi*alpha_em)

      amat2 = g_em**4*auxil_gamgam

      elseif(imethod.eq.1)then
c     =================================================================
c     Wolfgang's formulae
c     =================================================================

      ak1_x = z1m*pt1x-z1p*pt2x
      ak1_y = z1m*pt1y-z1p*pt2y

      ak2_x = z2m*pt1x-z2p*pt2x
      ak2_y = z2m*pt1y-z2p*pt2y

      t1abs = (q1t2+x1*(amx**2-am_p**2)+x1**2*am_p**2)/(1.d0-x1)
      t2abs = (q2t2+x2*(amy**2-am_p**2)+x2**2*am_p**2)/(1.d0-x2)

c      t1abs = (q1t2+x1**2*am_p**2)/(1.d0-x1)
c      t2abs = (q2t2+x2**2*am_p**2)/(1.d0-x2)

	t_max = max(pt1**2,pt2**2)

cc      if(t1abs.gt.t_max.or.t2abs.gt.t_max) then
cc        aintegrand = 0.d0
cc        goto 100
cc      endif



c      eps12 = am_l**2 + z1p*z1m*q1t2
c      eps22 = am_l**2 + z2p*z2m*q2t2
      eps12 = am_l**2 + z1p*z1m*t1abs
      eps22 = am_l**2 + z2p*z2m*t2abs

      Phi10 = 1.d0/((ak1_x+z1p*q2tx)**2+(ak1_y+z1p*q2ty)**2+eps12)
     2      - 1.d0/((ak1_x-z1m*q2tx)**2+(ak1_y-z1m*q2ty)**2+eps12)
      Phi11_x = (ak1_x+z1p*q2tx)/
     2          ((ak1_x+z1p*q2tx)**2+(ak1_y+z1p*q2ty)**2+eps12)
     3        - (ak1_x-z1m*q2tx)/
     4          ((ak1_x-z1m*q2tx)**2+(ak1_y-z1m*q2ty)**2+eps12)
      Phi11_y = (ak1_y+z1p*q2ty)/
     2          ((ak1_x+z1p*q2tx)**2+(ak1_y+z1p*q2ty)**2+eps12)
     3        - (ak1_y-z1m*q2ty)/
     4          ((ak1_x-z1m*q2tx)**2+(ak1_y-z1m*q2ty)**2+eps12)

      Phi102 = Phi10*Phi10
      Phi112 = Phi11_x**2+Phi11_y**2

      Phi20 = 1.d0/((ak2_x+z2p*q1tx)**2+(ak2_y+z2p*q1ty)**2+eps22)
     2      - 1.d0/((ak2_x-z2m*q1tx)**2+(ak2_y-z2m*q1ty)**2+eps22)

      Phi21_x = (ak2_x+z2p*q1tx)/
     2          ((ak2_x+z2p*q1tx)**2+(ak2_y+z2p*q1ty)**2+eps22)
     3        - (ak2_x-z2m*q1tx)/
     4          ((ak2_x-z2m*q1tx)**2+(ak2_y-z2m*q1ty)**2+eps22)
      Phi21_y = (ak2_y+z2p*q1ty)/
     2          ((ak2_x+z2p*q1tx)**2+(ak2_y+z2p*q1ty)**2+eps22)
     3        - (ak2_y-z2m*q1ty)/
     4          ((ak2_x-z2m*q1tx)**2+(ak2_y-z2m*q1ty)**2+eps22)

      Phi202 = Phi20*Phi20
      Phi212 = Phi21_x**2+Phi21_y**2

      Phi11_dot_e = (Phi11_x*q1tx + Phi11_y*q1ty)/dsqrt(q1t2)
      Phi11_cross_e = (Phi11_x*q1ty - Phi11_y*q1tx)/dsqrt(q1t2)

      Phi21_dot_e = (Phi21_x*q2tx +Phi21_y*q2ty)/dsqrt(q2t2)
      Phi21_cross_e = (Phi21_x*q2ty -Phi21_y*q2tx)/dsqrt(q2t2)

      aux2_1 = iterm11*(am_l**2+4.d0*z1p**2*z1m**2*t1abs)*Phi102
     1      +iterm22*( (z1p**2 + z1m**2)*(Phi11_dot_e**2 +
     2      Phi11_cross_e**2) )
     3      + itermtt*( Phi11_cross_e**2 - Phi11_dot_e**2)
     4      - iterm12*4.d0*z1p*z1m*(z1p-z1m)*Phi10
     5      *(q1tx*Phi11_x+q1ty*Phi11_y)

      aux2_2 = iterm11*(am_l**2+4.d0*z2p**2*z2m**2*t2abs)*Phi202
     1     +iterm22*( (z2p**2 + z2m**2)*(Phi21_dot_e**2 +
     2     Phi21_cross_e**2) )
     3     + itermtt*( Phi21_cross_e**2 - Phi21_dot_e**2)
     4     - iterm12*4.d0*z2p*z2m*(z2p-z2m)*Phi20
     5     *(q2tx*Phi21_x+q2ty*Phi21_y)


c     =================================================================
c     convention of matrix element as in our kt-factorization
c     for heavy flavours
c     =================================================================
      amat2_1 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
     2        * aux2_1 * 2.*z1p*z1m*q1t2 / (q1t2*q2t2)
      amat2_2 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
     2        * aux2_2 * 2.*z2p*z2m*q2t2 / (q1t2*q2t2)


cc      amat2_1 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
cc     2        * aux2_1 * 2.*z1p*z1m*t1abs / (q1t2*q2t2)*t2abs/q2t2

c      amat2_1 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
c     2        * aux2_1 * 2.*z1p*z1m*t1abs / (q1t2*q2t2)

cc      amat2_2 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
cc     2        * aux2_2 * 2.*z2p*z2m*t2abs / (q1t2*q2t2)*t1abs/q1t2

c
ctest      write(*,*) amat2_1, amat2_2
c     =================================================================
c     symmetrization
c     =================================================================

      amat2 = (imat1*amat2_1 + imat2*amat2_2)/2.d0

      xx1 = alpha1 + alpha2
      xx2 = beta1  + beta2


      sudakov_2 = (amx**2 - am_p**2 + q2t2 + xx2*am_p**2)/((1.d0-xx2)*s)
      sudakov_1 = (q1t2 + xx1*am_p**2)/((1.d0-xx1)*s)
      ratio1 = sudakov_1 / xx1
      ratio2 = sudakov_2 / xx2

c      if(ratio1.gt.0.01) then
c        aintegrand = 0.d0
c        goto 100
c      endif

      endif

c     ============================================
c     unintegrated photon distributions
c     interpolation on double logarithmic grid
c     of inelastic distributions
c     ============================================

      if(icontri.eq.1) then
      f1 = f_ela_Budnev(x1,q1t2)
      f2 = f_ela_nucleus(x2,q2t2)
      elseif(icontri.eq.2) then
      f1 = f_ela_Budnev(x1,q1t2)
      if(imode.eq.1) then
c      f2 = ff_int_xixxikt(imode,x2,q2t2,1)
      f2 = f_ine(x2,q2t2,amy)
c      f2 = f_ine_old(x2,q2t2,amy)
      else
      f2 = f_ine_Budnev(imode,x2,q2t2,amy)
c      f2 = f_ine_MRST(x2,q2t2,amy)
      endif
      elseif(icontri.eq.3) then
      if(imode.eq.1) then
      f1 = f_ine(x1,q1t2,amx)
cc      f1 = f_ine_old(x1,q1t2,amx)
      else
      f1 = f_ine_Budnev(imode,x1,q1t2,amx)
cc      f1 = f_ine_MRST(x1,q1t2,amx)
      endif
      f2 = f_ela_nucleus(x2,q2t2)
c      f1 = ff_int_xixxikt(imode,x1,q1t2,1)
c      f2 = f_ela(x2,q2t2)
       elseif(icontri.eq.4) then
      if(imode.eq.1) then
      f1 = f_ine(x1,q1t2,amx)
      f2 = f_ine(x2,q2t2,amy)
cc      f1 = f_ine_old(x1,q1t2,amx)
c        f1 = 1.d0
cc      f2 = f_ine_old(x2,q2t2,amy)
c        f2 = 1.d0
      else
      f1 = f_ine_Budnev(imode,x1,q1t2,amx)
c      f1 = f_ine_fit(imode,x1,q1t2,amx)
c      f2 = f_ine_fit(imode,x2,q2t2,amy)
      f2 = f_ine_Budnev(imode,x2,q2t2,amy)
cc      f1 = f_ine_MRST(x1,q1t2,amx)
cc      f2 = f_ine_MRST(x2,q2t2,amy)


      endif
c      f1 = ff_int_xixxikt(imode,x1,q1t2,1)
c      f2 = ff_int_xixxikt(imode,x2,q2t2,1)
      endif

ctest      write(*,*) x1,q1t2,f1,x2,q2t2,f2

      if(f1.lt.1.d-20) f1 = 0.0d0
      if(f2.lt.1.d-20) f2 = 0.0d0
c     =================================================================
c     factor 2.*pi below from integration over phi_sum
c     factor 1/4 below from jacobian of transformations
c     factors 1/pi and 1/pi due to integration
c     over d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
c     =================================================================

      aintegral = (2.d0*pi)*1.d0/(16.d0*pi**2*(x1*x2*s)**2) * amat2
     2          * f1/pi * f2/pi * (1.d0/4.d0) * units
     3          * 0.5d0 * 4.0d0 / (4.d0*pi)

c     *****************************************************************
c     =================================================================
      aintegrand = aintegral*q1t*q2t*ptdiff
c     =================================================================
c     *****************************************************************

c     -----------------------------------------------------------------
  100 continue
c     -----------------------------------------------------------------
      return
      end

c     =================================================================
c     =================================================================
      function f_ela_nucleus(x,akt2)
c     =================================================================
      double precision f_ela_nucleus,x,akt2
      double precision tau,q2_ela

        A = 208
        Z = 82
        RA = 1.1*A**(1./3.)
        R_ch = 5.5 ! Lead 208
        a0 = 0.7



cold      const = 0.004d0
cold      f_ela = const*(1.d0-x)/akt2**2

      am_p = 0.93827203d0
        am_A = A*am_p

      alpha_em = 1.d0/137.035d0
      pi = 4.d0*atan(1.d0)
c      am_p = 0.93827d0

      Q2_ela = (akt2 + x**2*am_A**2)/(1.d0-x)
        tau = dsqrt(Q2_ela)*RA/.1973
        tau1=dsqrt(Q2_ela)*a0/.1973
c      ''Realistic nuclear formfactor'' as used in STARLIGHT
        ff1 = 3.*(dsin(tau)-tau*dcos(tau))/(tau + 1.d-10)**3
        ff2 = 1.d0/(1.d0 + tau1**2)

      ela1 = (akt2/(akt2+x**2*am_A**2))**2
c      ela2 = (4.d0*am_p**2*G_E**2 + Q2_ela*G_M**2)/(4.d0*am_p**2+Q2_ela)
        ela2 = (ff1*ff2)**2
      ela3 = 1.d0-(Q2_ela-akt2)/Q2_ela
c        ela2 = 1.d0
c        ela3 = 1.d0 - x**2*am_p**2/Q2_ela/(1.d0-x)
c        ela3 = 1.d0
c        f_ela = alpha_em/pi*(1.d0-x+x**2/4.d0)*ela1*ela2*ela3 / akt2
      f_ela = Z**2*alpha_em/pi*ela1*ela2/Q2_ela
c        f_ela_Nucleus = alpha_em/pi*((1.d0-x)*ela1*ela2*ela3 + x**2/2. *G_M**2)
c     >     / akt2

      return
      end

