      subroutine pA_ll(aintegrand)
      implicit none
      double precision aintegrand
c     =================================================================
c     CepGen common blocks for kinematics definition
c     =================================================================
      include 'cepgen_blocks.inc'

c     =================================================================
c     local variables
c     =================================================================
      double precision s,s12
      double precision alpha1,alpha2,amt1,amt2
      double precision pt1,pt2,eta1,eta2,dely
      double precision pt1x,pt1y,pt2x,pt2y
      double precision ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z
      double precision beta1,beta2,x1,x2
      double precision z1p,z1m,z2p,z2m
      double precision q1tx,q1ty,q1z,q10,q1t2
      double precision q2tx,q2ty,q2z,q20,q2t2
      double precision ptsum
      double precision that1,that2,that,uhat1,uhat2,uhat
      double precision term1,term2,term3,term4,term5,term6,term7
      double precision term8,term9,term10,auxil_gamgam,g_em,amat2
      double precision ak1_x,ak1_y,ak2_x,ak2_y
      double precision eps12,eps22
      double precision aux2_1,aux2_2,f1,f2
      double precision Phi10,Phi102,Phi11_x,Phi11_y,Phi112
      double precision Phi20,Phi202,Phi21_x,Phi21_y,Phi212
      double precision Phi11_dot_e,Phi11_cross_e
      double precision Phi21_dot_e,Phi21_cross_e
      double precision aintegral
      integer imethod,imat1,imat2

      double precision px_plus,px_minus,py_plus,py_minus
      double precision r1,r2
      double precision ptdiffx,ptsumx,ptdiffy,ptsumy
      double precision invm,invm2,s1_eff,s2_eff
      double precision t1abs,t2abs
      double precision inp_A,am_nuc
      double precision p1_plus, p2_minus
      double precision amat2_1,amat2_2
      integer iterm11,iterm22,iterm12,itermtt

c     =================================================================
c     quarks production
c     =================================================================
#ifdef ALPHA_S
      double precision color,e_q,t_max,amu2,a_s,alphas
      double precision rx,rkt2,rmu2,parton
      logical first_init
      data first_init/.false./
      if(first_init.eqv..false.) then
        call f_inter_kmr_fg(rx,rkt2,rmu2,0,parton)
        call initAlphaS(2,1.d0,1.d0,0.5d0,1.4d0,4.75d0,1.d10)
        first_init = .true.
      endif
#endif

c     =================================================================
c     FIXME
c     =================================================================

      aintegrand = 0.d0
      q10 = 0.d0
      q1z = 0.d0
      q20 = 0.d0
      q2z = 0.d0

      call CepGen_print

c     =================================================================
c     go to energy of Nucleus = A*inp in order that generator puts out
c     proper momenta in LAB frame
c     =================================================================

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

c      s = 4.*inp1*inp_A
      s12 = dsqrt(s)

      p1_plus = (ak10+ak1z)/dsqrt(2.d0)
      p2_minus = (ak20-ak2z)/dsqrt(2.d0)

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
c     Outgoing proton final state's mass
c     =================================================================
      if((icontri.eq.1).or.(icontri.eq.2)) am_x = am_p
      if((icontri.eq.1).or.(icontri.eq.3)) am_y = am_p*a_nuc

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

      if(ipt.eq.1) then
        if(pt1.lt.pt_min.or.pt2.lt.pt_min) return
      endif

      amt1 = dsqrt(pt1**2+am_l**2)
      amt2 = dsqrt(pt2**2+am_l**2)

      invm2 = amt1**2 + amt2**2 + 2.d0*amt1*amt2*dcosh(y1-y2)
     >      -ptsum**2

      invm = dsqrt(invm2)

c      if(invm.lt.1500.0d0) return

c     =================================================================
      dely = dabs(y1-y2)
c     =================================================================
c     a window in rapidity distance
c     =================================================================
      if(idely.eq.1) then
        if(dely.lt.dely_min.or.dely.gt.dely_max) return
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

      x1 = alpha1 + alpha2
      x2 = beta1  + beta2

      z1p = alpha1/x1
      z1m = alpha2/x1
      z2p = beta1/x2
      z2m = beta2/x2

c     -----------------------------------------------------------------
      if(x1.gt.1.0.or.x2.gt.1.0) return
c     -----------------------------------------------------------------

      s1_eff = x1*s - q1t**2
      s2_eff = x2*s - q2t**2

c-------------------------------------------------------------------
c     Additional conditions for energy-momentum conservation
c     -----------------------------------------------------------------
      if(((icontri.eq.2).or.(icontri.eq.4))
     1       .and.(dsqrt(s1_eff).le.(am_y+invm))) return
      if(((icontri.eq.3).or.(icontri.eq.4))
     1       .and.(dsqrt(s2_eff).le.(am_x+invm))) return

c     =================================================================
c     >>> TO THE OUTPUT COMMON BLOCK
c     =================================================================

c     =================================================================
c     four-momenta of the outgoing protons (or remnants)
c     =================================================================

      px_plus = (1.d0-x1) * p1_plus
      px_minus = (am_x**2 + q1tx**2 + q1ty**2)/2.d0/px_plus

      px_0 = (px_plus + px_minus)/dsqrt(2.d0)
      px_z = (px_plus - px_minus)/dsqrt(2.d0)
      px_x = - q1tx
      px_y = - q1ty

      am_nuc = a_nuc*am_y
      py_minus = (1.d0-x2) * p2_minus
      py_plus =  (am_nuc**2 + q2tx**2 + q2ty**2)/2.d0/py_minus

      py_0 = (py_plus + py_minus)/dsqrt(2.d0)
      py_z = (py_plus - py_minus)/dsqrt(2.d0)
      py_x = - q2tx
      py_y = - q2ty

      q1t = dsqrt(q1t2)
      q2t = dsqrt(q2t2)

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

      eta1 = 0.5d0*dlog((dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2) +
     2       amt1*dsinh(y1))/(dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2)
     3     - amt1*dsinh(y1)))

      eta2 = 0.5d0*dlog((dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2) +
     2       amt2*dsinh(y2))/(dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2)
     3     - amt2*dsinh(y2)))

      if(ieta.eq.1) then
        if(eta1.lt.eta_min.or.eta1.gt.eta_max) return
        if(eta2.lt.eta_min.or.eta2.gt.eta_max) return
      endif

c     =================================================================
c     matrix element squared
c     averaged over initial spin polarizations
c     and summed over final spin polarizations
c     (--> see Wolfgang's notes
c     =================================================================

c     =================================================================
c     Mendelstam variables
c     =================================================================

      that1 = (q10-p10)**2 -(q1tx-p1x)**2 -(q1ty-p1y)**2 -(q1z-p1z)**2
      uhat1 = (q10-p20)**2 -(q1tx-p2x)**2 -(q1ty-p2y)**2 -(q1z-p2z)**2
      that2 = (q20-p20)**2 -(q2tx-p2x)**2 -(q2ty-p2y)**2 -(q2z-p2z)**2
      uhat2 = (q20-p10)**2 -(q2tx-p1x)**2 -(q2ty-p1y)**2 -(q2z-p1z)**2

      that = (that1+that2)/2.d0
      uhat = (uhat1+uhat2)/2.d0

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
     3             / ( (am_l**2-that)**2 * (am_l**2-uhat)**2 )

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

      t1abs = (q1t2+x1*(am_x**2-am_p**2)+x1**2*am_p**2)/(1.d0-x1)
      t2abs = (q2t2+x2*(am_y**2-am_p**2)+x2**2*am_p**2)/(1.d0-x2)

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

      if(pdg_l.le.6.and.imode.ge.20.and.imode.lt.40) then ! diquarks photo-production
#ifdef ALPHA_S
        color = 0.5d0
        e_Q = 2.d0/3.d0
        t_max = max(amt1**2,amt2**2)
        amu2 = max(eps12,t_max)
        a_S = alphaS(dsqrt(amu2))
        amat2_1 = (4.d0*pi)**2*e_Q**2*alpha_em*a_S*(x1*x2*s)**2
     2          * color* aux2_1 * 2.*z1p*z1m*q1t2 / (q1t2*q2t2)
        amat2_2 = (4.d0*pi)**2*e_Q**2*alpha_em*a_S* (x1*x2*s)**2
     2          * color* aux2_2 * 2.*z2p*z2m*q2t2 / (q1t2*q2t2)
        am_x = amu2
#else
        print *,'alphaS not linked to this instance!'
        stop
#endif
      else ! dilepton production
        amat2_1 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
     2          * aux2_1 * 2.*z1p*z1m*q1t2 / (q1t2*q2t2)
        amat2_2 = (4.d0*pi*alpha_em)**2 * (x1*x2*s)**2
     2          * aux2_2 * 2.*z2p*z2m*q2t2 / (q1t2*q2t2)
       endif

c     =================================================================
c     symmetrization
c     =================================================================

      amat2 = (imat1*amat2_1 + imat2*amat2_2)/2.d0

      endif

c     ============================================
c     unintegrated photon distributions
c     ============================================

      if(icontri.eq.1) then
        f1 = CepGen_kT_flux(imode,q1t2,x1,0,am_x)
        f2 = CepGen_kT_flux_HI(100,q2t2,x2,a_nuc,z_nuc)
      elseif(icontri.eq.2) then
        f1 = CepGen_kT_flux(imode,q1t2,x1,0,am_x)
        f2 = CepGen_kT_flux(imode,q2t2,x2,sfmod,am_y)
      elseif(icontri.eq.3) then
        f1 = CepGen_kT_flux(imode,q1t2,x1,sfmod,am_x)
        f2 = CepGen_kT_flux_HI(100,q2t2,x2,a_nuc,z_nuc)
      elseif(icontri.eq.4) then
        f1 = CepGen_kT_flux(imode,q1t2,x1,sfmod,am_x)
        f2 = CepGen_kT_flux(imode,q2t2,x2,sfmod,am_y)
      endif
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

      return
      end

