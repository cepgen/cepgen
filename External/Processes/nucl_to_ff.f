      function nucl_to_ff()
      implicit none
      double precision nucl_to_ff
c     =================================================================
c     CepGen common blocks for kinematics definition
c     =================================================================
      include 'KTBlocks.inc'
      data iflux1,iflux2,sfmod,pdg_l/10,100,11,13/
      data a_nuc1,z_nuc1,a_nuc2,z_nuc2/1,1,208,82/

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
      double precision term8,term9,term10,amat2
      double precision ak1_x,ak1_y,ak2_x,ak2_y
      double precision eps12,eps22
      double precision aux2_1,aux2_2,f1,f2
      double precision Phi10,Phi102,Phi11_x,Phi11_y,Phi112
      double precision Phi20,Phi202,Phi21_x,Phi21_y,Phi212
      double precision Phi11_dot_e,Phi11_cross_e
      double precision Phi21_dot_e,Phi21_cross_e
      double precision aintegral
      integer imethod,pdg_l,imat1,imat2

      double precision px_plus,px_minus,py_plus,py_minus
      double precision r1,r2
      double precision ptdiffx,ptsumx,ptdiffy,ptsumy
      double precision invm,invm2,s1_eff,s2_eff
      double precision t1abs,t2abs
      double precision am_l,q_l
      double precision inp_A,inp_B,am_A,am_B
      double precision p1_plus, p2_minus
      double precision amat2_1,amat2_2
      integer iterm11,iterm22,iterm12,itermtt

      double precision coupling

c     =================================================================
c     quarks production
c     =================================================================
#ifdef ALPHA_S
      double precision t_max,amu2,alphas
#endif
      logical first_init
      data first_init/.true./
      save first_init,imethod,pdg_l,am_l,q_l

      if(first_init) then
        call CepGen_set_process('nucl_to_ff')
        call CepGen_print
        imethod = CepGen_param_int('method', 1)
        pdg_l = CepGen_param_int('pair', 13)
        am_l = CepGen_particle_mass(pdg_l) ! central particles mass
        q_l = CepGen_particle_charge(pdg_l) ! central particles charge
        if(iflux1.ge.20.and.iflux1.lt.40) then
          if(icontri.eq.3.or.icontri.eq.4) then
            print *,'Invalid process mode for collinear gluon emission!'
            stop
          endif
#ifdef ALPHA_S
          print *,'Initialisation of the alpha(S) evolution algorithm..'
          call initAlphaS(0,1.d0,1.d0,0.5d0,
     &       CepGen_particle_mass(4), ! charm
     &       CepGen_particle_mass(5), ! bottom
     &       CepGen_particle_mass(6)) ! top
#endif
        endif
        first_init = .false.
      endif

c     =================================================================
c     FIXME
c     =================================================================

      nucl_to_ff = 0.d0
      amat2 = 0.d0
      eps12 = 0.d0
      eps22 = 0.d0
      q10 = 0.d0
      q1z = 0.d0
      q20 = 0.d0
      q2z = 0.d0

c     =================================================================
c     go to energy of Nucleus = A*inp in order that generator puts out
c     proper momenta in LAB frame
c     =================================================================

      inp_A = inp1*a_nuc1
      am_A = am_p*a_nuc1
      inp_B = inp2*a_nuc2
      am_B = am_p*a_nuc2

c     =================================================================
c     four-momenta for incoming beams in LAB !!!!
c     =================================================================

      r1 = dsqrt(1.d0+am_A**2/inp_A**2)
      r2 = dsqrt(1.d0+am_B**2/inp_B**2)

      ak10 = inp_A*r1
      ak1x = 0.d0
      ak1y = 0.d0
      ak1z = inp_A

      ak20 = inp_B*r2
      ak2x = 0.d0
      ak2y = 0.d0
      ak2z = -inp_B

      s = 4.*inp_A*inp_B*(1.d0 + r1*r2)/2.d0+(am_A**2+am_B**2)
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
c     two terms in the Wolfgang's formula for
c     off-shell gamma gamma --> l^+ l^-
c     =================================================================
      imat1 = 1
      imat2 = 1
c     =================================================================
c     Outgoing proton final state's mass
c     =================================================================
      if((icontri.eq.1).or.(icontri.eq.2)) am_x = am_A
      if((icontri.eq.1).or.(icontri.eq.3)) am_y = am_B

      q1tx = q1t*cos(phiq1t)
      q1ty = q1t*sin(phiq1t)

      q2tx = q2t*cos(phiq2t)
      q2ty = q2t*sin(phiq2t)

      ptsumx = q1tx+q2tx
      ptsumy = q1ty+q2ty

      ptsum = sqrt(ptsumx**2+ptsumy**2)

c     =================================================================
c     a window in final state transverse momentum
c     =================================================================

      if(iptsum) then
        if(ptsum.lt.ptsum_min.or.ptsum.gt.ptsum_max) return
      endif

c     =================================================================
c     compute the individual central particles momentum
c     =================================================================

      ptdiffx = ptdiff*cos(phiptdiff)
      ptdiffy = ptdiff*sin(phiptdiff)

      pt1x = 0.5*(ptsumx+ptdiffx)
      pt1y = 0.5*(ptsumy+ptdiffy)

      pt2x = 0.5*(ptsumx-ptdiffx)
      pt2y = 0.5*(ptsumy-ptdiffy)

      pt1 = sqrt(pt1x**2+pt1y**2)
      pt2 = sqrt(pt2x**2+pt2y**2)

      if(ipt) then
        if(pt1.lt.pt_min.or.pt2.lt.pt_min) return
        if(pt_max.gt.0d0.and.(pt2.gt.pt_max.or.pt2.gt.pt_max)) return
      endif

      amt1 = dsqrt(pt1**2+am_l**2)
      amt2 = dsqrt(pt2**2+am_l**2)

      invm2 = amt1**2 + amt2**2 + 2.d0*amt1*amt2*dcosh(y1-y2)
     &       -ptsum**2

      invm = dsqrt(invm2)

c     =================================================================
c     a window in final state invariant mass
c     =================================================================

      if(iinvm) then
        if(invm.lt.invm_min) return
        if(invm_max.gt.0d0.and.invm.gt.invm_max) return
      endif

c     =================================================================
c     a window in rapidity distance
c     =================================================================

      dely = dabs(y1-y2)
      if(idely) then
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

c     -----------------------------------------------------------------
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
      px_minus = ((a_nuc1*am_x)**2 + q1tx**2 + q1ty**2)/2.d0/px_plus

      px(1) = - q1tx
      px(2) = - q1ty
      px(3) = (px_plus - px_minus)/dsqrt(2.d0)
      px(4) = (px_plus + px_minus)/dsqrt(2.d0)

      py_minus = (1.d0-x2) * p2_minus
      py_plus =  ((a_nuc2*am_y)**2 + q2tx**2 + q2ty**2)/2.d0/py_minus

      py(1) = - q2tx
      py(2) = - q2ty
      py(3) = (py_plus - py_minus)/dsqrt(2.d0)
      py(4) = (py_plus + py_minus)/dsqrt(2.d0)

      q1t = dsqrt(q1t2)
      q2t = dsqrt(q2t2)

c     =================================================================
c     four-momenta of the outgoing central particles
c     =================================================================

      nout = 2

      ipdg(1) = pdg_l
      pc(1,1) = pt1x
      pc(1,2) = pt1y
      pc(1,3) = alpha1*ak1z + beta1*ak2z
      pc(1,4) = alpha1*ak10 + beta1*ak20

      ipdg(2) = -pdg_l
      pc(2,1) = pt2x
      pc(2,2) = pt2y
      pc(2,3) = alpha2*ak1z + beta2*ak2z
      pc(2,4) = alpha2*ak10 + beta2*ak20

      eta1 = 0.5d0*dlog((dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2) +
     2       amt1*dsinh(y1))/(dsqrt(amt1**2*(dcosh(y1))**2 - am_l**2)
     3     - amt1*dsinh(y1)))

      eta2 = 0.5d0*dlog((dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2) +
     2       amt2*dsinh(y2))/(dsqrt(amt2**2*(dcosh(y2))**2 - am_l**2)
     3     - amt2*dsinh(y2)))

      if(ieta) then
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

      that1 = (q10-pc(1,4))**2
     &       -(q1tx-pc(1,1))**2-(q1ty-pc(2,1))**2-(q1z-pc(3,1))**2
      uhat1 = (q10-pc(2,4))**2
     &       -(q1tx-pc(1,2))**2-(q1ty-pc(2,2))**2-(q1z-pc(3,2))**2
      that2 = (q20-pc(2,4))**2
     &       -(q2tx-pc(1,2))**2-(q2ty-pc(2,2))**2-(q2z-pc(3,2))**2
      uhat2 = (q20-pc(1,4))**2
     &       -(q2tx-pc(1,1))**2-(q2ty-pc(2,1))**2-(q2z-pc(3,1))**2

      that = (that1+that2)/2.d0
      uhat = (uhat1+uhat2)/2.d0

c     =================================================================
c     matrix elements
c     =================================================================
c     How matrix element is calculated
c         imethod = 0: on-shell formula
c         imethod = 1: off-shell formula
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

      amat2 = -2.d0*(  term1+term2+term3+term4+term5
     2                    +term6+term7+term8+term9+term10 )
     3             / ( (am_l**2-that)**2 * (am_l**2-uhat)**2 )

      elseif(imethod.eq.1)then
c     =================================================================
c     Wolfgang's formulae
c     =================================================================

      ak1_x = z1m*pt1x-z1p*pt2x
      ak1_y = z1m*pt1y-z1p*pt2y

      ak2_x = z2m*pt1x-z2p*pt2x
      ak2_y = z2m*pt1y-z2p*pt2y

      !FIXME FIXME FIXME am_p or am_A/B???
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

      amat2_1 = (x1*x2*s)**2 * aux2_1 * 2.*z1p*z1m*q1t2 / (q1t2*q2t2)
      amat2_2 = (x1*x2*s)**2 * aux2_2 * 2.*z2p*z2m*q2t2 / (q1t2*q2t2)

c     =================================================================
c     symmetrization
c     =================================================================

      amat2 = (imat1*amat2_1 + imat2*amat2_2)/2.d0

      endif

      coupling = 1.d0
c     =================================================================
c     first parton coupling
c     =================================================================
      if(iflux1.ge.20.and.iflux1.lt.40) then ! at least one gluon exchanged
#ifdef ALPHA_S
        t_max = max(amt1**2,amt2**2)
        amu2 = max(eps12,t_max)
        am_x = dsqrt(amu2)
        coupling = coupling * 4.d0*pi*alphaS(am_x)/2.d0 ! colour flow
#else
        print *,'alphaS not linked to this instance!'
        stop
#endif
      else ! photon exchange
        coupling = coupling * 4.d0*pi*alpha_em*q_l**2
      endif
c     =================================================================
c     second parton coupling
c     =================================================================
      coupling = coupling * 4.d0*pi*alpha_em*q_l**2 ! photon exchange
      coupling = coupling * 3.d0

c     ============================================
c     unintegrated parton distributions
c     ============================================

      if(a_nuc1.le.1) then
        f1 = CepGen_kT_flux(iflux1,x1,q1t2,sfmod,am_x)
      else
        f1 = CepGen_kT_flux_HI(iflux1,x1,q1t2,a_nuc1,z_nuc1)
      endif
      if(a_nuc2.le.1) then
        f2 = CepGen_kT_flux(iflux2,x2,q2t2,sfmod,am_y)
      else
        f2 = CepGen_kT_flux_HI(iflux2,x2,q2t2,a_nuc2,z_nuc2)
      endif

c     =================================================================
c     factor 2.*pi below from integration over phi_sum
c     factor 1/4 below from jacobian of transformations
c     factors 1/pi and 1/pi due to integration
c     over d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
c     =================================================================

      aintegral = (2.d0*pi)*1.d0/(16.d0*pi**2*(x1*x2*s)**2) * amat2
     &          * coupling
     &          * f1/pi * f2/pi * (1.d0/4.d0) * units
     &          * 0.5d0 * 4.0d0 / (4.d0*pi)

c     *****************************************************************
c     =================================================================
      nucl_to_ff = aintegral*q1t*q2t*ptdiff
c     =================================================================
c     *****************************************************************
c      print *,nucl_to_ff,aintegral,coupling

      return
      end
