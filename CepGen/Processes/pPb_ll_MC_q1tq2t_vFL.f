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

c      include 'ff_int_xixxikt.f'
c      include 'grv_lo.f'



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




c
c     =================================================================
c
      function f_ela_Budnev(x,akt2)

      double precision f_ela,x,akt2,f_ela_Budnev

cold      const = 0.004d0
cold      f_ela = const*(1.d0-x)/akt2**2

      am_p = 0.93827203d0

      alpha_em = 1.d0/137.035d0
      pi = 4.d0*atan(1.d0)
c      am_p = 0.93827d0

      Q2_ela = (akt2 + x**2*am_p**2)/(1.d0-x)
      Q2_min = x**2*am_p**2/(1.d0-x)
      G_dip = 1.d0/(1.d0+Q2_ela/0.71d0)**2
      G_E = G_dip      
      G_M = 2.79d0*G_dip 
      ela1 = (1.d0-x)*(1.d0-Q2_min/Q2_ela)
      ela2 = (4.d0*am_p**2*G_E**2 + Q2_ela*G_M**2)/(4.d0*am_p**2+Q2_ela)
      ela3 = 1.d0-(Q2_ela-akt2)/Q2_ela

      f_ela = alpha_em/pi*(ela1*ela2 + x**2/2.d0*G_M**2)/akt2
c	last factor below the jacobian from dQ^2/Q^2 --> dkT^2/kT^2*(kT^2/Q^2)

      f_ela_Budnev = f_ela*(1.d0-x)*akt2/Q2_ela

      return
      end
c
c     ====================================================================

c
c     =================================================================
c
      function f_ela(x,akt2)
      double precision f_ela,x,akt2

cold      const = 0.004d0
cold      f_ela = const*(1.d0-x)/akt2**2

      am_p = 0.93827203d0

      alpha_em = 1.d0/137.035d0
      pi = 4.d0*atan(1.d0)
c      am_p = 0.93827d0

      Q2_ela = (akt2 + x**2*am_p**2)/(1.d0-x)
      G_dip = 1.d0/(1.d0+Q2_ela/0.71d0)**2
      G_E = G_dip      
      G_M = 2.79d0*G_dip 
      ela1 = (akt2/(akt2+x**2*am_p**2))**2
      ela2 = (4.d0*am_p**2*G_E**2 + Q2_ela*G_M**2)/(4.d0*am_p**2+Q2_ela)
      ela3 = 1.d0-(Q2_ela-akt2)/Q2_ela

cc	f_Budnev = alpha_em/pi*(sqrt(ela1)*ela2 + x**2/2.d0*G_M**2)/akt2
c        ela2 = 1.d0
c        ela3 = 1.d0 - x**2*am_p**2/Q2_ela/(1.d0-x)
c        ela3 = 1.d0
cc        f_ela = alpha_em/pi*((1.d0-x+x**2/4.d0)*ela1*ela2*ela3
cc     >	+ x**2/2. *G_M**2)/akt2
cc      f_ela = alpha_em/pi*((1.d0-x+x**2/4.d0)*ela1*ela2*ela3)/akt2
        f_ela = alpha_em/pi*(1.d0-x)*ela1*ela2*ela3 / akt2
c       f_ela = alpha_em/pi*ela1*ela2/Q2_ela
cc        f_ela = alpha_em/pi*((1.d0-x)*ela1*ela2*ela3 + x**2/2. *G_M**2)
cc     >     / akt2

cc	f_ela = f_Budnev

      return
      end
c
c     ====================================================================
c
      function f_ine(x,akt2,amx)

      double precision f_ine,x,akt2,amx

      alpha_em = 1./137.035
      pi = 4.*atan(1.)
      am_p = 0.93827203
      am_pi = 0.1349766       ! mass of pi^0
c
      aMX2 = aMX*aMX
c
c     F2 structure function
c
      Q2min = 1./(1.-x)*(x*(aMX2-am_p**2)+x**2*am_p**2)
      Q2 = akt2/(1.d0-x) + Q2min
      x_Bj = Q2 / (Q2 + aMX2 - am_p**2)
cc      print *,'q2=',q2
cc      print *,'x_Bj=',x_Bj

      Q02 = 0.8d0

      amu2 = Q2 + Q02       ! scale is shifted


      if (amu2.gt.1d6) then
cc             write(*,*) 'amu2=',amu2 
              f2_aux = 0.d0
      else

      call grv95lo(x_Bj,amu2,xuv,xdv,xus,xds,xss,xg)

      F2_aux = 4./9.*(xuv + 2.*xus)
     2       + 1./9.*(xdv + 2.*xds)
     3       + 1./9.*2.*xss

      endif
c
c     F2 corrected for low Q^2 behaviour
c
      F2_corr = Q2 / (Q2 + Q02) * F2_aux

       term1 = (1.- x/2.d0 * (aMX2 - am_p**2 + Q2)/Q2)**2
c      term1 = (1.- x * (aMX2 - am_p**2 + Q2)/Q2)
c      term1 = (1.-(Q2-akt2)/Q2)
c      term1 = (1.- Q2min/Q2)
c       term1 = 1.d0
        term2 = (akt2/(akt2 + x*(aMX2 - am_p**2) + x**2*am_p**2))**2

      f_aux = F2_corr/(aMX2 + Q2 - am_p**2) * term1 * term2

      f_ine = alpha_em/pi*(1.-x)*f_aux/akt2

c      if(x.gt.x_Bj) then
c              f_ine= 0.d0
c      endif
c        if (Q2.gt.1.e5.or.x_Bj.gt.0.1) then
c                f_ine = 0.d0
c        endif

c        if (Q2.gt.1.e5) then
c                f_ine = 0.d0
c        endif

c        write(*,*) x_Bj,amu2,f_aux,f_ine

      return
      end

c
c     ====================================================================
c
c
c     ====================================================================
c
      function f_ine_Budnev(imode,x,akt2,amx)

c    implicit double precision (a-h,o-z)
      double precision f_ine_Budnev,x,x_Bj,akt2,amx,q2min,q2,tmp
      double precision FL

c      integer mode
      alpha_em = 1./137.035
      pi = 4.*atan(1.)
      am_p = 0.93827203
      am_pi = 0.1349766       ! mass of pi^0
c
      aMX2 = aMX*aMX
c
c     F2 structure function
c
      Q2min = 1./(1.-x)*(x*(aMX2-am_p**2)+x**2*am_p**2)
      Q2 = akt2/(1.d0-x) + Q2min
      x_Bj = Q2 / (Q2 + aMX2 - am_p**2)

cc    write(*,*) x,Q2,x_Bj


      if (imode.eq.2) then
      call F2_Fiore(x_Bj,Q2,tmp)
      elseif(imode.eq.3) then
c      call F2_SY(x_Bj,Q2,tmp)
        tmp = 0.d0
      elseif(imode.eq.4) then
      call F2_ALLM(x_Bj,Q2,tmp)
      elseif(imode.eq.5) then
      call F2_fit_luxlike(x_Bj,Q2,tmp,FL)
      endif

c    RL from Sibirtsev & Blunden Phys Rev C 88,065202 (2013)

      RL=.014*Q2*(dexp(-0.07d0*Q2)+41.d0*dexp(-0.8d0*Q2))
      F2=tmp
        !F2 = 1.0
c    for imode = 2,3,4 calculate F1 from Sibirtsev parametrization of RL

      F1 = (1.d0+4*x_Bj**2*am_p**2/Q2)*F2/(1.d0+RL)/(2.d0*x_Bj)

c    for imode = 5 (luxlike) calculate F1 directly from F2 and FL given by the fit
      if (imode.eq.5) then
cccc        FL=0.d0
              F1 = ((1+4.d0*x_Bj**2*am_p**2/Q2)*tmp - FL)/(2.d0*x_Bj)
      endif

c      term1 = (1.-(Q2-akt2)/Q2)
      term1 = (1.-x)*(1.d0-Q2min/Q2)

      f_D = F2/(aMX2 + Q2 - am_p**2) * term1
      f_C= 2.d0*F1/Q2
cc    f_C=0.d0
      f_ine_Budnev = alpha_em/pi*(1.-x)*akt2/Q2*
     &                (f_D+x**2/2.d0*f_C)/akt2

      return
      end
c
c     ====================================================================
c

c
c     ====================================================================
c
      function f_ine_fit(imode,x,akt2,amx)

      double precision f_ine_fit,x,x_Bj,akt2,amx,q2min,q2,tmp

      alpha_em = 1./137.035
      pi = 4.*atan(1.)
      am_p = 0.93827203
      am_pi = 0.1349766       ! mass of pi^0
c
      aMX2 = aMX*aMX
c
c     F2 structure function
c
      Q2min = 1./(1.-x)*(x*(aMX2-am_p**2)+x**2*am_p**2)
      Q2 = akt2/(1.d0-x) + Q2min
      x_Bj = Q2 / (Q2 + aMX2 - am_p**2)

      if (imode.eq.2) then
      call F2_Fiore(x_Bj,Q2,tmp)
      elseif(imode.eq.3) then
      call F2_SY(x_Bj,Q2,tmp)
      elseif(imode.eq.4) then
      call F2_ALLM(x_Bj,Q2,tmp)
      endif

c      term1 = (1.-(Q2-akt2)/Q2)
      term1 = (1.- x/2.d0 * (aMX2 - am_p**2 + Q2)/Q2)**2
      term2 = (akt2/(akt2 + x*(aMX2 - am_p**2) + x**2*am_p**2))**2

      f_aux = tmp/(aMX2 + Q2 - am_p**2) * term1 * term2
      f_ine_fit = alpha_em/pi*(1.-x)*f_aux/akt2

      return
      end
c
c     =================================================================
c
      function f_ine_old(x,akt2,amx)

      double precision f_ine_old,x,akt2,amx

c
c     a simple parametrization of unintegrated photon distribution
c
      const  = 0.004d0
      akt0 = 0.5d0        ! to be tested
      f_ine_old = const*(1.d0-x)/akt2 *(1.d0-dexp(-akt2/akt0))
     >   *1.d0/amx**4

      return
      end
c
c     =================================================================
c
      function f_ine_MRST(x,akt2,amx)

      double precision f_ine_MRST,x,x_Bj,akt2,amx,q2min,q2,tmp,q
      double precision upv,dnv,usea,dsea,str,chm,bot,glu,phot

      integer mode

      alpha_em = 1./137.035
      pi = 4.*atan(1.)
      am_p = 0.93827203
      am_pi = 0.1349766       ! mass of pi^0
c
      aMX2 = aMX*aMX
c
c     F2 structure function
c
      Q2min = 1./(1.-x)*(x*(aMX2-am_p**2)+x**2*am_p**2)
      Q2 = akt2/(1.d0-x) + Q2min
c      Q2 = x*aMX2/(1.d0-x)
      Q2=  Q2min
        x_Bj = Q2 / (Q2 + aMX2 - am_p**2)
      q = dsqrt(Q2)
c        if (q2.gt.1.d7) then
c                write(*,*) 'Q2=',Q2,'x_phot=',x
c        endif
      mode = 1
        
        if (q.lt.2.d0) then
                tmp=0.d0
        else
      call  mrstqed(x_Bj,q,mode,upv,dnv,usea,dsea,str,chm,bot,
     >  glu,phot)

      tmp = 4./9.*(upv + 2.*usea)
     2       + 1./9.*(dnv + 2.*dsea)
     3       + 1./9.*2.*str
        endif

      term1 = (1.-(Q2-akt2)/Q2)
      term2 = (akt2/(akt2 + x*(aMX2 - am_p**2) + x**2*am_p**2))**2
     
      f_aux = tmp/(aMX2 + Q2 - am_p**2) * term1 * term2
      f_ine_MRST = alpha_em/pi*(1.-x)*f_aux/akt2

      return
      end
c
c     =================================================================
c
c       ------------------------------------------------------

        subroutine F2_nasza(x,Q2,F2)
        implicit real*8 (a-h,o-z)
c       --------------------------------------------
        W2 = Q2*(1.d0-x)/x + .939**2
        W20 = 2.d0
        Q20 = 2.d0
        W2_large = 15.d0

        If (W2.gt.W2_large) then
                call F2_BDH(x,Q2,F2_block)
                F2 = F2_block
        endif        

        If ((W2.lt.W20).and.(Q2.lt.Q20)) then
                  call F2_Fiore(x,q2,F2_fi)
        endif                

        if((W2.lt.W20).and.(Q2.gt.Q20)) then
                  call F2_CTEQ(x,q2,F2_C)
                  F2 = F2_C
                endif
 
       if ((W2.gt.W20).and.(W2.lt.W2_large)) then
                call F2_SY(x,q2,F2_suri)
                F2 = F2_suri
        endif

       return 
       end

c       ------------------------------------------------------

        subroutine F2_val_CTEQ(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        dimension e(5)
        data (e(i),i=1,5) / 2.d0, -1.d0, -1.d0,
     >   2.d0, -1.d0 /

        nf = 4
        val = 0.d0
        Q = dsqrt(Q2)
        if(Q.lt.1.3d0) then
                F2 = 0.d0
                return
        endif
        do i = 1,nf
           val = val + e(i)**2/9.d0 *
     >   (Ctq6Pdf(i, X, Q) - Ctq6Pdf(-i, X, Q))
        enddo
        F2 = x*val

        return
        end

c       --------------------------------------------

        subroutine F2_CTEQ(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        dimension e(5)
        data (e(i),i=1,5) / 2.d0, -1.d0, -1.d0,
     >   2.d0, -1.d0 /

        nf = 4
        sf = 0.d0
        Q = dsqrt(Q2)
        if(Q.lt.1.3d0) then
                F2 = 0.d0
                return
        endif
        do i = 1,nf
           sf = sf + e(i)**2/9.d0 *
     >   (Ctq6Pdf(i, X, Q) + Ctq6Pdf(-i, X, Q))
        enddo
        F2 = x*sf

        return
        end

c       -------------------------------------------------------------

        subroutine F2_BDH(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        real*8 q2
c       From M.M. Block, L. Durand, P. Ha (2014)
        data a0,a1,a2 / 8.205d-4, -5.148d-2, -4.725d-3 /
        data b0,b1,b2 / 2.214d-3, 1.244d-2, 5.958d-4 /
        data c0,c1/ 0.255d0,1.475d-1/
        data xn/ 11.49d0/
        data xlambda,amu2,xm2/ 2.430d0, 2.82d0, 0.753d0/
        
        if(q2.eq.0.d0) then
        f2 = 0.d0
        return
        endif  

        tau = Q2/(Q2+amu2)
        xl = dlog(1.d0+Q2/amu2)
        xlx = dlog(tau/x)

        AA = a0 + a1*xl + a2*xl**2
        BB = b0 + b1*xl + b2*xl**2
        CC = c0  +c1*xl
        DD = Q2*(Q2+xlambda*xm2)/(Q2+xm2)**2

        F2 = DD*(1.d0-x)**xn *(CC + AA*xlx + BB*xlx**2)

        return
        end

c       -------------------------------------------------
c       -------------------------------------------------
c       --------------------------------------------------------- 
c       ========================================================
c              LUX
c       ========================================================

        subroutine F2_fit_luxlike(xbj,q2,F2,FL)
c       input: x,q2
c       output: F2,FL
        implicit real*8 (a-h,o-z)
c      double precision fmap_mstw_F2,fmap_mstw_FL
cc      common/transfer1/fmap_mstw_F2(120,120),
cc     2                fmap_mstw_FL(120,120)

cc        common/luxlike_params/amp,am_pi,alpha_em,
cc     &     q2_cut,w2_lo,w2_hi,
cc     &     ires_model,icont_model
c
c       -----------------------------
        ires_model = 1  ! Christy-Bosted
        icont_model= 1  ! GD11p
        iread = 1

c        write(*,*) "x=",xbj,"q2=",q2

        amp = 0.9382727
        am_pi = 0.135

        q2_cut = 9.d0
        W2_hi = 4.d0
        W2_lo = 3.d0
c       -----------------------------
        w_thr = amp+am_pi
        w2 = amp**2 + q2*(1.d0-xbj)/xbj
        w = dsqrt(w2)
c       -----------------------------
        omega = (w2-w2_lo)/(w2_hi-w2_lo)
        rho = 2.d0*omega**2 - omega**4 

c       write(*,*) 'x=',xbj,'q2=',q2

c        if (q2.lt.q2_cut) then
c        call F2_cont(xbj,q2,F2c,FLc)
c        F2 = F2c
c        FL = FLc
c        elseif (q2.gt.q2_cut) then
c        call F2_inter_MSTW(xbj,q2,iread,F2p,FLp)        
c        F2 = F2p
c        FL = FLp
c        endif
c        return

        if(q2.ge.q2_cut) then
           if(w2.gt.w2_hi) then ! MSTW grid
c       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                   
c        write(*,*) "going into F2_Pert!!"
         call F2_inter_MSTW(xbj,q2,iread,F2p,FLp)        
ccc             call F2_cont(xbj,q2,F2p,FLp)
ccc         call F2_inter_new(xbj,q2,F2p,FLp,fmap_mstw_F2,fmap_mstw_FL)        
              F2 = F2p
              FL = FLp
c                f2 =0.d0
c                fl =0.d0
           elseif(w2.lt.w2_hi) then
c        write(*,*) "going into F2_continuum"                   
             call F2_cont(xbj,q2,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        else ! Q2 < Q2cut
           if(w2.le.w2_lo) then
              if(ires_model.eq.1) then     ! Christy-Bosted
c        write(*,*) "going into Christy-Bosted"                   
                 call F2_res(xbj,q2,F2r,FLr)
              elseif(ires_model.eq.2) then ! Fiore-Brasse 
                 F2r = 0.d0
                 FLr = 0.d0
              endif
              F2 = F2r
              FL = FLr
c                F2 = 0.d0
c                FL= 0.d0
         elseif(w2.gt.w2_lo.and.w2.lt.w2_hi) then
c        write(*,*) "going into mixed F2"                   
              if(ires_model.eq.1) then     ! Christy-Bosted
                 call F2_res(xbj,q2,F2r,FLr)
              elseif(ires_model.eq.2) then ! Fiore-Brasse (not implemented)
                 F2r = 0.d0
                 FLr = 0.d0
              endif
               call F2_cont(xbj,q2,F2c,FLc)
c             F2 = F2c
c             FL = Flc
              F2 = (1.d0-rho)*F2r + rho*F2c
              FL = (1.d0-rho)*FLr + rho*FLc
           elseif(w2.ge.w2_hi) then
              call F2_cont(xbj,q2,F2c,FLc)
              F2 = F2c
              FL = FLc
           endif
        endif

c        write(*,*) "F2=",F2

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
        factor_mod = 1.d0
c       ------------------------------

c        if(q2.lt.q2_0) then
c                tau = 4.d0*xbj**2*amp**2/q2
c                prefac = 1.d0/(4.d0*pi**2*alpha_em)
c     >                          *q2*(1.d0-xbj)/(1+tau)
c                call christy507(w2,q2,f1,R,sigt,sigl)
c                F2 = prefac*(sigt+sigl)/0.3894e3
c                FL = F2*(1+tau)*R/(1.d0+R)
c        else
c                w2 = amp**2 + q2_mod*(1.d0-xbj)/xbj
c                tau = 4.d0*xbj**2*amp**2/q2_mod
c                prefac = 1.d0/(4.d0*pi**2*alpha_em)
c     >                          *q2_mod*(1.d0-xbj)/(1+tau)
c                call christy507(w2,q2_mod,f1,R,sigt,sigl)
c                F2 = prefac*(sigt+sigl)/0.3894e3*factor_mod
c                FL = F2*(1+tau)*R/(1.d0+R)*factor_mod
c        endif
                tau = 4.d0*xbj**2*amp**2/q2
                prefac = 1.d0/(4.d0*pi**2*alpha_em)
     >                          *q2*(1.d0-xbj)/(1+tau)
                call christy507(w2,q2,f1,R,sigt,sigl)
                F2 = prefac*(sigt+sigl)/0.3894e3
                FL = F2*(1+tau)*R/(1.d0+R)
c
        return
        end
c       ---------------------------------------------------------
        subroutine F2_CONT(xbj,q2,F2c,FLc)
        implicit real*8 (a-h,o-z)
        amp = 0.938d0

        call F2_GD11P(xbj,Q2,F2)
ccc        call F2_ALLM(xbj,Q2,F2)
        F2c = F2
        R = R_1998(xbj,Q2)
        tau = 4.d0*xbj**2*amp**2/q2
        FLc = F2c*(1+tau)*R/(1.d0+R)

        return
        end        
c       ---------------------------------------------------------       
     
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
c        if (x.lt.1.d-6.or.q2.gt.1.d5) then
c                F2 = 0.d0
c                return
c        endif

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
        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**(1.+bpom)
        F2_Reg = factor*creg*xr**areg *(1.d0-x)**(1.+breg)

cc        F2_Pom = factor*cpom*xp**apom *(1.d0-x)**bpom
cc        F2_Reg = factor*creg*xr**areg *(1.d0-x)**breg

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
        data ap1,ap2,ap3/-0.11895d0,-0.4783d0,1.353d0/
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

c====================================================================
c          subroutine F2_inter_MSTW(rx,rmu2,iread,F2,FL)
c====================================================================

      subroutine F2_inter_MSTW(xbj,q2,iread,F2,FL)

c     =================================================================
c     interpolation routine of the partonic(!) F2 or FL from grids  
c     =================================================================

cc        use f2_grid
      implicit real*8 (a-h,o-z)

c      real*8, allocatable, dimension(:,:) :: fmap_mstw_f2
c      real*8, allocatable, dimension(:,:) :: fmap_mstw_fl

c        allocate(fmap_mstw_f2(1:120,1:120))
c        allocate(fmap_mstw_fl(1:120,1:120))


c        dimension fmap_mstw_F2(120,120)
c        dimension fmap_mstw_FL(120,120)
c      integer iread
c      double precision xbj,q2,rx,rkt2,rmu2,F2,FL
c      double precision fmap_mstw_F2,fmap_mstw_FL
      common/transfer1/fmap_mstw_F2(120,120),
     2                fmap_mstw_FL(120,120)

c      save /transfer1/
c     -----------------------------------------------------------------
c     rx = log(x)
c     -----------------------------------------------------------------

c        write(*,*)"iread=",iread

        rx = dlog10(xbj)
        rmu2 = dlog10(q2)

      rx_min = -6.d0 
      rx_max =  0.d0
c      nrx = 120    
      nrx = 120    
      drx = (rx_max-rx_min)/nrx

c     -----------------------------------------------------------------
c     rmu2 = log(mu2)
c     -----------------------------------------------------------------

      rmu2_min =  0.d0 
      rmu2_max =  6.d0
c      nrmu2 = 120
      nrmu2 = 120
      drmu2 = (rmu2_max-rmu2_min)/nrmu2  

c     =================================================================
c     first time, only reading
c     =================================================================

      if(iread.eq.0) then

      do i = 1,120
      do j = 1,120

cc        write(*,*) fmap_mstw_F2(i,j)

      fmap_mstw_F2(i,j) = 0.d0
      fmap_mstw_FL(i,j) = 0.d0

      enddo
      enddo

c     =================================================================
c     reading files with grids
c     =================================================================
c     
      open(unit=11,file='F2_MSTW_nnlo.dat',
     &     status='unknown')

      open(unit=12,file='FL_MSTW_nnlo.dat',
     &     status='unknown') 

c     =================================================================

      do 100 irx=1,nrx

      rx = rx_min + irx*drx - drx/2.

      do 100 irmu2=1,nrmu2

      rmu2 = rmu2_min + irmu2*drmu2 - drmu2/2.

      read(11,*) alog10_x, alog10_mu2, F2
      read(12,*) alog10_x, alog10_mu2, FL

      fmap_mstw_F2(irx,irmu2) = dble(F2)
      fmap_mstw_FL(irx,irmu2) = dble(FL)

  100 continue

      close(unit=11)
      close(unit=12)

      elseif(iread.eq.1) then

c    ==================================================================
c     outside of the grid
c    ==================================================================

      if(rmu2.lt.(rmu2_min+drmu2/2.).or.rmu2.gt.(rmu2_max-drmu2/2.))
     & then
      F2 = 0.d0
      FL = 0.d0
c     goto 200       
        return
      endif

      if(rx.lt.(rx_min+drx/2.).or.rx.gt.(rx_max-drx/2.))
     & then
      F2 = 0.d0
      FL = 0.d0
c      goto 200       
        return
      endif

     
c     =================================================================
c     interpolation variables
c     =================================================================

      delrx = (rx-rx_min)
      delrmu2 = (rmu2-rmu2_min)
      
      srx = (delrx-drx/2.) / drx
      irx = (delrx-drx/2.) / drx
      irx_lo = irx+1
      irx_up = irx_lo+1

      srmu2 = (delrmu2-drmu2/2.) / drmu2
      irmu2 = (delrmu2-drmu2/2.) / drmu2
      irmu2_lo = irmu2+1
      irmu2_up = irmu2_lo+1

c     =================================================================
c     for testing the code only
c     =================================================================

      rx_lo = rx_min+irx_lo*drx-drx/2.
      rx_up = rx_min+irx_up*drx-drx/2.
      rmu2_lo = rmu2_min+irmu2_lo*drmu2-drmu2/2.
      rmu2_up = rmu2_min+irmu2_up*drmu2-drmu2/2.

c     =================================================================
c     neighbouring grid points
c     =================================================================

      F2_11 = fmap_mstw_F2(irx_lo,irmu2_lo)
      F2_12 = fmap_mstw_F2(irx_lo,irmu2_up)
      F2_21 = fmap_mstw_F2(irx_up,irmu2_lo)
      F2_22 = fmap_mstw_F2(irx_up,irmu2_up)

      FL_11 = fmap_mstw_FL(irx_lo,irmu2_lo)
      FL_12 = fmap_mstw_FL(irx_lo,irmu2_up)
      FL_21 = fmap_mstw_FL(irx_up,irmu2_lo)
      FL_22 = fmap_mstw_FL(irx_up,irmu2_up)

c        write(*,*) "map F2:",F2_11,F2_12,F2_21,F2_22
c        write(*,*) "map FL:",FL_11,FL_12,FL_21,FL_22


c     =================================================================
c     weights for neighbouring points
c     =================================================================

      ddrx = (srx - float(irx))*drx
      ddrmu2  = (srmu2 - float(irmu2))*drmu2

      w11  =   (drx-ddrx)      /drx 
     *       * (drmu2-ddrmu2)/drmu2
      w12  =   (drx-ddrx)      /drx 
     *       *   ddrmu2        /drmu2
      w21  =     ddrx            /drx 
     *       * (drmu2-ddrmu2)/drmu2
      w22  =     ddrx            /drx 
     *       *   ddrmu2         /drmu2

c     =========================================
c     interpolation
c     =========================================

      f_int_F2 =  w11*F2_11 + w12*F2_12 
     *          + w21*F2_21 + w22*F2_22

      f_int_FL =  w11*FL_11 + w12*FL_12 
     *          + w21*FL_21 + w22*FL_22

      F2 = f_int_F2
      FL = f_int_FL

  200 continue   

c        deallocate(fmap_mstw_f2(1:120,1:120))
c        deallocate(fmap_mstw_fl(1:120,1:120))
      endif

      return

      end

c     ====================================================

      subroutine F2_inter_new(xbj,q2,F2,FL,fmap_mstw_F2,fmap_mstw_FL)

c     =================================================================
c     interpolation routine of the partonic(!) F2 or FL from grids  
c     =================================================================

      implicit real*8 (a-h,o-z)
      dimension fmap_mstw_F2(120,120),fmap_mstw_FL(120,120)

c      save fmap_mstw_F2, fmap_mstw_FL

c      integer iread
c      double precision xbj,q2,rx,rkt2,rmu2,F2,FL
c      double precision fmap_mstw_F2,fmap_mstw_FL
cc      common/transfer/fmap_mstw_F2(120,120)
cc     2                fmap_mstw_FL(120,120)

c      save fmap_mstw_F2,fmap_mstw_FL

c     -----------------------------------------------------------------
c     rx = log(x)
c     -----------------------------------------------------------------

c        write(*,*)"iread=",iread

        rx = dlog10(xbj)
        rmu2 = dlog10(q2)

      rx_min = -6.d0 
      rx_max =  0.d0
      nrx = 120    
      drx = (rx_max-rx_min)/nrx

c     -----------------------------------------------------------------
c     rmu2 = log(mu2)
c     -----------------------------------------------------------------

      rmu2_min =  0.d0 
      rmu2_max =  6.d0
      nrmu2 = 120
      drmu2 = (rmu2_max-rmu2_min)/nrmu2  
c    ==================================================================
c     outside of the grid
c    ==================================================================

      if(rmu2.lt.(rmu2_min+drmu2/2.).or.rmu2.gt.(rmu2_max-drmu2/2.))
     & then
      F2 = 0.d0
      FL = 0.d0
      goto 200       
      endif

      if(rx.lt.(rx_min+drx/2.).or.rx.gt.(rx_max-drx/2.))
     & then
      F2 = 0.d0
      FL = 0.d0
      goto 200       
      endif

c     =================================================================
c     interpolation variables
c     =================================================================

      delrx = (rx-rx_min)
      delrmu2 = (rmu2-rmu2_min)
      
      srx = (delrx-drx/2.) / drx
      irx = (delrx-drx/2.) / drx
      irx_lo = irx+1
      irx_up = irx_lo+1

      srmu2 = (delrmu2-drmu2/2.) / drmu2
      irmu2 = (delrmu2-drmu2/2.) / drmu2
      irmu2_lo = irmu2+1
      irmu2_up = irmu2_lo+1

c     =================================================================
c     for testing the code only
c     =================================================================

      rx_lo = rx_min+irx_lo*drx-drx/2.
      rx_up = rx_min+irx_up*drx-drx/2.
      rmu2_lo = rmu2_min+irmu2_lo*drmu2-drmu2/2.
      rmu2_up = rmu2_min+irmu2_up*drmu2-drmu2/2.

c     =================================================================
c     neighbouring grid points
c     =================================================================

      F2_11 = fmap_mstw_F2(irx_lo,irmu2_lo)
      F2_12 = fmap_mstw_F2(irx_lo,irmu2_up)
      F2_21 = fmap_mstw_F2(irx_up,irmu2_lo)
      F2_22 = fmap_mstw_F2(irx_up,irmu2_up)

      FL_11 = fmap_mstw_FL(irx_lo,irmu2_lo)
      FL_12 = fmap_mstw_FL(irx_lo,irmu2_up)
      FL_21 = fmap_mstw_FL(irx_up,irmu2_lo)
      FL_22 = fmap_mstw_FL(irx_up,irmu2_up)

cc        write(*,*) "map F2:",F2_11,F2_12,F2_21,F2_22
cc        write(*,*) "map FL:",FL_11,FL_12,FL_21,FL_22


c     =================================================================
c     weights for neighbouring points
c     =================================================================

      ddrx = (srx - float(irx))*drx
      ddrmu2  = (srmu2 - float(irmu2))*drmu2

      w11  =   (drx-ddrx)      /drx 
     *       * (drmu2-ddrmu2)/drmu2
      w12  =   (drx-ddrx)      /drx 
     *       *   ddrmu2        /drmu2
      w21  =     ddrx            /drx 
     *       * (drmu2-ddrmu2)/drmu2
      w22  =     ddrx            /drx 
     *       *   ddrmu2         /drmu2

c     =========================================
c     interpolation
c     =========================================

      f_int_F2 =  w11*F2_11 + w12*F2_12 
     *          + w21*F2_21 + w22*F2_22

      f_int_FL =  w11*FL_11 + w12*FL_12 
     *          + w21*FL_21 + w22*FL_22

      F2 = f_int_F2
      FL = f_int_FL

  200 continue   


      return

      end



c       -------------------------------------------------
c       -------------------------------------------------
c       -------------------------------------------------

        subroutine F2_Fiore(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        real*8 q2

        dimension iswitch(4)
        dimension fit(2,4,5),s0(2),anorm(2),spin(4)

        amp= 0.939d0
        alpha_em = 1.d0/137.d0
        pi = dacos(-1.d0)
        akin = 1.d0 + 4.d0*amp**2*x**2/q2

        prefactor = q2*(1.d0-x)
     >  /(4.d0*pi*alpha_em*akin)

c        write(*,*) 'prefactor=',prefactor

        s = Q2*(1.d0-x)/x + amp**2 

        s0(1) = 1.14d0
        s0(2) = 1.2871d0

        anorm(1) = 0.021d0
        anorm(2) = 0.0207d0

        spin(1) = 3.d0/2.d0
        spin(2) = 5.d0/2.d0
        spin(3) = 3.d0/2.d0
        spin(4) = 1.d0

        iswitch(1) = 1
        iswitch(2) = 1
        iswitch(3) = 1
        iswitch(4) = 1

c       --------------------------

        fit(1,1,1) = -0.8377d0
        fit(1,1,2) =  0.95d0
        fit(1,1,3) =  0.1473d0
        fit(1,1,4) =  1.d0
        fit(1,1,5) =  2.4617d0

        fit(1,2,1) = -0.37d0
        fit(1,2,2) =  0.95d0
        fit(1,2,3) =  0.1471d0
        fit(1,2,4) =  0.5399d0
        fit(1,2,5) =  2.4617d0

        fit(1,3,1) =  0.0038d0
        fit(1,3,2) =  0.85d0
        fit(1,3,3) =  0.1969d0
        fit(1,3,4) =  4.2225d0
        fit(1,3,5) =  1.5722d0

        fit(1,4,1) =  0.5645d0
        fit(1,4,2) =  0.1126d0
        fit(1,4,3) =  1.3086d0
        fit(1,4,4) =  19.2694d0
        fit(1,4,5) =  4.5259d0

c       --------------------------

        fit(2,1,1) = -0.8070d0
        fit(2,1,2) =  0.9632d0
        fit(2,1,3) =  0.1387d0
        fit(2,1,4) =  1.d0
        fit(2,1,5) =  2.6066d0

        fit(2,2,1) = -0.3640d0
        fit(2,2,2) =  0.9531d0
        fit(2,2,3) =  0.1239d0
        fit(2,2,4) =  0.6086d0
        fit(2,2,5) =  2.6066d0

        fit(2,3,1) =  -0.0065d0
        fit(2,3,2) =  0.8355d0
        fit(2,3,3) =  0.2320d0
        fit(2,3,4) =  4.7279d0
        fit(2,3,5) =  1.4828d0

        fit(2,4,1) =  0.5484d0
        fit(2,4,2) =  0.1373d0
        fit(2,4,3) =  1.3139d0
        fit(2,4,4) =  14.7267d0
        fit(2,4,5) =  4.6041d0

c       --------------------------------     

c        ifit = 2
        ifit = 1
        ampli_res = 0.d0
        do ires = 1,3

        alpha_0 = fit(ifit,ires,1) 
        alpha_1 = fit(ifit,ires,2) 
        alpha_2 = fit(ifit,ires,3) 

        s_thr = s0(ifit)

        if (s.gt.s_thr) then

        alpha_Re = alpha_0 + alpha_2*dsqrt(s_thr) + alpha_1*s
        alpha_Im = alpha_2*dsqrt(s-s_thr)

        else

        alpha_Re = alpha_0 + alpha_1*s + 
     >          alpha_2*(dsqrt(s_thr) - dsqrt(s_thr - s))
        
        alpha_Im = 0.d0

        endif


        Q02 = fit(ifit,ires,5) 
        formfactor = 1.d0/(1.d0 + Q2/Q02)**2

        a = fit(ifit,ires,4)
        sp = spin(ires)
        denom = (sp-alpha_Re)**2 + alpha_Im**2

        ampli_imag = a*formfactor**2*alpha_Im/denom


        j = iswitch(ires)

        ampli_res = ampli_res + j* ampli_imag
        enddo

        s_E = fit(ifit,4,3)
        aBG_0 = fit(ifit,4,1)
        aBG_2 = fit(ifit,4,2)

        if (s.gt.s_E) then

        alpha_Re = aBG_0 + aBG_2*dsqrt(s_E)
        alpha_Im = aBG_2*dsqrt(s-s_E)

        else

        alpha_Re = aBG_0 +
     >          aBG_2*(dsqrt(s_E) - dsqrt(s_E - s))
        
        alpha_Im = 0.d0

        endif

        a = fit(ifit,4,4)
        Q02 =  fit(ifit,4,5)

        formfactor = 1.d0/(1.d0 + Q2/Q02)**2

        sp = 1.5*spin(4)
        denom = (sp-alpha_Re)**2 + alpha_Im**2

        j = iswitch(4)
        ampli_bg = a*formfactor**2*alpha_Im/denom
        an = anorm(ifit)
        ampli_tot = an*(ampli_res+j*ampli_bg)


        F2 = prefactor*ampli_tot

        return
        end

c       --------------------------------------------------------

        subroutine F2_SU(x_val,Q2_val,F2)
        real*8 x_val,q2_val,F2
        data q02/ 0.8 /

        x = real(x_val)
        q2 = real(q2_val)

        amu2 = q2 + q02

        call grv95lo(x,amu2,xuv,xdv,xus,xds,xss,xg)


        F2_aux = 4./9.*(xuv + 2.*xus)
     >        + 1./9.*(xdv + 2.*xds)
     >        + 1./9.*2.*xss

        F2 = dble((Q2 / amu2) * F2_aux)

        return
        end

c       ---------------------------------------------------------       

        subroutine F2_SY(x,Q2,F2)
        implicit real*8 (a-h,o-z)
        amp= 0.939

        rS = Q2*(1.d0-x)/x + amp**2
        rQs = Q2
        rNu = (Q2 + rS - amp**2)/2.d0/amp

	call SuriYennieSF(rS, rQs, rW1, rW2) 	

        F2 = rNu/amp * rW2

        return
        end


c       -------------------------------------------------
c-----------------------------------------------------------------------
*  Suri & Yennie structure functions.
*
*  This is a piece of ZLPAIR code rewritten by A.Poluektov, BINP.
*
*  Alternative set of parameters mentioned in the paper of Suri & Yennie
*  is added as a comment.
*
*  Thanks to S.P.Baranov, O.Dinger, H.Shooshtari and J.A.M.Vermaseren.
*-----------------------------------------------------------------------
	subroutine SuriYennieSF(rS, rQs, rW1, rW2) 	
	implicit double precision (r)

	parameter (rMp = 0.9383,rMps = rMp**2)
	
	common /SuriYennieSFcb/ rC1,rC2,rD1,rMrhos,rCp,rBp

	data rC1/0.86926/,rC2/2.23422/,rD1/0.12549/  ! LPAIR default
	data rMrhos/0.585/,rCp/0.96/,rBp/0.63/

!	data rC1/0.6303/, rC2/2.3049/, rD1/0.04681/  ! alternative
!	data rMrhos/1.05/,rCp/1.23/,rBp/0.61/

	save /SuriYennieSFcb/ 

	rSMs = rS - rMps
	
	rXpr = rQs/(rQs + rS)
	rEn = rSMs + rQs
	rTau = rQs / (4.*rMps)
	rMQ = (rMrhos + rQs)
	
	rC = (rC1*rSMs*(rMrhos/rMQ)**2 + 
	1	(rC2*rMps*(1.-rXpr)**4)/
	2	(1.+rXpr*(rXpr*rCp-2.*rBp)))/rQs
!	2	(1.+rXpr*(rXpr*rCp+2.*rBp)))/rQs

	rD=(rTau*rC + 
	1	rD1*rSMs*rQs*rMrhos/rMps*(rSMs/rMQ/rEn)**2)/ 
	2	( 1.+(rEn**2)/(4.*rMps*rQs) )

	rW2 = 2.*rMp*rD
	rW1 = rC*rQs/2./rMp

	end

c---------------------------------------------------------------------------

      subroutine ampli_WW(shat,that,uhat,lam1,lam2,lam3,lam4,ampli)
      
      double precision  shat,that,uhat,ampli
      double precision  beta,gam,am_l,alpha_em
      double precision  e2,pi,costheta,sintheta,A
      double precision  term1,term2,term3,term4
      integer ilam1,ilam2,ilam3,ilam4
c
c     fundamental constants
c
      pi = 4.*atan(1.)
      alpha_em = 1./137.035
      am_l = 80.385
      e2 = 4.d0*pi*alpha_em
      t1 = q12
      t2 = q22
      
c      costheta = (that-uhat)/(shat-2.d0*am_l**2-t1-t2)
c      costheta = (1.d0+ 2.d0*(that-am_l**2)/shat)
c     &  /sqrt(1.d0+1.d-10-4.d0*am_l**2/shat)
        costheta = (that-uhat)/shat
     >       /sqrt(1.d0+1.d-10-4.d0*am_l**2/shat)

cc        uu = (1.d0 - costheta)*(1.d0 + costheta)
       sintheta = dsqrt(1.d0 -costheta**2)
cc        sintheta = dsqrt(1.d0-uu)

      beta = dsqrt(1.d0-4.d0*am_l**2/shat)  
      gam = 1.d0/dsqrt(1.d0- beta**2)
c        gam = sqrt(shat)/2./am_l
      A = (1.d0-beta**2*costheta**2)
c      B = shat-am_H**2


      term1 = 1.d0/(gam**2*A)*((gam**2+1.d0)*(1.d0-lam1*lam2)
     &        *sintheta**2
     &        -(1.d0+lam1*lam2))
      term2 = (-1.d0)*sqrt(2.)/(gam*A)*(lam1-lam2)
     &        *(1.d0+lam1*lam3*costheta)*sintheta
      term3 = (-1.d0)/(2.d0*A)*(2.d0*beta*(lam1+lam2)*(lam3+lam4) 
     &		- (1.d0/gam**2)
     &        *(1.d0+lam3*lam4)*(2.d0*lam1*lam2+(1.d0-lam1*lam2)
     &         *costheta**2)
     &        +(1.d0+lam1*lam2*lam3*lam4)*(3.d0+lam1*lam2)
     &        + 2.d0*(lam1-lam2)*(lam3-lam4)*costheta
     &        + (1.d0-lam1*lam2)*(1.d0-lam3*lam4)*costheta**2) 
      term4 = (-1.d0)*sqrt(2.)/(gam*A)*(lam2-lam1)
     &        *(1.d0+lam2*lam4*costheta)*sintheta



      if(lam3.eq.0.and.lam4.eq.0) then
		ampli = e2*term1
      else if(lam4.eq.0) then
		ampli = e2*term2
      else if(lam3.eq.0) then
		ampli = e2*term4
      else if(lam3.ne.0.and.lam4.ne.0) then
		ampli = e2*term3
      endif


      return
      end
*










