c     ================================================================
c     CepGen interface
c     ================================================================

c     =================================================================
c     collection of fundamental physics constants
c     =================================================================
      common/constants/am_p,units,pi,alpha_em
      double precision am_p,units,pi,alpha_em

c     =================================================================
c     information on the full process
c       inp1 = proton energy in lab frame
c       inp2 = nucleus energy **per nucleon** in LAB frame
c       Collision is along z-axis
c     =================================================================
      common/params/mode,sfmod,pdg_l,a_nuc,z_nuc,am_l,inp1,inp2
      integer mode,sfmod,pdg_l,a_nuc,z_nuc
      double precision am_l,inp1,inp2

c     =================================================================
c     kt-factorisation kinematics
c     =================================================================
      common/ktkin/q1t,q2t,phiq1t,phiq2t,y1,y2,ptdiff,phiptdiff,
     &     am_x,am_y
      double precision q1t,q2t,phiq1t,phiq2t,y1,y2,ptdiff,phiptdiff,
     &     am_x,am_y

c     =================================================================
c     phase space cuts
c     =================================================================
      common/kincuts/ipt,ieta,idely,pt_min,pt_max,eta_min,eta_max,
     &      dely_min,dely_max
      integer ipt,ieta,idely
      double precision pt_min,pt_max,eta_min,eta_max,dely_min,dely_max

c     =================================================================
c     generated event kinematics
c     =================================================================
      common/evtkin/p10,p1x,p1y,p1z,p20,p2x,p2y,p2z,
     &     px_0,px_x,px_y,px_z,py_0,py_x,py_y,py_z
      double precision p10,p1x,p1y,p1z,p20,p2x,p2y,p2z,
     &     px_0,px_x,px_y,px_z,py_0,py_x,py_y,py_z

c     =================================================================
c     helpers for the evaluation of unintegrated parton fluxes
c     =================================================================
      external CepGen_kT_flux,CepGen_kT_flux_HI
      double precision CepGen_kT_flux,CepGen_kT_flux_HI

