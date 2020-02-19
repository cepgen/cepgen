      !> \file KTBblocks.inc
c     ================================================================
c     F77 process-to-CepGen helper
c     > this header file defines a set of common blocks useful for
c     > the interfacing of any kt-factorised matrix element computation
c     > to a mother CepGen instance
c     ================================================================

c     =================================================================
      !> collection of fundamental physics constants
      common/constants/am_p,units,pi,alpha_em
      double precision am_p,units,pi,alpha_em

c     =================================================================
c     information on the full process
c       inp1 = proton energy in lab frame
c       inp2 = nucleus energy **per nucleon** in LAB frame
c       Collision is along z-axis
c     =================================================================
      common/genparams/icontri,iflux1,iflux2,sfmod,
     &     a_nuc1,z_nuc1,a_nuc2,z_nuc2,inp1,inp2
      integer icontri,iflux1,iflux2,sfmod,
     &     a_nuc1,z_nuc1,a_nuc2,z_nuc2
      double precision inp1,inp2

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
      common/kincuts/ipt,iene,ieta,iinvm,iptsum,idely,
     &     pt_min,pt_max,ene_min,ene_max,eta_min,eta_max,
     &     invm_min,invm_max,ptsum_min,ptsum_max,
     &     dely_min,dely_max
      logical ipt,iene,ieta,iinvm,iptsum,idely
      double precision pt_min,pt_max,ene_min,ene_max,eta_min,eta_max,
     &     invm_min,invm_max,ptsum_min,ptsum_max,
     &     dely_min,dely_max

c     =================================================================
c     generated event kinematics
c     =================================================================
      common/evtkin/nout,ipdg,idum,pc,px,py
      integer nout,idum,ipdg(10)
      double precision pc(4,10),px(4),py(4)

c     =================================================================
c     helpers for the evaluation of unintegrated parton fluxes
c     =================================================================
      external CepGen_kT_flux,CepGen_kT_flux_HI
      double precision CepGen_kT_flux,CepGen_kT_flux_HI
      external CepGen_particle_charge,CepGen_particle_mass
      external CepGen_particle_colour
      double precision CepGen_particle_charge,CepGen_particle_mass
      double precision CepGen_particle_colour

c     =================================================================
c     helpers for the retrieval of input parameters from the process
c     =================================================================
      external CepGen_param_int,CepGen_param_real
      integer CepGen_param_int
      double precision CepGen_param_real

