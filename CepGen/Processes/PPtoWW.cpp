#include "PPtoWW.h"

using namespace CepGen::Process;

PPtoWW::PPtoWW() :
  GenericKTProcess( "pptoww", "gamma,gamma->W+,W-", 0 /*FIXME*/, { Particle::Photon, Particle::Photon }, { Particle::W, Particle::W } )
{}

void
PPtoWW::setExtraContent()
{
  Particles& wbosons = event_->getByRole( Particle::CentralSystem );

  //--- add the decay products for the Ws

  Particle lep1( Particle::CentralSystem, Particle::Muon );
  lep1.addMother( wbosons[0] );
  event_->addParticle( lep1 );

  Particle nu1( Particle::CentralSystem, Particle::MuonNeutrino );
  nu1.addMother( wbosons[0] );
  event_->addParticle( nu1 );

  Particle lep2( Particle::CentralSystem, Particle::Electron );
  lep2.addMother( wbosons[1] );
  event_->addParticle( lep2 );

  Particle nu2( Particle::CentralSystem, Particle::ElectronNeutrino );
  nu2.addMother( wbosons[1] );
  event_->addParticle( nu2 );

  //--- reinitialise the event content to add these extra decay products by default

  event_->init();
}

void
PPtoWW::prepareKTKinematics()
{}

double
PPtoWW::computeJacobian()
{
  double jac = GenericKTProcess::minimalJacobian();
  jac *= cuts_.central_cuts[Cuts::rapidity_single].range(); // d(y1)
  jac *= cuts_.central_cuts[Cuts::rapidity_single].range(); // d(y2)
  jac *= cuts_.central_cuts[Cuts::pt_diff].range(); // d(Dpt)
  jac *= 2.*M_PI; // d(phiDpt)

  return jac;
}

double
PPtoWW::computeKTFactorisedMatrixElement()
{
  //=================================================================
  //return aintegral*_q1t*_q2t*pt_diff_;
  //=================================================================
  return 0.;
}

void
PPtoWW::fillCentralParticlesKinematics()
{
  // randomise the charge of the outgoing W boson
  int sign = ( drand()>.5 ) ? +1 : -1;

  Particles& fs = event_->getByRole( Particle::CentralSystem );

  //=================================================================
  //     outgoing Ws
  //=================================================================

  fs[0].setPdgId( fs[0].pdgId(), sign );
  fs[0].setStatus( Particle::Undecayed );
  fs[0].setMomentum( p_w1_ );

  fs[1].setPdgId( fs[1].pdgId(), -sign);
  fs[1].setStatus( Particle::Undecayed );
  fs[1].setMomentum( p_w2_ );

  //=================================================================
  //     final state particles
  //=================================================================

  fs[2].setPdgId( fs[2].pdgId(), sign );
  fs[3].setPdgId( fs[3].pdgId(), -sign );
  fs[4].setPdgId( fs[4].pdgId(), -sign );
  fs[5].setPdgId( fs[5].pdgId(), sign );

  event_->dump();
}
