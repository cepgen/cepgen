#include "PPtoWW.h"

using namespace CepGen::Process;

PPtoWW::PPtoWW() :
  GenericKTProcess( "pptoww", "gamma,gamma->W+,W-", 0 /*FIXME*/, Particle::Photon, Particle::WPlus )
{}

void
PPtoWW::setExtraContent()
{
  Particles& wbosons = event_->getByRole( Particle::CentralSystem );

  //--- add the decay products for the Ws

  Particle lep1( Particle::CentralSystem, Particle::Muon ), nu1( Particle::CentralSystem, Particle::MuonNeutrino );
  lep1.addMother( wbosons[0] );
  nu1.addMother( wbosons[0] );
  event_->addParticle( lep1 );
  event_->addParticle( nu1 );

  Particle lep2( Particle::CentralSystem, Particle::Electron ), nu2( Particle::CentralSystem, Particle::ElectronNeutrino );
  lep2.addMother( wbosons[1] );
  nu2.addMother( wbosons[1] );
  event_->addParticle( lep2 );
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

  //=================================================================
  //     outgoing Ws
  //=================================================================

  /*Particles& wbosons = event_->getByRole( Particle::CentralSystem );

  wbosons[0].setPdgId( wbosons[0].pdgId(), sign );
  wbosons[0].setStatus( Particle::Undecayed );
  wbosons[0].setMomentum( p_w1_ );

  wbosons[1].setPdgId( wbosons[1].pdgId(), -sign);
  wbosons[1].setStatus( Particle::Undecayed );
  wbosons[1].setMomentum( p_w2_ );*/
}
