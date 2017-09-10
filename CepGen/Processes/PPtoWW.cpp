#include "PPtoWW.h"

using namespace CepGen::Process;

PPtoWW::PPtoWW() : GenericKTProcess("gamma,gamma->W+,W-", 0 /*FIXME*/, Particle::Photon, Particle::WPlus)
{}

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

  Particles& wbosons = event_->getByRole( Particle::CentralSystem );

  wbosons[0].setPdgId( wbosons[0].pdgId(), sign );
  wbosons[0].setStatus( Particle::Undecayed );
  wbosons[0].setMomentum( p_w1_ );

  wbosons[1].setPdgId( wbosons[1].pdgId(), -sign);
  wbosons[1].setStatus( Particle::Undecayed );
  wbosons[1].setMomentum( p_w2_ );
}
