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
  jac *= ( y_max_-y_min_ ); // d(y1)
  jac *= ( y_max_-y_min_ ); // d(y2)
  jac *= ( cuts_.ptdiff_max-cuts_.ptdiff_min ); // d(Dpt)
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
  //     first outgoing W
  //=================================================================
  Particle& w1 = event_->getOneByRole( Particle::CentralParticle1 );
  w1.setPdgId( w1.pdgId(), sign );
  w1.setStatus( Particle::Undecayed );
  w1.setMomentum( p_w1_ );

  //=================================================================
  //     second outgoing W
  //=================================================================
  Particle& w2 = event_->getOneByRole( Particle::CentralParticle2 );
  w2.setPdgId( w2.pdgId(), -sign);
  w2.setStatus( Particle::Undecayed );
  w2.setMomentum( p_w2_ );
}
