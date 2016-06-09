#include "PPtoWW.h"

PPtoWW::PPtoWW() : GenericKTProcess("gamma,gamma->W+,W-", 0 /*FIXME*/, Particle::Photon, Particle::WPlus)
{}

void
PPtoWW::PrepareKTKinematics()
{}

double
PPtoWW::ComputeJacobian()
{
  double jac = GenericKTProcess::MinimalJacobian();
  jac *= (fYmax-fYmin); // d(y1)
  jac *= (fYmax-fYmin); // d(y2)
  jac *= (fCuts.ptdiffmax-fCuts.ptdiffmin); // d(Dpt)
  jac *= 2.*Constants::Pi; // d(phiDpt)
  
  return jac;
}

double
PPtoWW::ComputeKTFactorisedMatrixElement()
{
  //=================================================================
  //return aintegral*_q1t*_q2t*_ptdiff;
  //=================================================================
  return 0.;
}

void
PPtoWW::FillCentralParticlesKinematics()
{
  // randomise the charge of the outgoing W boson
  int sign = (drand()>.5) ? +1 : -1;

  //=================================================================
  //     first outgoing W
  //=================================================================
  Particle* w1 = GetParticle(Particle::CentralParticle1);
  w1->SetPDGId(w1->GetPDGId(), sign);
  w1->status = Particle::Undecayed;
  if (!w1->SetMomentum(fPw1)) { InError("Invalid outgoing W1"); }

  //=================================================================
  //     second outgoing W
  //=================================================================
  Particle* w2 = GetParticle(Particle::CentralParticle2);
  w2->SetPDGId(w2->GetPDGId(), -sign);
  w2->status = Particle::Undecayed;
  if (!w2->SetMomentum(fPw2)) { InError("Invalid outgoing W2"); }
}
