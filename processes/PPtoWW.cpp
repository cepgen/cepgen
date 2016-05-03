#include "PPtoWW.h"

PPtoWW::PPtoWW() : GenericKTProcess("gamma,gamma->W+,W-", 0 /*FIXME*/, Particle::Photon, Particle::WPlus)
{}

void
PPtoWW::PrepareKTKinematics()
{}

double
PPtoWW::ComputeJacobian()
{
  double jac = 1.;
  jac *= (fLogQmax-fLogQmin)*fQT1; // d(q1t) . q1t
  jac *= (fLogQmax-fLogQmin)*fQT2; // d(q2t) . q2t
  jac *= 2.*Constants::Pi; // d(phi1)
  jac *= 2.*Constants::Pi; // d(phi2)
  jac *= (fYmax-fYmin); // d(y1)
  jac *= (fYmax-fYmin); // d(y2)
  switch (fCuts.kinematics) {
    case 1: default: break;
    case 2: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*fMY; break;
    case 3: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*fMX; break;
    case 4: jac *= (fCuts.mxmax-fCuts.mxmin)*2.*fMX;
            jac *= (fCuts.mxmax-fCuts.mxmin)*2.*fMY; break;
  } // d(mx/y**2)
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
  if (!w1->SetMomentum(fPw1)) { Error("Invalid outgoing W1"); }

  //=================================================================
  //     second outgoing W
  //=================================================================
  Particle* w2 = GetParticle(Particle::CentralParticle2);
  w2->SetPDGId(w2->GetPDGId(), -sign);
  w2->status = Particle::Undecayed;
  if (!w2->SetMomentum(fPw2)) { Error("Invalid outgoing W2"); }
}
