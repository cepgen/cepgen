#include "GenericProcess.h"

GenericProcess::GenericProcess(std::string name_) :
  fX(0), fNumDimensions(0), fEvent(new Event), fIsPointSet(false),
  fIsInStateSet(false), fIsOutStateSet(false), fIsKinematicSet(false),
  fName(name_)
{}

GenericProcess::~GenericProcess()
{
  if (fIsPointSet) delete[] fX;
  delete fEvent;
}

void
GenericProcess::SetPoint(const unsigned int ndim_,double x_[])
{
  // Number of dimensions on which the integration will be performed
  fNumDimensions = ndim_;

  // Phase space coordinate becomes a protected attribute
  if (!fX) fX = new double[ndim_];

  std::copy(x_, x_+ndim_, fX);  
  fIsPointSet = true;
  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) { DumpPoint(DebugMessage); }
}

void
GenericProcess::PrepareKinematics()
{
  if (!IsKinematicsDefined()) return; // FIXME dump some information...
  fSqS = CMEnergy(*GetParticle(Particle::IncomingBeam1),
                  *GetParticle(Particle::IncomingBeam2));
  fS = pow(fSqS, 2);
  
  Debugging(Form("Kinematics successfully prepared! sqrt(s) = %.2f", fSqS));
}

void
GenericProcess::DumpPoint(const ExceptionType& et=Information)
{
  std::ostringstream os;
  for (unsigned int i=0; i<(unsigned int)fNumDimensions; i++) {
    os << Form("  x(%2d) = %8.6f\n\t", i, fX[i]);
  }
  if (et<DebugMessage) { Information(Form("Number of integration parameters: %d\n\t"
                                          "%s", fNumDimensions, os.str().c_str())); }
  else              { Debugging(Form("Number of integration parameters: %d\n\t"
                                     "%s", fNumDimensions, os.str().c_str())); }
}

void
GenericProcess::SetEventContent(IncomingState is, OutgoingState os)
{  
  for (IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++) { fEvent->AddParticle(Particle(ip->first, ip->second)); }

  // Prepare the central system if not already there
  IncomingState::const_iterator central_system = is.find(Particle::CentralSystem);
  if (central_system==is.end()) {
    Particle* moth = GetParticle(Particle::Parton1);
    Particle cs(Particle::CentralSystem, moth->GetPDGId());
    cs.SetMother(moth);
    fEvent->AddParticle(cs);
  }

  for (OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++) { fEvent->AddParticle(Particle(op->first, op->second)); }
  
  // Incoming particles (incl. eventual partons)
  for (IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++) {
    Particle* p = GetParticle(ip->first);
    p->status = Particle::Undefined;
    switch (ip->first) {
      case Particle::IncomingBeam1:
      case Particle::IncomingBeam2: break;
      case Particle::Parton1:       p->SetMother(GetParticle(Particle::IncomingBeam1)); break;
      case Particle::Parton2:       p->SetMother(GetParticle(Particle::IncomingBeam2)); break;
      case Particle::CentralSystem: p->SetMother(GetParticle(Particle::Parton1)); break;
      default: break;
    }
  }
  // Outgoing particles (central, and outgoing primary particles or remnants)
  for (OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++) {
    Particle* p = GetParticle(op->first);
    p->status = Particle::Undefined;
    switch (op->first) {
      case Particle::OutgoingBeam1:    p->SetMother(GetParticle(Particle::IncomingBeam1)); break;
      case Particle::OutgoingBeam2:    p->SetMother(GetParticle(Particle::IncomingBeam2)); break;
      case Particle::CentralParticle1: p->SetMother(GetParticle(Particle::CentralSystem)); break;
      case Particle::CentralParticle2: p->SetMother(GetParticle(Particle::CentralSystem)); break;
      default: break;
    }
  }
  fEvent->Init();
}

void
GenericProcess::SetIncomingKinematics(Particle::Momentum p1, Particle::Momentum p2)
{
  if (!GetParticle(Particle::IncomingBeam1)->SetMomentum(p1)) { InError("Invalid incoming beam 1"); }
  if (!GetParticle(Particle::IncomingBeam2)->SetMomentum(p2)) { InError("Invalid incoming beam 2"); }
}

std::ostream&
operator<<(std::ostream& os, const GenericProcess::StructureFunctions& sf)
{
  switch (sf) {
    case GenericProcess::Electron:            os << "electron"; break;
    case GenericProcess::ElasticProton:       os << "elastic proton"; break;
    case GenericProcess::SuriYennie:          os << "dissociating proton [SY structure functions]"; break;
    case GenericProcess::SuriYennieLowQ2:     os << "dissociating proton [SY structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
    case GenericProcess::SzczurekUleshchenko: os << "dissociating proton [SU structure functions]"; break;
    case GenericProcess::FioreVal:            os << "dissociating proton [parton model, only valence quarks]"; break;
    case GenericProcess::FioreSea:            os << "dissociating proton [parton model, only sea quarks]"; break;
    case GenericProcess::Fiore:               os << "dissociating proton [parton model, valence and sea quarks]"; break;
  }
  return os;
}
