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
  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop)
    DumpPoint(Debugging);
}

void
GenericProcess::ClearEvent()
{
  fEvent->Restore();
  //AddEventContent();
  AddEventKinematics();
}

void
GenericProcess::PrepareKinematics()
{
  Particle *ip1 = fEvent->GetOneByRole(IncomingBeam1), *ip2 = fEvent->GetOneByRole(IncomingBeam2);
  double *p1 = ip1->P4(), *p2 = ip2->P4(), k = 0.;

  for (unsigned int i=0; i<3; i++) k += p1[i]*p2[i];
  fS = ip1->M2()+ip2->M2()+2.*(ip1->E()*ip2->E()-k);
  fSqS = sqrt(fS);
}

void
GenericProcess::DumpPoint(const ExceptionType& et=Info)
{
  std::ostringstream os;
  for (unsigned int i=0; i<(unsigned int)fNumDimensions; i++) {
    os << Form("  x(%2d) = %8.6f\n\t", i, fX[i]);
  }
  if (et<Debugging) { Info(Form("Number of integration parameters: %d\n\t"
                                "%s", fNumDimensions, os.str().c_str())); }
  else              { Debug(Form("Number of integration parameters: %d\n\t"
                                 "%s", fNumDimensions, os.str().c_str())); }
}

void
GenericProcess::SetEventContent(IncomingState is, OutgoingState os)
{
  bool has_cs = false;
  // Incoming particles (incl. eventual partons)
  for (IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++) {
    fEvent->AddParticle(Particle(ip->first, ip->second));
    Particle* p = fEvent->GetOneByRole(ip->first);
    //p->SetM(ip-second);
    p->status = 0;
    switch (ip->first) {
      case IncomingBeam1: case IncomingBeam2: break;
      case Parton1: p->SetMother(fEvent->GetOneByRole(IncomingBeam1)); break;
      case Parton2: p->SetMother(fEvent->GetOneByRole(IncomingBeam2)); break;
      case CentralSystem: p->SetMother(fEvent->GetOneByRole(Parton1)); has_cs = true; break;
      default: break;
    }
  }
  // Prepare the central system if not already there
  if (!has_cs) {
    Particle* moth = fEvent->GetOneByRole(Parton1);
    Particle cs(4, moth->GetPDGId());
    cs.SetMother(moth);
    fEvent->AddParticle(cs);
  }
  // Outgoing particles (central, and outgoing primary particles or remnants)
  for (OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++) {
    fEvent->AddParticle(Particle(op->first, op->second));
    Particle* p = fEvent->GetOneByRole(op->first);
    //p->SetM(ip->second);
    p->status = 0;
    switch (op->first) {
      case OutgoingBeam1: p->SetMother(fEvent->GetOneByRole(IncomingBeam1)); break;
      case OutgoingBeam2: p->SetMother(fEvent->GetOneByRole(IncomingBeam2)); break;
      case CentralParticle1: p->SetMother(fEvent->GetOneByRole(CentralSystem)); break;
      case CentralParticle2: p->SetMother(fEvent->GetOneByRole(CentralSystem)); break;
      default: break;
    }
  }
  fEvent->Init();
}

std::ostream&
operator<<(std::ostream& os, const GenericProcess::ProcessMode& pm)
{
  switch (pm) {
    case GenericProcess::ElasticElastic:      os << "Elastic/Elastic"; break;
    case GenericProcess::InelasticElastic:    os << "Inelastic/Elastic"; break;
    case GenericProcess::ElasticInelastic:    os << "Elastic/Inelastic"; break;
    case GenericProcess::InelasticInelastic:  os << "Inelastic/Inelastic"; break;    
  }
  return os;
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
