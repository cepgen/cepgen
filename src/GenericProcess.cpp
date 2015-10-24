#include "GenericProcess.h"

GenericProcess::GenericProcess(std::string name_) :
  fX(0), fNumDimensions(0), fIsPointSet(false),
  fIsInStateSet(false), fIsOutStateSet(false), fIsKinematicSet(false),
  fEvent(new Event), fName(name_)
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
GenericProcess::DumpPoint(const ExceptionType& et=Info)
{
  std::ostringstream os;
  for (unsigned int i=0; i<(unsigned int)fNumDimensions; i++) {
    os << Form("  x(%2d) = %8.6f\n\t", i, fX[i]);
  }
  if (et<Debugging) { Info(Form("Number of integration parameters: %d\n\t"
                                "%s", fNumDimensions, os.str().c_str())); }
  else           {   Debug(Form("Number of integration parameters: %d\n\t"
                                "%s", fNumDimensions, os.str().c_str())); }
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
