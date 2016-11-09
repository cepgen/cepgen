#include "StructureFunctions.h"

std::ostream&
operator<<( std::ostream& os, const StructureFunctions& sf )
{
  switch (sf) {
    case Electron:            os << "electron"; break;
    case ElasticProton:       os << "elastic proton"; break;
    case SuriYennie:          os << "dissociating proton [SY structure functions]"; break;
    case SuriYennieLowQ2:     os << "dissociating proton [SY structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
    case SzczurekUleshchenko: os << "dissociating proton [SU structure functions]"; break;
    case FioreVal:            os << "dissociating proton [parton model, only valence quarks]"; break;
    case FioreSea:            os << "dissociating proton [parton model, only sea quarks]"; break;
    case Fiore:               os << "dissociating proton [parton model, valence and sea quarks]"; break;
  }
  return os;
}
