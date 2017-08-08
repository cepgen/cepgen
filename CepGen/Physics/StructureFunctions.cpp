#include "StructureFunctions.h"

namespace CepGen
{
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    switch ( sf ) {
      case Electron:            os << "electron"; break;
      case ElasticProton:       os << "elastic proton"; break;
      case SuriYennie:          os << "Suri-Yennie (diss.proton)"; break;
      case SuriYennieLowQ2:     os << "Suri-Yennie (diss.proton, MX < 2 GeV, Q^2 < 5 GeV^2)"; break;
      case SzczurekUleshchenko: os << "Szczurek-Uleshchenko (diss.proton)"; break;
      case FioreVal:            os << "Fiore (diss.proton, only valence quarks)"; break;
      case FioreSea:            os << "Fiore (diss.proton, only sea quarks)"; break;
      case Fiore:               os << "Fiore (diss.proton, valence and sea quarks)"; break;
    }
    return os;
  }
}
