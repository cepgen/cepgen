#include "CepGen/IO/ExportHandler.h"

#include <iostream>

namespace cepgen
{
  namespace output
  {
    ExportHandler::ExportHandler( const OutputType& type ) :
      type_( type ), event_num_( 0. )
    {}
  }

  std::ostream&
  operator<<( std::ostream& os, const output::ExportHandler::OutputType& type )
  {
    switch ( type ) {
      case output::ExportHandler::HepMC:
        return os << "HepMC ASCII";
      case output::ExportHandler::LHE:
        return os << "LHEF";
    }
    return os;
  }
}

