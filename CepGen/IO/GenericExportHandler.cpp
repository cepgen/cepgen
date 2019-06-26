#include "CepGen/IO/GenericExportHandler.h"

#include <iostream>

namespace cepgen
{
  namespace output
  {
    GenericExportHandler::GenericExportHandler( const OutputType& type ) :
      type_( type ), event_num_( 0. )
    {}
  }

  std::ostream&
  operator<<( std::ostream& os, const output::GenericExportHandler::OutputType& type )
  {
    switch ( type ) {
      case output::GenericExportHandler::HepMC:
        return os << "HepMC ASCII";
      case output::GenericExportHandler::LHE:
        return os << "LHEF";
      case output::GenericExportHandler::DOT:
        return os << "DOT graphics";
    }
    return os;
  }
}

