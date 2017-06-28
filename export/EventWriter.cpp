#include "EventWriter.h"

OutputHandler::EventWriter::EventWriter( const OutputHandler::ExportHandler::OutputType& type, const char* filename ) :
  type_( type )
{
  switch ( type_ ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC: { file_handler_ = std::make_unique<OutputHandler::HepMCHandler>( filename ); } break;
    case OutputHandler::ExportHandler::LHE:   { file_handler_ = std::make_unique<OutputHandler::LHEFHandler>( filename ); } break;
#endif
    default: return;
  }
}

OutputHandler::EventWriter::~EventWriter()
{}

void
OutputHandler::EventWriter::operator<<( const Event* evt )
{
  switch ( type_ ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC:
    case OutputHandler::ExportHandler::LHE: {
      ( *file_handler_ ) << evt;
    } break;
#endif
    default: return;
  }
}

