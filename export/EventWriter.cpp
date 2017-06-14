#include "EventWriter.h"

OutputHandler::EventWriter::EventWriter( const OutputHandler::ExportHandler::OutputType& type, const char* filename ) :
  file_handler_( 0 ), type_( type )
{
  switch ( type_ ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC: { file_handler_ = new OutputHandler::HepMCHandler( filename ); } break;
    case OutputHandler::ExportHandler::LHE:   { file_handler_ = new OutputHandler::LHEFHandler( filename ); } break;
#endif
    default: return;
  }
}

OutputHandler::EventWriter::~EventWriter()
{
  // HepMC persistent objects
  if ( file_handler_ ) delete file_handler_;
}

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

