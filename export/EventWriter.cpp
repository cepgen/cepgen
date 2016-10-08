#include "EventWriter.h"

EventWriter::EventWriter( const OutputHandler::ExportHandler::OutputType& type, const char* filename ) :
  fFileHandler( 0 ), fType( type )
{
  switch ( fType ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC: { fFileHandler = new OutputHandler::HepMCHandler( filename ); } break;
    case OutputHandler::ExportHandler::LHE:   { fFileHandler = new OutputHandler::LHEFHandler( filename ); } break;
#endif
    default: return;
  }
}

EventWriter::~EventWriter()
{
  // HepMC persistent objects
  if ( fFileHandler ) delete fFileHandler;
}

void
EventWriter::operator<<( const Event* evt )
{
  switch ( fType ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC:
    case OutputHandler::ExportHandler::LHE: {
      (*fFileHandler) << evt; 
    } break;
#endif
    default: return;
  }
}

