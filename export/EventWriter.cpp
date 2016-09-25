#include "EventWriter.h"

EventWriter::EventWriter( const OutputType& type, const char* filename ) :
  fType( type ),
  fHepMCHandler( 0 ), fLHEFHandler( 0 )
{
  switch ( fType ) {
    case HepMC: { fHepMCHandler = new OutputHandler::HepMCHandler( filename ); } break;
    case LHE:   { fLHEFHandler = new OutputHandler::LHEFHandler( filename ); } break;
    default: return;
  }
}

EventWriter::~EventWriter()
{
  // HepMC persistent objects
  if ( fHepMCHandler ) delete fHepMCHandler;
  if ( fLHEFHandler ) delete fLHEFHandler;
}

void
EventWriter::operator<<( const Event* evt )
{
  switch ( fType ) {
    case HepMC: { (*fHepMCHandler) << evt; } break;
    case LHE:   { (*fLHEFHandler) << evt; } break;
    default: return;
  }
}

