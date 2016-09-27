#include "EventWriter.h"

EventWriter::EventWriter( const OutputType& type, const char* filename ) :
  fType( type ),
  fHepMCHandler( 0 ), fLHEFHandler( 0 )
{
  switch ( fType ) {
#ifdef HEPMC_LINKED
    case HepMC: { fHepMCHandler = new OutputHandler::HepMCHandler( filename ); } break;
    case LHE:   { fLHEFHandler = new OutputHandler::LHEFHandler( filename ); } break;
#endif
    default: return;
  }
}

EventWriter::~EventWriter()
{
#ifdef HEPMC_LINKED
  // HepMC persistent objects
  if ( fHepMCHandler ) delete fHepMCHandler;
  if ( fLHEFHandler ) delete fLHEFHandler;
#endif
}

void
EventWriter::operator<<( const Event* evt )
{
  switch ( fType ) {
#ifdef HEPMC_LINKED
    case HepMC: { (*fHepMCHandler) << evt; } break;
    case LHE:   { (*fLHEFHandler) << evt; } break;
#endif
    default: return;
  }
}

