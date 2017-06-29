#include "EventWriter.h"

OutputHandler::EventWriter::EventWriter( const OutputHandler::ExportHandler::OutputType& type, const char* filename ) :
  type_( type )
{
  switch ( type_ ) {
#ifdef HEPMC_LINKED
    case OutputHandler::ExportHandler::HepMC: { file_handler_ = std::unique_ptr<OutputHandler::HepMCHandler>( new OutputHandler::HepMCHandler( filename ) ); } break;
#ifdef HEPMC_VERSION3
    case OutputHandler::ExportHandler::LHE:   { file_handler_ = std::unique_ptr<OutputHandler::LHEFHandler>( new OutputHandler::LHEFHandler( filename ) ); } break;
#endif
#endif
    default: {
      std::ostringstream os; os << type_;
      Exception( __PRETTY_FUNCTION__, Form( "Unsupported output mode: %s", os.str().c_str() ), FatalError ).dump();
    }
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

