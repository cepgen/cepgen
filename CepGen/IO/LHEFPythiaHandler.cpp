#include "CepGen/IO/ExportHandler.h"
#include "CepGen/IO/PythiaEventInterface.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include "Pythia8/Pythia.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFPythiaHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        explicit LHEFPythiaHandler( const ParametersList& );
        ~LHEFPythiaHandler();

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::unique_ptr<Pythia8::CepGenEvent> lhaevt_;
    };

    LHEFPythiaHandler::LHEFPythiaHandler( const ParametersList& params ) :
      GenericExportHandler( "lhef" ),
      pythia_( new Pythia8::Pythia ), lhaevt_( new Pythia8::CepGenEvent )
    {
      lhaevt_->openLHEF( params.get<std::string>( "filename", "output.lhe" ) );
    }

    LHEFPythiaHandler::~LHEFPythiaHandler()
    {
      if ( lhaevt_ )
        lhaevt_->closeLHEF( false ); // we do not want to rewrite the init block
    }

    void
    LHEFPythiaHandler::initialise( const Parameters& params )
    {
      std::ostringstream oss_init;
      oss_init
        << "<!--\n" << banner( params ) << "\n-->";
      oss_init << std::endl; // LHEF is usually not as beautifully parsed as a standard XML...
                             // we're physicists, what do you expect?
      lhaevt_->addComments( oss_init.str() );
      lhaevt_->initialise( params );
      pythia_->setLHAupPtr( lhaevt_.get() );
      pythia_->settings.flag( "ProcessLevel:all", false ); // we do not want Pythia to interfere...
      pythia_->settings.flag( "PartonLevel:all", false ); // we do not want Pythia to interfere...
      pythia_->settings.flag( "HadronLevel:all", false ); // we do not want Pythia to interfere...
      pythia_->settings.mode( "Beams:frameType", 5 ); // LHEF event readout
      pythia_->settings.mode( "Next:numberCount", 0 ); // remove some of the Pythia output
      //pythia_->settings.flag( "Check:event", false );
      pythia_->init();
      lhaevt_->initLHEF();
    }

    void
    LHEFPythiaHandler::operator<<( const Event& ev )
    {
      lhaevt_->feedEvent( ev, Pythia8::CepGenEvent::Type::centralAndFullBeamRemnants );
      pythia_->next();
      lhaevt_->eventLHEF();
    }

    void
    LHEFPythiaHandler::setCrossSection( double xsect, double xsect_err )
    {
      lhaevt_->setCrossSection( 0, xsect, xsect_err );
    }
  }
}

REGISTER_IO_MODULE( lhef, LHEFPythiaHandler )
