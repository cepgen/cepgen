#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Hadronisers/PythiaEventInterface.h"
#include "Pythia8/Pythia.h"

#include <sstream>

namespace cepgen
{
  namespace output
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
        /// \param[in] filename Output file path
        explicit LHEFPythiaHandler( const char* filename );
        explicit LHEFPythiaHandler( const ParametersList& );
        ~LHEFPythiaHandler() override;

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::unique_ptr<Pythia8::CepGenEvent> lhaevt_;
    };

    LHEFPythiaHandler::LHEFPythiaHandler( const char* filename ) :
      GenericExportHandler( GenericExportHandler::LHE ),
      pythia_( new Pythia8::Pythia ), lhaevt_( new Pythia8::CepGenEvent )
    {
      lhaevt_->openLHEF( filename );
    }

    LHEFPythiaHandler::LHEFPythiaHandler( const ParametersList& params ) :
      GenericExportHandler( GenericExportHandler::LHE ),
      pythia_( new Pythia8::Pythia ), lhaevt_( new Pythia8::CepGenEvent )
    {
      lhaevt_->openLHEF( params.get<std::string>( "filename" ) );
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
        << "<!--\n"
        << "  ***** Sample generated with CepGen v" << version() << " *****\n"
        << "  * process: " << params.processName() << " (" << params.kinematics.mode << ")\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic ) {
        oss_init << "  * structure functions: " << params.kinematics.structure_functions->type << "\n";
        if ( !params.hadroniserName().empty() )
          oss_init << "  * hadroniser: " << params.hadroniserName() << "\n";
      }
      oss_init
        << "  *--- incoming state\n";
      if ( params.kinematics.cuts.initial.q2.valid() )
        oss_init
          << "  * Q2 range (GeV2): "
          << params.kinematics.cuts.initial.q2 << "\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic
        && params.kinematics.cuts.remnants.mass_single.valid() )
        oss_init
          << "  * remnants mass range (GeV/c2): "
          << params.kinematics.cuts.remnants.mass_single << "\n";
      oss_init << "  *--- central system\n";
      if ( params.kinematics.cuts.central.pt_single.valid() )
        oss_init
          << "  * single particle pt (GeV/c): "
          << params.kinematics.cuts.central.pt_single << "\n";
      if ( params.kinematics.cuts.central.energy_single.valid() )
        oss_init
          << "  * single particle energy (GeV): "
          << params.kinematics.cuts.central.energy_single << "\n";
      if ( params.kinematics.cuts.central.eta_single.valid() )
        oss_init
          << "  * single particle eta: "
          << params.kinematics.cuts.central.eta_single << "\n";
      if ( params.kinematics.cuts.central.pt_sum.valid() )
        oss_init
          << "  * total pt (GeV/c): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      if ( params.kinematics.cuts.central.mass_sum.valid() )
        oss_init
          << "  * total invariant mass (GeV/c2): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      oss_init
        << "  **************************************************\n"
        << "-->";
      oss_init << std::endl; // LHEF is usually not as beautifully parsed as a standard XML...
                             // we're physicists, what do you expect?
      lhaevt_->addComments( oss_init.str() );
      lhaevt_->initialise( params );
      pythia_->settings.mode( "Beams:frameType", 5 ); // LHEF event readout
      pythia_->settings.mode( "Next:numberCount", 0 ); // remove some of the Pythia output
      pythia_->settings.flag( "ProcessLevel:all", false ); // we do not want Pythia to interfere...
      pythia_->setLHAupPtr( lhaevt_.get() );
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
