#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>

#ifndef LIBHEPMC
# ifndef PYTHIA8
#  pragma message( "HepMC/Pythia8 are not linked to this instance! You will not able to export in LHEF." )
# endif
#else
# ifndef HEPMC_VERSION3
#  ifdef PYTHIA8
#   include "CepGen/Hadronisers/PythiaEventInterface.h"
#   include "Pythia8/Pythia.h"
namespace Pythia8 { class CepGenEvent; }
#   define PYTHIA_LHEF 1
#  else
#   pragma message( "HepMC v3 or Pythia8 are required for the LHEF export!" )
#  endif
# else
#  include "CepGen/IO/HepMCHandler.h"
#  include "HepMC/LHEF.h"
#  define HEPMC_LHEF 1
# endif
#endif

namespace cepgen
{
  namespace output
  {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        /// \param[in] filename Output file path
        explicit LHEFHandler( const char* filename );
        explicit LHEFHandler( const ParametersList& );
        ~LHEFHandler() override;

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
#if defined ( HEPMC_LHEF )
        /// Writer object (from HepMC)
        std::unique_ptr<LHEF::Writer> lhe_output_;
        LHEF::HEPRUP run_;
#elif defined ( PYTHIA_LHEF )
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::unique_ptr<Pythia8::CepGenEvent> lhaevt_;
#endif
    };

    LHEFHandler::LHEFHandler( const char* filename ) :
      GenericExportHandler( GenericExportHandler::LHE )
#if defined ( HEPMC_LHEF )
      , lhe_output_( new LHEF::Writer( filename ) )
#elif defined ( PYTHIA_LHEF )
      , pythia_( new Pythia8::Pythia ), lhaevt_( new Pythia8::CepGenEvent )
#endif
    {
#if defined ( PYTHIA_LHEF )
      lhaevt_->openLHEF( filename );
#endif
    }

    LHEFHandler::LHEFHandler( const ParametersList& params ) :
      GenericExportHandler( GenericExportHandler::LHE )
#if defined ( HEPMC_LHEF )
      , lhe_output_( new LHEF::Writer( params.get<std::string>( "filename" ) ) )
#elif defined ( PYTHIA_LHEF )
      , pythia_( new Pythia8::Pythia ), lhaevt_( new Pythia8::CepGenEvent )
#endif
    {
#if defined ( PYTHIA_LHEF )
      lhaevt_->openLHEF( params.get<std::string>( "filename" ) );
#endif
    }

    LHEFHandler::~LHEFHandler()
    {
#if defined ( PYTHIA_LHEF )
      if ( lhaevt_ )
        lhaevt_->closeLHEF( false ); // we do not want to rewrite the init block
#endif
    }

    void
    LHEFHandler::initialise( const Parameters& params )
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
#if defined ( HEPMC_LHEF )
      lhe_output_->headerBlock() << oss_init.str();
      //--- first specify information about the run
      LHEF::HEPRUP run = lhe_output_->heprup;
      run.IDBMUP = { (int)params.kinematics.incoming_beams.first.pdg, (int)params.kinematics.incoming_beams.second.pdg };
      run.EBMUP = { (double)params.kinematics.incoming_beams.first.pz, (double)params.kinematics.incoming_beams.second.pz };
      run.NPRUP = 1;
      run.resize();
      run.XSECUP[0] = params.integration().result;
      run.XERRUP[0] = params.integration().err_result;
      run.XMAXUP[0] = 1.;
      run.LPRUP[0] = 1;
      lhe_output_->heprup = run;
      //--- ensure everything is properly parsed
      lhe_output_->init();
#elif defined ( PYTHIA_LHEF )
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
#endif
    }

    void
    LHEFHandler::operator<<( const Event& ev )
    {
#if defined ( HEPMC_LHEF )
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double,double>( 0., 0. );
      out.SCALUP = 0.;
      out.AQEDUP = constants::ALPHA_EM;
      out.AQCDUP = constants::ALPHA_QCD;
      out.NUP = ev.numParticles();
      out.resize();
      for ( unsigned short ip = 0; ip < ev.numParticles(); ++ip ) {
        const Particle part = ev[ip];
        out.IDUP[ip] = part.integerPdgId(); // PDG id
        out.ISTUP[ip] = (short)part.status(); // status code
        out.PUP[ip] = part.momentum().pVector(); // momentum
        out.MOTHUP[ip] = { // mothers
          part.mothers().size() > 0 ? *part.mothers(). begin()+1 : 0,
          part.mothers().size() > 1 ? *part.mothers().rbegin()+1 : 0
        };
        out.ICOLUP[ip] = { 0, 0 };
        out.VTIMUP[ip] = 0.; // invariant lifetime
        out.SPINUP[ip] = 0.;
      }
      //lhe_output_->eventComments() << "haha";
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
#elif defined ( PYTHIA_LHEF )
      lhaevt_->feedEvent( ev, Pythia8::CepGenEvent::Type::centralAndFullBeamRemnants );
      pythia_->next();
      lhaevt_->eventLHEF();
#endif
    }

    void
    LHEFHandler::setCrossSection( double xsect, double xsect_err )
    {
#if defined ( PYTHIA_LHEF )
      lhaevt_->setCrossSection( 0, xsect, xsect_err );
#endif
    }
  }
}

REGISTER_IO_MODULE( lhef, LHEFHandler )
