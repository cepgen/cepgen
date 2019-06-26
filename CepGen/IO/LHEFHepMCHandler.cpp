#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>

#include "HepMC/Version.h"
#ifndef HEPMC_VERSION_CODE // HepMC v2
#  error "HepMC v3 is required for the LHEF export!"
#else // HepMC v3+
#  include "HepMC/LHEF.h"
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
    class LHEFHepMCHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        /// \param[in] filename Output file path
        explicit LHEFHepMCHandler( const char* filename );
        explicit LHEFHepMCHandler( const ParametersList& );
        ~LHEFHepMCHandler() override;

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override {}

      private:
        /// Writer object (from HepMC)
        std::unique_ptr<LHEF::Writer> lhe_output_;
        LHEF::HEPRUP run_;
    };

    LHEFHepMCHandler::LHEFHepMCHandler( const char* filename ) :
      GenericExportHandler( GenericExportHandler::LHE ),
      lhe_output_( new LHEF::Writer( filename ) )
    {}

    LHEFHepMCHandler::LHEFHepMCHandler( const ParametersList& params ) :
      GenericExportHandler( GenericExportHandler::LHE ),
      lhe_output_( new LHEF::Writer( params.get<std::string>( "filename" ) ) )
    {}

    void
    LHEFHepMCHandler::initialise( const Parameters& params )
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
    }

    void
    LHEFHepMCHandler::operator<<( const Event& ev )
    {
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
    }
  }
}

REGISTER_IO_MODULE( lhef, LHEFHepMCHandler )
