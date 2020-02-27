#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/ExportModule.h"

#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/EventModifier.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Limits.h"

#include "CepGen/Core/Integrator.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <fstream>

namespace cepgen
{
  namespace card
  {
    /// Command line parser
    class CommandLineHandler : public Handler
    {
      public:
        explicit CommandLineHandler( const ParametersList& );

      private:
        typedef std::vector<std::string> Args;

        static ParametersList vectorise( const Args& );
        static const double INVALID;

        const std::string filename_;
        Args argv_;
    };

    const double CommandLineHandler::INVALID = -999.999;

    CommandLineHandler::CommandLineHandler( const ParametersList& params ) :
      filename_( params.get<std::string>( FILENAME_KEY ) ),
      argv_( params.get<std::vector<std::string> >( "args" ) )
    {
      if ( !filename_.empty() ) {
        std::ifstream file( filename_ );
        std::string line;
        while ( getline( file, line ) )
          argv_.emplace_back( line );
        file.close();
      }
      const auto pars = vectorise( argv_ );

      //----- process definition
      const auto& proc = pars.get<ParametersList>( "process" );
      if ( proc.empty() )
        throw CG_FATAL( "CommandLineHandler" )
          << "Failed to retrieve a process in the configuration!";
      params_.setProcess( proc::ProcessesFactory::get().build( proc ) );

      //----- structure functions
      auto strfun = pars.get<ParametersList>( "strfun" ); // copy
      if ( strfun.name<int>( -999 ) == -999 )
        strfun.setName<int>( 11 ); // default is Suri-Yennie
      params_.kinematics.structure_functions = strfun::StructureFunctionsFactory::get().build( strfun );

      //----- phase space definition
      const auto& kin = pars.get<ParametersList>( "kinematics" );
      params_.kinematics.mode = (KinematicsMode)kin.get<int>( "mode", 1 );
      //--- incoming beams
      params_.kinematics.incoming_beams.first.pdg = kin.get<int>( "beam1id", 2212 );
      params_.kinematics.incoming_beams.first.pz = kin.get<double>( "beam1pz", 6500. );
      const int hi_A1 = kin.get<int>( "beam1A", 1 ), hi_Z1 = kin.get<int>( "beam1Z", 0 );
      if ( hi_Z1 != 0 )
        params_.kinematics.incoming_beams.first.pdg = HeavyIon( hi_A1, (Element)hi_Z1 );
      params_.kinematics.incoming_beams.second.pdg = kin.get<int>( "beam2id", 2212 );
      params_.kinematics.incoming_beams.second.pz = kin.get<double>( "beam2pz", 6500. );
      const int hi_A2 = kin.get<int>( "beam2A", 1 ), hi_Z2 = kin.get<int>( "beam2Z", 0 );
      if ( hi_Z2 != 0 )
        params_.kinematics.incoming_beams.second.pdg = HeavyIon( hi_A2, (Element)hi_Z2 );
      const double sqrt_s = kin.get<double>( "sqrts", INVALID );
      if ( sqrt_s != INVALID )
        params_.kinematics.setSqrtS( sqrt_s );
      //--- incoming partons
      params_.kinematics.cuts.central.q2.min() = kin.get<double>( "q2min", INVALID );
      params_.kinematics.cuts.central.q2.max() = kin.get<double>( "q2max", INVALID );
      params_.kinematics.cuts.central.qt.min() = kin.get<double>( "qtmin", INVALID );
      params_.kinematics.cuts.central.qt.max() = kin.get<double>( "qtmax", INVALID );
      //--- outgoing remnants
      params_.kinematics.cuts.remnants.mass_single.min() = kin.get<double>( "mxmin", INVALID );
      params_.kinematics.cuts.remnants.mass_single.max() = kin.get<double>( "mxmax", INVALID );
      const Limits xi_rng = { kin.get<double>( "ximin", 0. ), kin.get<double>( "ximax", 1. ) };
      if ( xi_rng.valid() )
        params_.kinematics.cuts.remnants.energy_single = -( xi_rng-1. )*0.5*params_.kinematics.sqrtS();
      //--- central system
      params_.kinematics.cuts.central.pt_single.min() = kin.get<double>( "ptmin", INVALID );
      params_.kinematics.cuts.central.pt_single.max() = kin.get<double>( "ptmax", INVALID );
      params_.kinematics.cuts.central.eta_single.min() = kin.get<double>( "etamin", INVALID );
      params_.kinematics.cuts.central.eta_single.max() = kin.get<double>( "etamax", INVALID );
      params_.kinematics.cuts.central.rapidity_single.min() = kin.get<double>( "rapmin", INVALID );
      params_.kinematics.cuts.central.rapidity_single.max() = kin.get<double>( "rapmax", INVALID );
      params_.kinematics.cuts.central.pt_sum.min() = kin.get<double>( "ptsummin", INVALID );
      params_.kinematics.cuts.central.pt_sum.max() = kin.get<double>( "ptsummax", INVALID );
      params_.kinematics.cuts.central.mass_sum.min() = kin.get<double>( "msummin", INVALID );
      params_.kinematics.cuts.central.mass_sum.max() = kin.get<double>( "msummax", INVALID );

      //----- integration
      const auto& integ = pars.get<ParametersList>( "integrator" );
      const std::string integ_name = integ.name<std::string>( "vegas" );
      if ( integ_name == "plain" )
        params_.integration().type = IntegratorType::plain;
      else if ( integ_name == "vegas" )
        params_.integration().type = IntegratorType::Vegas;
      else if ( integ_name == "miser" )
        params_.integration().type = IntegratorType::MISER;

      //----- events generation
      const auto& gen = pars.get<ParametersList>( "generation" );
      params_.generation().maxgen = (unsigned long)gen.get<int>( "ngen", 10000 );
      params_.generation().enabled = params_.generation().maxgen > 1;
      params_.generation().num_threads = gen.get<int>( "nthreads", 2 );
      params_.generation().gen_print_every = gen.get<int>( "nprn", 1000 );
      params_.generation().treat = gen.get<bool>( "treat", true );
      params_.integration().rng_seed = gen.get<int>( "seed" );

      //----- event modification modules
      const auto& mod = pars.get<ParametersList>( "eventmod" );
      if ( !mod.keys().empty() ) {
        params_.addModifier( EventModifierFactory::get().build( mod ) );
        (*params_.eventModifiersSequence().rbegin())->setParameters( params_ );
      }

      //----- output modules definition
      const auto& out = pars.get<ParametersList>( "output" );
      if ( !out.keys().empty() )
        params_.addOutputModule( io::ExportModuleFactory::get().build( out ) );
    }

    ParametersList
    CommandLineHandler::vectorise( const Args& args )
    {
      ParametersList params;
      for ( const auto& arg : args ) {
        auto cmd = utils::split( arg, ':' );
        auto& plist = cmd.size() < 2
          ? params
          : params.operator[]<ParametersList>( cmd.at( 0 ) );
        if ( cmd.size() > 2 ){ // sub-parameters word found
          plist += vectorise( Args( 1,
            utils::merge( std::vector<std::string>( cmd.begin()+1, cmd.end() ), ":" ) ) );
          continue;
        }
        const auto word = cmd.size() < 2
          ? utils::split( cmd.at( 0 ), '=' )
          : utils::split( cmd.at( 1 ), '=' );
        auto key = word.at( 0 );
        if ( key == "name" )
          key = ParametersList::MODULE_NAME;
        if ( word.size() == 1 ) // basic key=true
          plist.set<bool>( key, true );
        else if ( word.size() == 2 ) { // basic key=value
          try {
            if ( word.at( 1 ).find( '.' ) != std::string::npos )
              // double if contains a '.'
              plist.set<double>( key, std::stod( word.at( 1 ) ) );
            else
              plist.set<int>( key, std::stod( word.at( 1 ) ) );
          } catch ( const std::invalid_argument& ) {
            plist.set<std::string>( key, word.at( 1 ) );
          }
        }
      }
      CG_INFO("")<<params;
      return params;
    }
  }
}

REGISTER_CARD_HANDLER( "cmd", CommandLineHandler )
