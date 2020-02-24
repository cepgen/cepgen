#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/ExportModule.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Limits.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

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

        const Args argv_;
    };

    CommandLineHandler::CommandLineHandler( const ParametersList& params ) :
      argv_( params.get<std::vector<std::string> >( "args" ) )
    {
      auto pars = vectorise( argv_ );
      const auto& proc = pars.get<ParametersList>( "process" );
      if ( proc.empty() )
        throw CG_FATAL( "CommandLineHandler" )
          << "Failed to retrieve a process in the configuration!";
      params_.setProcess( proc::ProcessesFactory::get().build( proc ) );
      params_.kinematics.mode = (KinematicsMode)proc.get<int>( "mode", 1 );
      const auto& kin = pars.get<ParametersList>( "kinematics" );
      params_.kinematics.cuts.central.pt_single.min() = kin.get<double>( "ptmin", 0. );
      params_.kinematics.cuts.central.pt_single.max() = kin.get<double>( "ptmax", 9999. );
      const auto& gen = pars.get<ParametersList>( "generation" );
      params_.generation().maxgen = (unsigned long)gen.get<int>( "ngen", 10000 );
      params_.generation().enabled = params_.generation().maxgen > 1;
      params_.generation().num_threads = gen.get<int>( "nthreads", 2 );
      params_.generation().gen_print_every = gen.get<int>( "nprn", 1000 );
      params_.generation().treat = gen.get<bool>( "treat", true );
      params_.integration().rng_seed = gen.get<int>( "seed" );
      const auto& out = pars.get<ParametersList>( "output" );
      if ( !out.empty() )
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
            if ( strtok( (char*)word.at( 1 ).data(), "." ) )
              plist.set<double>( key, std::stod( word.at( 1 ) ) );
            else
              plist.set<int>( key, std::stod( word.at( 1 ) ) );
          } catch ( const std::invalid_argument& ) {
            plist.set<std::string>( key, word.at( 1 ) );
          }
        }
        else // sub-parameters word found
          plist += vectorise( Args( 1,
            utils::merge( std::vector<std::string>( word.begin()+1, word.end() ), "=" ) ) );
      }
      CG_INFO("")<<params;
      return params;
    }
  }
}

REGISTER_CARD_HANDLER( "cmd", CommandLineHandler )
