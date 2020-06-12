#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Limits.h"

#include "CepGen/Integration/Integrator.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

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

        Parameters* parse( const std::string&, Parameters* ) override;

      private:
        typedef std::vector<std::string> Args;

        static ParametersList vectorise( const std::string& );
        static const double INVALID;

        Args argv_;
    };

    const double CommandLineHandler::INVALID = -999.999;

    CommandLineHandler::CommandLineHandler( const ParametersList& params ) :
      Handler( params ),
      argv_( params.get<std::vector<std::string> >( "args" ) )
    {
      if ( !filename_.empty() )
        parse( filename_, params_ );
    }

    Parameters*
    CommandLineHandler::parse( const std::string& filename, Parameters* params )
    {
      if ( !filename.empty() ) {
        std::ifstream file( filename );
        std::string line;
        while ( getline( file, line ) )
          argv_.emplace_back( line );
        file.close();
      }

      ParametersList pars;
      for ( const auto& arg : argv_ )
        pars += vectorise( arg );
      CG_INFO( "CommandLineHandler" )
        << "Arguments list: " << argv_ << " unpacked to:\n\t"
        << pars << ".";

      params_ = params;

      //----- timer definition
      if ( pars.get<bool>( "timer", false ) )
        params_->setTimeKeeper( new utils::TimeKeeper );

      //----- process definition
      auto proc = pars.get<ParametersList>( "process" );
      if ( !proc.empty() ) {
        if ( params_->hasProcess() )
          proc = ParametersList( params_->process().parameters() )+proc;
        params_->setProcess( proc::ProcessesFactory::get().build( proc ) );
      }

      //----- phase space definition
      auto kin = pars.get<ParametersList>( "kinematics" )
        .set<ParametersList>( "structureFunctions", pars.get<ParametersList>( "strfun" ) );
      params_->kinematics = Kinematics( params_->kinematics.parameters()+kin );

      //----- integration
      pars.fill<ParametersList>( "integrator", *params_->integrator );

      //----- events generation
      const auto& gen = pars.get<ParametersList>( "generation" );
      params_->generation().maxgen = (unsigned long)gen.get<int>( "ngen", 10000 );
      params_->generation().enabled = params_->generation().maxgen > 1;
      if ( gen.has<int>( "nthreads" ) )
        params_->generation().num_threads = gen.get<int>( "nthreads" );
      if ( gen.has<int>( "nprn" ) )
        params_->generation().gen_print_every = gen.get<int>( "nprn" );
      if ( gen.has<int>( "seed" ) )
        params_->integrator->set<int>( "seed", gen.get<int>( "seed" ) );

      //----- event modification modules
      const auto& mod = pars.get<ParametersList>( "eventmod" );
      if ( !mod.keys( true ).empty() )
        params_->addModifier( EventModifierFactory::get().build( mod ) );

      //----- output modules definition
      const auto& out = pars.get<ParametersList>( "output" );
      if ( !out.keys( true ).empty() )
        params_->addOutputModule( io::ExportModuleFactory::get().build( out ) );
      return params_;
    }

    ParametersList
    CommandLineHandler::vectorise( const std::string& arg )
    {
      ParametersList params;
      auto cmd = utils::split( arg, '/' );
      auto& plist = cmd.size() > 1
        ? params.operator[]<ParametersList>( cmd.at( 0 ) )
        : params;
      if ( cmd.size() > 1 ) // sub-parameters word found
        plist += vectorise( utils::merge(
          std::vector<std::string>( cmd.begin()+1, cmd.end() ), "/" ) );
      else {
        const auto word = cmd.at( 0 );
        auto words = utils::split( word, '=' );
        auto key = words.at( 0 );
        if ( key == "name" )
          key = ParametersList::MODULE_NAME;
        if ( words.size() == 1 ) // basic key=true
          plist.set<bool>( key, true );
        else if ( words.size() == 2 ) { // basic key=value
          const auto value = words.at( 1 );
          try {
            if ( value.find( '.' ) != std::string::npos
              || value.find( 'e' ) != std::string::npos
              || value.find( 'E' ) != std::string::npos )
              // double if contains a '.'/'e'
              plist.set<double>( key, std::stod( value ) );
            else
              plist.set<int>( key, std::stod( value ) );
          } catch ( const std::invalid_argument& ) {
            if ( value == "off" || value == "no" || value == "false" )
              plist.set<bool>( key, false );
            else if ( value == "on" || value == "yes" || value == "true" )
              plist.set<bool>( key, true );
            else
              plist.set<std::string>( key, value );
          }
        }
        else
          throw CG_FATAL( "CommandLineHandler" )
            << "Invalid key=value unpacking: " << word << "!";
      }
      return params;
    }
  }
}

REGISTER_CARD_HANDLER( "cmd", CommandLineHandler )
