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

        static ParametersList vectorise( const Args& );
        static const double INVALID;

        Args argv_;
    };

    const double CommandLineHandler::INVALID = -999.999;

    CommandLineHandler::CommandLineHandler( const ParametersList& params ) :
      argv_( params.get<std::vector<std::string> >( "args" ) )
    {
      const auto filename = params.get<std::string>( FILENAME_KEY );
      if ( !filename.empty() )
        parse( filename, params_ );
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

      const auto pars = vectorise( argv_ );

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

      //----- structure functions
      auto strfun = pars.get<ParametersList>( "strfun" ); // copy
      if ( !strfun.empty() || !params_->kinematics.structure_functions ) {
        if ( strfun.name<int>( -999 ) == -999 )
          strfun.setName<int>( 11 ); // default is Suri-Yennie
        params_->kinematics.structure_functions = strfun::StructureFunctionsFactory::get().build( strfun );
      }

      //----- phase space definition
      params_->kinematics = Kinematics( pars.get<ParametersList>( "kinematics" ) );

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
    CommandLineHandler::vectorise( const Args& args )
    {
      const char delim_block = ':', delim_eq1 = '[', delim_eq2 = ']';
      ParametersList params;
      for ( const auto& arg : args ) {
        auto cmd = utils::split( arg, delim_block );
        auto& plist = cmd.size() < 2
          ? params
          : params.operator[]<ParametersList>( cmd.at( 0 ) );
        if ( cmd.size() > 2 ){ // sub-parameters word found
          plist += vectorise( Args( 1,
            utils::merge( std::vector<std::string>( cmd.begin()+1, cmd.end() ), std::string( delim_block, 1 ) ) ) );
          continue;
        }
        const auto word = cmd.size() < 2 ? cmd.at( 0 ) : cmd.at( 1 );
        auto words = utils::split( word, delim_eq1 );
        auto key = words.at( 0 );
        if ( key == "name" )
          key = ParametersList::MODULE_NAME;
        if ( words.size() == 1 ) // basic key=true
          plist.set<bool>( key, true );
        else if ( words.size() == 2 ) { // basic key=value
          words = utils::split( words.at( 1 ), delim_eq2 );
          if ( words.size() != 1 )
            throw CG_FATAL( "CommandLineHandler" ) << "Invalid syntax for key \"" << key << "\"!";
          const auto value = words.at( 0 );
          try {
            if ( value.find( '.' ) != std::string::npos
              || value.find( 'e' ) != std::string::npos
              || value.find( 'E' ) != std::string::npos )
              // double if contains a '.'/'e'
              plist.set<double>( key, std::stod( value ) );
            else
              plist.set<int>( key, std::stod( value ) );
          } catch ( const std::invalid_argument& ) {
            if ( value == "off" || value == "no" )
              plist.set<bool>( key, false );
            else if ( value == "on" || value == "yes" )
              plist.set<bool>( key, true );
            else
              plist.set<std::string>( key, value );
          }
        }
      }
      CG_INFO("")<<params;
      return params;
    }
  }
}

REGISTER_CARD_HANDLER( "cmd", CommandLineHandler )
