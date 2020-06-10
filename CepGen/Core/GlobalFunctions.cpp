#include "CepGen/Generator.h"
#include "CepGen/Version.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/AlphaS.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <fstream>
#include <atomic>
#include <utility>

#ifdef _WIN32
# include <libloaderapi.h>
#else
# include <dlfcn.h>
#endif

namespace cepgen
{
  namespace utils
  {
    const std::vector<std::string> libraries{
      "CepGenProcesses",
      "CepGenAddOns",
      "CepGenRoot",
      "CepGenPythia",
      "CepGenLHAPDF",
      "CepGenHepMC",
      "CepGenBoost",
      "CepGenRivet",
    };
    std::atomic<int> gSignal; ///< Abort signal handler
  }

  void
  loadLibrary( const std::string& path, bool match )
  {
#ifdef _WIN32
    const auto fullpath = match ? path+".dll" : path;
    if ( LoadLibraryA( fullpath.c_str() ) == nullptr )
      throw CG_WARNING( "loadLibrary" )
        << "Failed to load library \"" << path << "\".\n\t"
        << "Error code #" << GetLastError() << ".";
#else
    const auto fullpath = match ? "lib"+path+".so" : path;
    if ( dlopen( fullpath.c_str(), RTLD_LAZY | RTLD_LOCAL ) == nullptr )
      throw CG_WARNING( "loadLibrary" )
        << "Failed to load library \"" << path << "\".\n\t"
        << dlerror();
#endif
    CG_DEBUG( "loadLibrary" ) << "Loaded library \"" << path << "\".";
    loaded_libraries.emplace_back( path );
  }

  void
  initialise()
  {
    //--- parse all particles properties
    static const std::string pdg_file = "";
    search_paths = std::vector<std::string>{
      utils::environ( "CEPGEN_PATH", "." ),
      "/usr/share/CepGen"
    };

    //--- header message
    try { printHeader(); } catch ( const Exception& e ) { e.dump(); }

    //--- particles table parsing
    for ( const auto& path : search_paths )
      if ( (bool)std::ifstream( path+"/mass_width_2020.mcd" ) ) {
        pdg::MCDFileParser::parse( path+"/mass_width_2020.mcd" );
        break;
      }
    if ( PDG::get().size() < 10 )
      CG_WARNING( "init" )
        << "Only " << utils::s( "particle", PDG::get().size(), true )
        << " are defined in the runtime environment.\n\t"
        << "Make sure the path to the MCD file is correct.";

    //--- load all necessary modules
    for ( const auto& lib : utils::libraries )
      try { loadLibrary( lib, true ); } catch ( const Exception& e ) {
        e.dump(); //FIXME temporary
      }

    //--- greetings message
    CG_INFO( "init" )
      << "CepGen " << version::tag << " (" << version::extended << ") "
      << "initialised with the following add-ons:\n\t"
      << loaded_libraries << ".\n\t"
      << "Greetings!";
  }

  void
  printHeader()
  {
    for ( const auto& path : search_paths ) {
      std::ifstream hf( path+"/README" );
      if ( hf.good() ) {
        CG_LOG( "printHeader" )
          << std::string( std::istreambuf_iterator<char>( hf ),
                          std::istreambuf_iterator<char>() );
        return;
      }
    }
    throw CG_WARNING( "printHeader" ) << "Failed to open README file.";
  }

  void
  dumpModules()
  {
    const std::string sep_mid( 80, '-' );
    std::ostringstream oss;
    oss
      << "List of modules registered in the runtime database:\n";
    { oss << sep_mid << "\n"
        << utils::boldify( "Steering cards parsers" );
      if ( card::CardsHandlerFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : card::CardsHandlerFactory::get().modules() )
        oss << "\n> ." << utils::colourise( mod, utils::Colour::green )
          << " extension";
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Integration algorithms" );
      if ( IntegratorFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : IntegratorFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Physics processes" );
      if ( proc::ProcessesFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : proc::ProcessesFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green )
          << ": " << proc::ProcessesFactory::get().build( mod )->description();
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Structure functions modellings" );
      if ( strfun::StructureFunctionsFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : strfun::StructureFunctionsFactory::get().modules() )
        oss << "\n> " << utils::colourise( std::to_string( mod ), utils::Colour::green )
          << ": " << (strfun::Type)mod;
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Cross section ratios modellings" );
      if ( sigrat::SigmaRatiosFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : sigrat::SigmaRatiosFactory::get().modules() )
        oss << "\n> " << utils::colourise( std::to_string( mod ), utils::Colour::green )
          << ": " << (sigrat::Type)mod;
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Event modification modules" );
      if ( EventModifierFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : EventModifierFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Export modules" );
      if ( io::ExportModuleFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : io::ExportModuleFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "Functional evaluators" );
      if ( utils::FunctionalFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : utils::FunctionalFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    { oss << "\n" << sep_mid << "\n"
        << utils::boldify( "alpha(s) evolution algorithms" );
      if ( AlphaSFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : AlphaSFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    CG_INFO( "dumpModules" ) << oss.str();
  }
}
