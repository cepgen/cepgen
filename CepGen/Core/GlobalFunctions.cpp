#include "CepGen/Generator.h"
#include "CepGen/Version.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/String.h"

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

#include <fstream>
#include <atomic>

#ifdef WIN32
# include <libloaderapi.h>
#else
# include <dlfcn.h>
#endif

namespace cepgen
{
  namespace utils
  {
    std::atomic<int> gSignal; ///< Abort signal handler
  }

  void
  loadLibrary( const std::string& path )
  {
#ifdef WIN32
    if ( LoadLibraryA( path.c_str() ) == nullptr )
      CG_WARNING( "loadLibrary" )
        << "Failed to load library \"" << path << "\".\n\t"
        << "Error code #" << GetLastError() << ".";
#else
    if ( dlopen( path.c_str(), RTLD_LAZY | RTLD_LOCAL ) == nullptr )
      CG_WARNING( "loadLibrary" )
        << "Failed to load library \"" << path << "\".\n\t"
        << dlerror();
#endif
  }

  void
  initialise()
  {
    //--- parse all particles properties
    static const std::string pdg_file = "External/mass_width_2019.mcd";
    pdg::MCDFileParser::parse( pdg_file.c_str() );
    //--- load all necessary modules
    loadLibrary( "CepGenAddOns/ROOTWrapper/libCepGenRoot.so" );
    loadLibrary( "CepGenAddOns/PythiaWrapper/libCepGenPythia.so" );
    loadLibrary( "CepGenAddOns/LHAPDFWrapper/libCepGenLHAPDF.so" );
    loadLibrary( "CepGenAddOns/HepMCWrapper/libCepGenHepMC.so" );
    loadLibrary( "CepGenAddOns/BoostWrapper/libCepGenBoost.so" );
    //--- header message
    try { printHeader(); } catch ( const Exception& e ) { e.dump(); }
    //--- greetings message
    CG_INFO( "init" ) << "CepGen v" << version() << " initialised. Greetings!";
  }

  void
  printHeader()
  {
    std::string tmp;
    std::ostringstream os;
    std::ifstream hf( "README" );
    if ( !hf.good() )
      throw CG_WARNING( "printHeader" ) << "Failed to open README file.";
    while ( true ) {
      if ( !hf.good() ) break;
      getline( hf, tmp );
      os << "\n " << tmp;
    }
    hf.close();
    CG_LOG( "printHeader" ) << os.str();
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
