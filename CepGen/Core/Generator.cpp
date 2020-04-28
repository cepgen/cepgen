#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Timer.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/GridParameters.h"

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

#include "CepGen/Event/Event.h"

#include <fstream>
#include <chrono>
#include <atomic>

namespace cepgen
{
  namespace utils
  {
    std::atomic<int> gSignal; ///< Abort signal handler
  }

  Generator::Generator() :
    parameters_( new Parameters ),
    result_( -1. ), result_error_( -1. )
  {
    static const std::string pdg_file = "External/mass_width_2019.mcd";
    pdg::MCDFileParser::parse( pdg_file.c_str() );
    CG_DEBUG( "Generator:init" ) << "Generator initialized";
    try {
      printHeader();
    } catch ( const Exception& e ) {
      e.dump();
    }
    //--- random number initialization
    std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
    srandom( time.time_since_epoch().count() );
  }

  Generator::Generator( Parameters* ip ) :
    parameters_( ip ),
    result_( -1. ), result_error_( -1. )
  {}

  Generator::~Generator()
  {}

  void
  Generator::clearRun()
  {
    generator_.reset( new GeneratorWorker( parameters_.get() ) );
    result_ = result_error_ = -1.;
    parameters_->prepareRun();
    tmr_.clear();
  }

  Parameters&
  Generator::parameters()
  {
    return *parameters_;
  }

  void
  Generator::setParameters( Parameters ip )
  {
    parameters_.reset( new Parameters( ip ) ); // copy constructor
    if ( parameters_->hasProcess() )
      parameters_->process().setKinematics( parameters_->kinematics );
  }

  void
  Generator::printHeader() const
  {
    std::string tmp;
    std::ostringstream os; os << "version " << version() << std::endl;
    std::ifstream hf( "README" );
    if ( !hf.good() )
      throw CG_WARNING( "Generator" ) << "Failed to open README file.";
    while ( true ) {
      if ( !hf.good() ) break;
      getline( hf, tmp );
      os << "\n " << tmp;
    }
    hf.close();
    CG_LOG( "Generator" ) << os.str();
  }

  double
  Generator::computePoint( double* x )
  {
    clearRun();
    if ( !parameters_->hasProcess() )
      throw CG_FATAL( "Generator:computePoint" )
        << "Trying to compute a point with no process specified!";
    const size_t ndim = parameters_->process().ndim();
    if ( ndim < 1 )
      throw CG_FATAL( "Generator:computePoint" )
        << "Invalid phase space dimension (ndim=" << ndim << ")!";
    std::vector<double> coord( x, x+ndim );
    double res = generator_->integrand().eval( coord );
    CG_DEBUG( "Generator:computePoint" )
      << "Result for x[" << ndim << "] = " << coord << ":\n\t" << res << ".";
    return res;
  }

  void
  Generator::computeXsection( double& xsec, double& err )
  {
    CG_INFO( "Generator" )
      << "Starting the computation of the process cross-section.";

    integrate();

    xsec = result_;
    err = result_error_;

    if ( xsec < 1.e-2 )
      CG_INFO( "Generator" ) << "Total cross section: "
        << xsec*1.e3 << " +/- " << err*1.e3 << " fb.";
    else if ( xsec < 0.5e3 )
      CG_INFO( "Generator" ) << "Total cross section: "
        << xsec << " +/- " << err << " pb.";
    else if ( xsec < 0.5e6 )
      CG_INFO( "Generator" ) << "Total cross section: "
        << xsec*1.e-3 << " +/- " << err*1.e-3 << " nb.";
    else if ( xsec < 0.5e9 )
      CG_INFO( "Generator" ) << "Total cross section: "
        << xsec*1.e-6 << " +/- " << err*1.e-6 << " µb.";
    else
      CG_INFO( "Generator" ) << "Total cross section: "
        << xsec*1.e-9 << " +/- " << err*1.e-9 << " mb.";
  }

  void
  Generator::setIntegrator( std::unique_ptr<Integrator> integ )
  {
    if ( !integ ) {
      if ( !parameters_->integrator )
        throw CG_FATAL( "Generator:integrate" ) << "No integrator parameters found!";
      if ( parameters_->integrator->name<std::string>().empty() )
        parameters_->integrator->setName<std::string>( "Vegas" );
      integ = IntegratorFactory::get().build( *parameters_->integrator );
    }
    integrator_ = std::move( integ );
    integrator_->setIntegrand( generator_->integrand() );
    generator_->setIntegrator( integrator_.get() );
    CG_INFO( "Generator:integrator" )
      << "Generator will use a " << integrator_->name() << "-type integrator.";
  }

  void
  Generator::integrate()
  {
    clearRun();
    if ( !parameters_->hasProcess() )
      throw CG_FATAL( "Generator:integrate" )
        << "Trying to integrate while no process is specified!";
    const size_t ndim = parameters_->process().ndim();
    if ( ndim < 1 )
      throw CG_FATAL( "Generator:computePoint" )
        << "Invalid phase space dimension (ndim=" << ndim << ")!";

    // first destroy and recreate the integrator instance
    setIntegrator( nullptr );

    CG_DEBUG( "Generator:integrate" )
      << "New integrator instance created for " << ndim << "-dimensional integration.";

    integrator_->integrate( result_, result_error_ );

    for ( auto& mod : parameters_->eventModifiersSequence() )
      mod->setCrossSection( result_, result_error_ );
    for ( auto& mod : parameters_->outputModulesSequence() )
      mod->setCrossSection( result_, result_error_ );
  }

  const Event&
  Generator::generateOneEvent( Event::callback callback )
  {
    generate( 1, callback );
    return parameters_->process().event();
  }

  void
  Generator::generate( size_t num_events, Event::callback callback )
  {
    if ( !parameters_ )
      throw CG_FATAL( "Generator:generate" ) << "No steering parameters specified!";

    for ( auto& mod : parameters_->outputModulesSequence() )
      mod->initialise( *parameters_ );

    //--- if invalid argument, retrieve from runtime parameters
    if ( num_events < 1 )
      num_events = parameters_->generation().maxgen;

    CG_INFO( "Generator" )
      << utils::s( "event", num_events, true ) << " will be generated.";

    const utils::Timer tmr;

    //--- launch the event generation

    generator_->generate( num_events, callback );

    const double gen_time_s = tmr.elapsed();
    const double rate_ms = ( parameters_->numGeneratedEvents() > 0 )
      ? gen_time_s/parameters_->numGeneratedEvents()*1.e3 : 0.;
    const double equiv_lumi = parameters_->numGeneratedEvents()/crossSection();
    CG_INFO( "Generator" )
      << utils::s( "event", parameters_->numGeneratedEvents() )
      << " generated in " << gen_time_s << " s "
      << "(" << rate_ms << " ms/event).\n\t"
      << "Equivalent luminosity: " << utils::format( "%g", equiv_lumi ) << " pb^-1.";
  }

  void
  Generator::dumpModules() const
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
        << utils::boldify( "alpha(s) evolution algorithms" );
      if ( AlphaSFactory::get().modules().empty() )
        oss << "\n>>> " << utils::colourise( "none found", utils::Colour::red ) << " <<<";
      for ( const auto& mod : AlphaSFactory::get().modules() )
        oss << "\n> " << utils::colourise( mod, utils::Colour::green );
    }
    CG_INFO( "Generator:dumpModules" ) << oss.str();
  }
}
