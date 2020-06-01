#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Processes/Process.h"

#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/GridParameters.h"

#include "CepGen/Modules/IntegratorFactory.h"

#include <chrono>

namespace cepgen
{
  Generator::Generator() :
    parameters_( new Parameters ),
    result_( -1. ), result_error_( -1. )
  {
    CG_DEBUG( "Generator:init" ) << "Generator initialized";
    static bool init = false;
    if ( !init ) {
      initialise();
      init = true;
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
  {
    if ( parameters_->timeKeeper() )
      CG_INFO( "Generator:destructor" )
        << parameters_->timeKeeper()->summary();
  }

  void
  Generator::clearRun()
  {
    generator_.reset( new GeneratorWorker( parameters_.get() ) );
    result_ = result_error_ = -1.;
    parameters_->prepareRun();
  }

  Parameters&
  Generator::parameters()
  {
    return *parameters_;
  }

  void
  Generator::setParameters( Parameters* ip )
  {
    parameters_.reset( ip );
    if ( parameters_->hasProcess() )
      parameters_->process().setKinematics( parameters_->kinematics );
  }

  double
  Generator::computePoint( const std::vector<double>& coord )
  {
    if ( !generator_ )
      clearRun();
    if ( !parameters_->hasProcess() )
      throw CG_FATAL( "Generator:computePoint" )
        << "Trying to compute a point with no process specified!";
    const size_t ndim = parameters_->process().ndim();
    if ( coord.size() != ndim )
      throw CG_FATAL( "Generator:computePoint" )
        << "Invalid phase space dimension (ndim=" << ndim << ")!";
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
    CG_TICKER( parameters_->timeKeeper() );

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
    CG_TICKER( parameters_->timeKeeper() );

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
    CG_TICKER( parameters_->timeKeeper() );

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
}
