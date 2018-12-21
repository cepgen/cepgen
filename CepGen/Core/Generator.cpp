#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/Timer.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Event/Event.h"

#include <fstream>
#include <chrono>
#include <atomic>

namespace cepgen
{
  namespace utils { std::atomic<int> gSignal; }
  Generator::Generator() :
    parameters( std::unique_ptr<Parameters>( new Parameters ) ), result_( -1. ), result_error_( -1. )
  {
    CG_DEBUG( "Generator:init" ) << "Generator initialized";
    try {
      printHeader();
    } catch ( const Exception& e ) {
      e.dump();
    }
    // Random number initialization
    std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
    srandom( time.time_since_epoch().count() );
  }

  Generator::Generator( Parameters* ip ) :
    parameters( ip ), result_( -1. ), result_error_( -1. )
  {}

  Generator::~Generator()
  {
    if ( parameters->generation.enabled
      && parameters->process() && parameters->numGeneratedEvents() > 0 ) {
      CG_INFO( "Generator" )
        << "Mean generation time / event: "
        << parameters->totalGenerationTime()*1.e3/parameters->numGeneratedEvents()
        << " ms.";
    }
  }

  size_t
  Generator::numDimensions() const
  {
    if ( !parameters->process() )
      return 0;

    parameters->process()->addEventContent();
    parameters->process()->setKinematics( parameters->kinematics );
    return parameters->process()->numDimensions();
  }

  void
  Generator::clearRun()
  {
    integrator_.reset();
    parameters->process()->first_run = true;
    result_ = result_error_ = -1.;
    {
      std::ostringstream os;
      for ( const auto& pr : cepgen::proc::ProcessesHandler::get().modules() )
        os << " " << pr;
      CG_DEBUG( "Generator:clearRun" ) << "Processes handled:" << os.str() << ".";
    }
  }

  void
  Generator::setParameters( const Parameters& ip )
  {
    parameters.reset( new Parameters( (Parameters&)ip ) ); // copy constructor
  }

  void
  Generator::printHeader()
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
    CG_INFO( "Generator" ) << os.str();
  }

  double
  Generator::computePoint( double* x )
  {
    double res = integrand::eval( x, numDimensions(), (void*)parameters.get() );
    std::ostringstream os;
    for ( unsigned int i = 0; i < numDimensions(); ++i )
      os << x[i] << " ";
    CG_DEBUG( "Generator:computePoint" )
      << "Result for x[" << numDimensions() << "] = { " << os.str() << "}:\n\t"
      << res << ".";
    return res;
  }

  void
  Generator::computeXsection( double& xsec, double& err )
  {
    CG_INFO( "Generator" ) << "Starting the computation of the process cross-section.";

    clearRun();
    integrate();

    xsec = result_;
    err = result_error_;

    if ( xsec < 1.e-2 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e3 << " +/- " << err*1.e3 << " fb.";
    else if ( xsec > 5.e2 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e-3 << " +/- " << err*1.e-3 << " nb.";
    else
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec << " +/- " << err << " pb.";
  }

  void
  Generator::integrate()
  {
    // first destroy and recreate the integrator instance
    if ( !integrator_ )
      integrator_ = std::unique_ptr<Integrator>( new Integrator( numDimensions(), integrand::eval, parameters.get() ) );
    else if ( integrator_->dimensions() != numDimensions() )
      integrator_.reset( new Integrator( numDimensions(), integrand::eval, parameters.get() ) );

    CG_DEBUG( "Generator:newInstance" )
      << "New integrator instance created\n\t"
      << "Considered topology: " << parameters->kinematics.mode << " case\n\t"
      << "Will proceed with " << numDimensions() << "-dimensional integration.";

    const int res = integrator_->integrate( result_, result_error_ );
    if ( res != 0 )
      throw CG_FATAL( "Generator" )
        << "Error while computing the cross-section!\n\t"
        << "GSL error: " << gsl_strerror( res ) << ".";
  }

  std::shared_ptr<Event>
  Generator::generateOneEvent()
  {
    integrator_->generateOne();

    parameters->addGenerationTime( parameters->process()->last_event->time_total );
    return parameters->process()->last_event;
  }

  void
  Generator::generate( std::function<void( const Event&, unsigned long )> callback )
  {
    const utils::Timer tmr;

    CG_INFO( "Generator" )
      << parameters->generation.maxgen << " events will be generated.";

    integrator_->generate( parameters->generation.maxgen, callback );

    const double gen_time_s = tmr.elapsed();
    CG_INFO( "Generator" )
      << parameters->generation.ngen << " events generated "
      << "in " << gen_time_s << " s "
      << "(" << gen_time_s/parameters->generation.ngen*1.e3 << " ms/event).";
  }
}
