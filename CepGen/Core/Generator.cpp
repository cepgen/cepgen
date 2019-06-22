#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/Timer.h"
#include "CepGen/Core/utils.h"

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
    parameters_( new Parameters ), result_( -1. ), result_error_( -1. )
  {
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
    parameters_( ip ), result_( -1. ), result_error_( -1. )
  {}

  Generator::~Generator()
  {
    if ( parameters_->generation().enabled
      && parameters_->process()
      && parameters_->numGeneratedEvents() > 0 )
      CG_INFO( "Generator" )
        << "Mean generation time / event: "
        << parameters_->totalGenerationTime()*1.e3/parameters_->numGeneratedEvents()
        << " ms.";
  }

  size_t
  Generator::numDimensions() const
  {
    if ( !parameters_->process() )
     return 0;
    return parameters_->process()->numDimensions();
  }

  void
  Generator::clearRun()
  {
    if ( parameters_->process() ) {
      parameters_->process()->first_run = true;
      parameters_->process()->addEventContent();
      parameters_->process()->setKinematics( parameters_->kinematics );
    }
    result_ = result_error_ = -1.;
    {
      std::ostringstream os;
      for ( const auto& pr : cepgen::proc::ProcessesHandler::get().modules() )
        os << " " << pr;
      CG_DEBUG( "Generator:clearRun" ) << "Processes handled:" << os.str() << ".";
    }
  }

  Parameters&
  Generator::parameters()
  {
    return *parameters_;
  }

  void
  Generator::setParameters( Parameters& ip )
  {
    parameters_.reset( new Parameters( ip ) ); // copy constructor
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
    CG_LOG( "Generator" ) << os.str();
  }

  double
  Generator::computePoint( double* x )
  {
    clearRun();

    double res = integrand::eval( x, numDimensions(), (void*)parameters_.get() );
    std::ostringstream os;
    for ( size_t i = 0; i < numDimensions(); ++i )
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

    integrate();

    xsec = result_;
    err = result_error_;

    if ( xsec < 1.e-2 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e3 << " +/- " << err*1.e3 << " fb.";
    else if ( xsec < 0.5e3 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec << " +/- " << err << " pb.";
    else if ( xsec < 0.5e6 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e-3 << " +/- " << err*1.e-3 << " nb.";
    else if ( xsec < 0.5e9 )
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e-6 << " +/- " << err*1.e-6 << " µb.";
    else
      CG_INFO( "Generator" )
        << "Total cross section: " << xsec*1.e-9 << " +/- " << err*1.e-9 << " mb.";
  }

  void
  Generator::integrate()
  {
    clearRun();

    // first destroy and recreate the integrator instance
    if ( !integrator_ || integrator_->dimensions() != numDimensions() )
      integrator_.reset( new Integrator( numDimensions(), integrand::eval, *parameters_ ) );

    CG_DEBUG( "Generator:newInstance" )
      << "New integrator instance created\n\t"
      << "Considered topology: " << parameters_->kinematics.mode << " case\n\t"
      << "Will proceed with " << numDimensions() << "-dimensional integration.";

    integrator_->integrate( result_, result_error_ );
  }

  std::shared_ptr<Event>
  Generator::generateOneEvent()
  {
    integrator_->generateOne();
    return parameters_->process()->last_event;
  }

  void
  Generator::generate( std::function<void( const Event&, unsigned long )> callback )
  {
    const utils::Timer tmr;

    CG_INFO( "Generator" )
      << parameters_->generation().maxgen << " events will be generated.";

    integrator_->generate( parameters_->generation().maxgen, callback );

    const double gen_time_s = tmr.elapsed();
    CG_INFO( "Generator" )
      << parameters_->numGeneratedEvents() << " event" << utils::s( parameters_->numGeneratedEvents() )
      << " generated "
      << "in " << gen_time_s << " s "
      << "(" << gen_time_s/parameters_->numGeneratedEvents()*1.e3 << " ms/event).";
  }
}
