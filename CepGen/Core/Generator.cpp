#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/Timer.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/Process.h"

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
    parameters_( ip ), result_( -1. ), result_error_( -1. )
  {}

  Generator::~Generator()
  {}

  void
  Generator::clearRun( bool clear_proc )
  {
    if ( parameters_->process() ) {
      if ( clear_proc )
        parameters_->clearProcess();
      else {
        parameters_->process()->first_run = true;
        parameters_->process()->addEventContent();
        parameters_->process()->setKinematics( parameters_->kinematics );
      }
    }
    result_ = result_error_ = -1.;
  }

  Parameters&
  Generator::parameters()
  {
    return *parameters_.get();
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

    if ( !parameters_->process() )
      throw CG_FATAL( "Generator:computePoint" )
        << "Trying to compute a point with no process specified!";
    const size_t ndim = parameters_->process()->ndim();
    double res = integrand::eval( x, ndim, (void*)parameters_.get() );
    std::ostringstream os;
    std::string sep;
    for ( size_t i = 0; i < ndim; ++i )
      os << sep << x[i], sep = ", ";
    CG_DEBUG( "Generator:computePoint" )
      << "Result for x[" << ndim << "] = {" << os.str() << "}:\n\t"
      << res << ".";
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
  Generator::integrate()
  {
    clearRun();
    result_ = result_error_ = 0.;

    // first destroy and recreate the integrator instance
    if ( !parameters_->process() )
      throw CG_FATAL( "Generator:integrate" )
        << "Trying to integrate while no process is specified!";

    const size_t ndim = parameters_->process()->ndim();
    if ( !integrator_ || integrator_->dimensions() != ndim )
      integrator_.reset( new Integrator( ndim, integrand::eval, *parameters_ ) );

    CG_DEBUG( "Generator:integrate" )
      << "New integrator instance created\n\t"
      << "Considered topology: " << parameters_->kinematics.mode << " case\n\t"
      << "Will proceed with " << ndim << "-dimensional integration.";

    integrator_->integrate( result_, result_error_ );
  }

  const Event&
  Generator::generateOneEvent()
  {
    integrator_->generateOne();
    return parameters_->process()->event();
  }

  void
  Generator::generate( std::function<void( const Event&, unsigned long )> callback )
  {
    const utils::Timer tmr;

    CG_INFO( "Generator" )
      << parameters_->generation().maxgen << " events will be generated.";

    integrator_->generate( parameters_->generation().maxgen, callback );

    const double gen_time_s = tmr.elapsed();
    const double rate_ms = ( parameters_->numGeneratedEvents() > 0 )
      ? gen_time_s/parameters_->numGeneratedEvents()*1.e3 : 0.;
    CG_INFO( "Generator" )
      << utils::s( "event", parameters_->numGeneratedEvents() )
      << " generated in " << gen_time_s << " s "
      << "(" << rate_ms << " ms/event).";
  }
}
