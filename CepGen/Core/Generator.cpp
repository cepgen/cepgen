#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Event/Event.h"

#include <fstream>

namespace CepGen
{
  volatile int gSignal;
  Generator::Generator() :
    parameters( std::unique_ptr<Parameters>( new Parameters ) ),
    cross_section_( -1. ), cross_section_error_( -1. )
  {
    if ( Logger::get().level > Logger::Nothing ) {
      Debugging( "Generator" ) << "Generator initialized";
      try {
        printHeader();
      } catch ( Exception& e ) {
        e.dump();
      }
    }
    // Random number initialization
    struct timespec ts;
    if ( timespec_get( &ts, TIME_UTC ) != 0 )
      srandom( ts.tv_nsec ^ ts.tv_sec );
  }

  Generator::Generator( Parameters* ip ) :
    parameters( ip ),
    cross_section_( -1. ), cross_section_error_( -1. )
  {}

  Generator::~Generator()
  {
    if ( parameters->generation.enabled
      && parameters->process() && parameters->numGeneratedEvents() > 0 ) {
      Information( "Generator" )
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

    return parameters->process()->numDimensions( parameters->kinematics.mode );
  }

  void
  Generator::clearRun()
  {
    integrator_.reset();
    parameters->process()->first_run = true;
    cross_section_ = cross_section_error_ = -1.;
  }

  void
  Generator::setParameters( Parameters& ip )
  {
    parameters = std::unique_ptr<Parameters>( new Parameters( ip ) ); // copy constructor
  }

  void
  Generator::printHeader()
  {
    std::string tmp;
    std::ostringstream os; os << "version " << version() << std::endl;
    std::ifstream hf( "README" );
    if ( !hf.good() )
      throw InWarning( "Generator" ) << "Failed to open README file.";
    while ( true ) {
      if ( !hf.good() ) break;
      getline( hf, tmp );
      os << "\n " << tmp;
    }
    hf.close();
    Information( "Generator" ) << os.str();
  }

  double
  Generator::computePoint( double* x )
  {
    prepareFunction();

    double res = Integrand::eval( x, numDimensions(), (void*)parameters.get() );
    std::ostringstream os;
    for ( unsigned int i = 0; i < numDimensions(); ++i )
      os << x[i] << " ";
    Debugging( "Generator" )
      << "Result for x[" << numDimensions() << "] = ( " << os.str() << "):\n\t"
      << res << ".";
    return res;
  }

  void
  Generator::computeXsection( double& xsec, double& err )
  {
    Information( "Generator" ) << "Starting the computation of the process cross-section.";

    try {
      prepareFunction();
    } catch ( Exception& e ) {
      e.dump();
    }

    // first destroy and recreate the integrator instance
    if ( !integrator_ )
      integrator_ = std::unique_ptr<Integrator>( new Integrator( numDimensions(), Integrand::eval, parameters.get() ) );
    else if ( integrator_->dimensions() != numDimensions() )
      integrator_.reset( new Integrator( numDimensions(), Integrand::eval, parameters.get() ) );

    if ( Logger::get().level >= Logger::Debug )
      Debugging( "Generator" )
        << "New integrator instance created\n\t"
        << "Considered topology: " << parameters->kinematics.mode << " case\n\t"
        << "Will proceed with " << numDimensions() << "-dimensional integration.";

    const int res = integrator_->integrate( cross_section_, cross_section_error_ );
    if ( res != 0 )
      throw FatalError( "Generator" ) << "Error while computing the cross-section: return value = " << res << ".";

    xsec = cross_section_;
    err = cross_section_error_;

    if ( xsec < 1.e-2 )
      Information( "Generator" ) << "Total cross section: " << xsec*1.e3 << " +/- " << err*1.e3 << " fb.";
    else if ( xsec > 5.e2 )
      Information( "Generator" ) << "Total cross section: " << xsec*1.e-3 << " +/- " << err*1.e-3 << " nb.";
    else
      Information( "Generator" ) << "Total cross section: " << xsec << " +/- " << err << " pb.";
  }

  std::shared_ptr<Event>
  Generator::generateOneEvent()
  {
    if ( cross_section_ < 0. )
      computeXsection( cross_section_, cross_section_error_ );

    integrator_->generate( 1 );

    parameters->addGenerationTime( parameters->process()->last_event->time_total );
    return parameters->process()->last_event;
  }

  void
  Generator::generate( std::function<void( const Event&, unsigned long )> callback )
  {
    if ( cross_section_ < 0. )
      computeXsection( cross_section_, cross_section_error_ );

    Information( "Generator" )
      << parameters->generation.maxgen << " events will be generated.";

    integrator_->generate( parameters->generation.maxgen, callback );

    Information( "Generator" )
      << parameters->generation.ngen << " events generated.";
  }

  void
  Generator::prepareFunction()
  {
    if ( !parameters->process() )
      throw FatalError( "Generator" ) << "No process defined!";

    Kinematics kin = parameters->kinematics;
    parameters->process()->addEventContent();
    parameters->process()->setKinematics( kin );
    Debugging( "Generator" ) << "Function prepared to be integrated!";
  }
}

