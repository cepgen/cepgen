#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include <fstream>

namespace CepGen
{
  Generator::Generator() :
    parameters( std::unique_ptr<Parameters>( new Parameters ) ),
    cross_section_( -1. ), cross_section_error_( -1. ), has_cross_section_( false )
  {
    if ( Logger::get().level > Logger::Nothing ) {
      Debugging( "Generator initialized" );
      try { printHeader(); } catch ( Exception& e ) { e.dump(); }
    }
    srand( time( 0 ) ); // Random number initialization
  }

  Generator::Generator( Parameters* ip ) :
    parameters( ip )
  {}

  Generator::~Generator()
  {
    if ( parameters->generation.enabled && parameters->process() && parameters->process()->numGeneratedEvents()>0 ) {
      Information( Form( "Mean generation time / event: %.3f ms", parameters->process()->totalGenerationTime()*1.e3/parameters->process()->numGeneratedEvents() ) );
    }
  }

  size_t
  Generator::numDimensions() const
  {
    if ( !parameters->process() ) return 0;
    return parameters->process()->numDimensions( parameters->kinematics.mode );
  }

  void
  Generator::clearRun()
  {
    integrator_.reset();
    parameters->integrator.first_run = true;
    has_cross_section_ = false; // force the recreation of the integrator instance
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
    std::ostringstream os; os << std::endl;
    std::ifstream hf( "README" );
    if ( !hf.good() ) throw Exception( __PRETTY_FUNCTION__, "Failed to open README file", JustWarning );
    while ( true ) {
      if ( !hf.good() ) break;
      getline( hf, tmp );
      os << "\n " << tmp;
    }
    hf.close();
    Information( os.str().c_str() );
  }

  double
  Generator::computePoint( double* x )
  {
    prepareFunction();
    double res = f( x, numDimensions(), (void*)parameters.get() );
    std::ostringstream os;
    for ( unsigned int i=0; i<numDimensions(); i++ ) { os << x[i] << " "; }
    Debugging( Form( "Result for x[%zu] = ( %s):\n\t%10.6f", numDimensions(), os.str().c_str(), res ) );
    return res;
  }

  void
  Generator::computeXsection( double& xsec, double& err )
  {
    // first destroy and recreate the integrator instance
    if ( !integrator_ ) {
      integrator_ = std::unique_ptr<Integrator>( new Integrator( numDimensions(), f, parameters.get() ) );
    }
    else if ( integrator_->dimensions() != numDimensions() ) {
      integrator_.reset( new Integrator( numDimensions(), f, parameters.get() ) );
    }

    if ( Logger::get().level>=Logger::Debug ) {
      std::ostringstream topo; topo << parameters->kinematics.mode;
      Debugging( Form( "New integrator instance created\n\t"
                       "Considered topology: %s case\n\t"
                       "Will proceed with %d-dimensional integration", topo.str().c_str(), numDimensions() ) );
    }

    Information( "Starting the computation of the process cross-section" );

    try { prepareFunction(); } catch ( Exception& e ) { e.dump(); }

    has_cross_section_ = ( integrator_->integrate( cross_section_, cross_section_error_ ) == 0 );

    xsec = cross_section_;
    err = cross_section_error_;

    Information( Form( "Total cross section: %f +/- %f pb", xsec, err ) );
  }

  Event*
  Generator::generateOneEvent()
  {
    bool good = false;
    if ( !has_cross_section_ ) {
      computeXsection( cross_section_, cross_section_error_ );
    }
    while ( !good ) { good = integrator_->generateOneEvent(); }

    last_event = this->parameters->generation.last_event;
    return last_event.get();
  }

  void
  Generator::prepareFunction()
  {
    if ( !parameters->process() ) {
      throw Exception( __PRETTY_FUNCTION__, "No process defined!", FatalError );
    }
    Kinematics kin = parameters->kinematics;
    parameters->process()->addEventContent();
    parameters->process()->setKinematics( kin );
    Debugging( "Function prepared to be integrated!" );
  }
}

