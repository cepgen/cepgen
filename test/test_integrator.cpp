#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Functional.h"

#include "CepGen/Processes/GenericProcess.h"

#include <iostream>

using namespace std;

/// Generic process to test the integrator instance
template<size_t N=3> class TestProcess : public cepgen::proc::GenericProcess
{
  public:
    TestProcess( const cepgen::ParametersList& params = cepgen::ParametersList() ) :
      cepgen::proc::GenericProcess( params, "test", ".oO TEST PROCESS Oo.", false ),
      funct_( "1./(1.-cos(x*_pi)*cos(y*_pi)*cos(z*_pi))", { "x", "y", "z" } ) {}
    TestProcess( const char* formula, const std::vector<std::string>& args ) :
      cepgen::proc::GenericProcess( cepgen::ParametersList(), "test", cepgen::Form( ".oO TEST PROCESS (%s) Oo.", formula ), false ),
      funct_( formula, args ) {}

    cepgen::proc::ProcessPtr clone( const cepgen::ParametersList& params ) const override {
      return cepgen::proc::ProcessPtr( new TestProcess<N>( *this ) );
    }

    void addEventContent() override {}
    /// Number of dimensions on which to perform the integration
    unsigned int numDimensions() const override { return N; }
    /// Generic formula to compute a weight out of a point in the phase space
    double computeWeight() override {
      std::array<double,N> args;
      std::copy_n( x_.begin(), N, args.begin() );
      return funct_.eval( args );
    }
    /// Dummy function to be called on events generation
    void fillKinematics( bool ) override {}

  private:
    cepgen::utils::Functional<N> funct_;
};

int
main( int argc, char* argv[] )
{
  const bool debug = ( argc > 2 && string( argv[2] ) == "debug" );
  if ( !debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;
  const double max_sigma = 3.0;

  cepgen::Generator mg;
  auto& params = mg.parameters();
  auto& integ = params.integration();
  if ( argc > 1 ) {
    const std::string int_type( argv[1] );
    if ( int_type == "plain" )
      integ.type = cepgen::IntegratorType::plain;
    if ( int_type == "vegas" )
      integ.type = cepgen::IntegratorType::Vegas;
    if ( int_type == "miser" )
      integ.type = cepgen::IntegratorType::MISER;
  }

  { // test 1
    const char* test = "Test 1";
    const double exact = 1.3932039296856768591842462603255;
    params.setProcess( new TestProcess<3> );
    mg.integrate();
    const double result = mg.crossSection(), error = mg.crossSectionError();
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << test << ": pull = " << fabs( exact-result )/error << ".";
    cout << test << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << test << ": ref.: " << exact << ", result: " << result << " +/- " << error << ".";
  }
  { // test 2
    const char* test = "Test 2";
    const double exact = 2./3.;
    params.setProcess( new TestProcess<2>( "x^2+y^2", { "x", "y" } ) );
    mg.integrate();
    const double result = mg.crossSection(), error = mg.crossSectionError();
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << test << ": pull = " << fabs( exact-result )/error << ".";
    cout << test << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << test << ": ref.: " << exact << ", result: " << result << " +/- " << error << ".";
  }
  { // test 3
    const char* test = "Test 3";
    const double exact = 13./12.;
    params.setProcess( new TestProcess<3>( "x+y^2+z^3", { "x", "y", "z" } ) );
    mg.integrate();
    const double result = mg.crossSection(), error = mg.crossSectionError();
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << test << ": pull = " << fabs( exact-result )/error << ".";
    cout << test << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << test << ": ref.: " << exact << ", result: " << result << " +/- " << error << ".";
  }

  return 0;
}

