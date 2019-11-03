#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Functional.h"

#include "CepGen/Modules/Process.h"

#include "ArgumentsParser.h"

#include <iostream>

using namespace std;

/// Generic process to test the integrator instance
template<size_t N=3> class TestProcess : public cepgen::proc::Process
{
  public:
    TestProcess( const cepgen::ParametersList& params = cepgen::ParametersList() ) :
      cepgen::proc::Process( params, "test", ".oO TEST PROCESS Oo.", false ),
      funct_( "1./(1.-cos(x*_pi)*cos(y*_pi)*cos(z*_pi))", { "x", "y", "z" } ) {}
    TestProcess( const char* formula, const std::vector<std::string>& args ) :
      cepgen::proc::Process( cepgen::ParametersList(), "test", cepgen::Form( ".oO TEST PROCESS (%s) Oo.", formula ), false ),
      funct_( formula, args ) {}

    cepgen::proc::ProcessPtr clone( const cepgen::ParametersList& params ) const override {
      return cepgen::proc::ProcessPtr( new TestProcess<N>( *this ) );
    }

    void prepareKinematics() override {
      for ( auto& var : variables_ )
        defineVariable( var, cepgen::proc::Process::Mapping::linear );
    }
    /// Generic formula to compute a weight out of a point in the phase space
    double computeWeight() override {
      std::array<double,N> args;
      std::copy_n( variables_.begin(), N, args.begin() );
      return funct_.eval( args );
    }
    /// Dummy function to be called on events generation
    void fillKinematics( bool ) override {}

  private:
    std::array<double,N> variables_;
    cepgen::utils::Functional<N> funct_;
};

int
main( int argc, char* argv[] )
{
  bool debug;
  double num_sigma;
  string integrator;

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "num-sigma", "max. number of std.dev.", 3., &num_sigma, 'n' )
    .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
    .addOptionalArgument( "integrator", "type of integrator used", "vegas", &integrator, 'i' )
    .parse();

  if ( !debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  cepgen::Generator mg;
  auto& params = mg.parameters();
  auto& integ = params.integration();
  if ( integrator == "plain" )
    integ.type = cepgen::IntegratorType::plain;
  else if ( integrator == "vegas" )
    integ.type = cepgen::IntegratorType::Vegas;
  else if ( integrator == "miser" )
    integ.type = cepgen::IntegratorType::MISER;

  { // test 1
    const char* test = "Test 1";
    const double exact = 1.3932039296856768591842462603255;
    params.setProcess( new TestProcess<3> );
    mg.integrate();
    const double result = mg.crossSection(), error = mg.crossSectionError();
    if ( fabs( exact - result ) > num_sigma * error )
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
    if ( fabs( exact - result ) > num_sigma * error )
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
    if ( fabs( exact - result ) > num_sigma * error )
      throw CG_FATAL( "main" ) << test << ": pull = " << fabs( exact-result )/error << ".";
    cout << test << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << test << ": ref.: " << exact << ", result: " << result << " +/- " << error << ".";
  }

  return 0;
}

