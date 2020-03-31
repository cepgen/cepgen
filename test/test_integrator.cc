#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Integration/Integrator.h"

#include "CepGen/Processes/Process.h"

#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include <iostream>

using namespace std;

/// Generic process to test the integrator instance
template<size_t N> class TestProcess : public cepgen::proc::Process
{
  public:
    /// Test process constructor
    inline explicit TestProcess( const string& formula, const vector<string>& args ) :
      cepgen::proc::Process( cepgen::ParametersList().set<string>( "description", formula ), false ),
      funct_( formula, args ) {}
    /// Dummy function to be called on phase space definition
    void prepareKinematics() override {
      for ( auto& var : variables_ )
        defineVariable( var, cepgen::proc::Process::Mapping::linear );
    }
    /// Dummy function to be called on events generation
    void fillKinematics( bool ) override {}
    /// Generic formula to compute a weight out of a point in the phase space
    double computeWeight() override {
      array<double,N> args;
      copy_n( variables_.begin(), N, args.begin() );
      return funct_.eval( args );
    }

  private:
    array<double,N> variables_;
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
    .addOptionalArgument( "integrator", "type of integrator used", "Vegas", &integrator, 'i' )
    .parse();

  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;
  else
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  //--- integrator definition
  cepgen::Generator gen;
  gen.setIntegrator( cepgen::IntegratorFactory::get().build( integrator ) );

  //--- tests definition
  struct test_t
  {
    cepgen::proc::Process* process;
    double result;
  };
  vector<test_t> tests = {
    { new TestProcess<2>( "x^2+y^2", { "x", "y" } ), 2./3 },
    { new TestProcess<3>( "x+y^2+z^3", { "x", "y", "z" } ), 13./12. },
    { new TestProcess<3>( "1./(1.-cos(x*3.141592654)*cos(y*3.141592654)*cos(z*3.141592654))", { "x", "y", "z" } ), 1.3932039296856768591842462603255 },
  };

  //--- integration part
  size_t i = 0;
  for ( const auto& test : tests ) {
    gen.parameters().setProcess( test.process );
    gen.integrate();
    const double result = gen.crossSection(), error = gen.crossSectionError();
    if ( fabs( test.result - result ) > num_sigma * error )
      throw CG_FATAL( "main" ) << "Test " << i << ": pull = " << fabs( test.result-result )/error << ".";
    cout << "Test " << i << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << "Test " << i << ": ref.: " << test.result << ", result: " << result << " +/- " << error << ".";
    ++i;
  }

  return 0;
}

