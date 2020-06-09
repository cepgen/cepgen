#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"

#include "CepGen/Processes/Process.h"

#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Core/Exception.h"

#include <iostream>

using namespace std;

/// Generic process to test the integrator instance
class TestProcess : public cepgen::proc::Process
{
  public:
    /// Test process constructor
    inline explicit TestProcess( const string& func_mod, const string& formula, const vector<string>& args ) :
      cepgen::proc::Process( cepgen::ParametersList().set<string>( "description", formula ), false ),
      variables_( args.size() ) {
      funct_ = cepgen::utils::FunctionalFactory::get().build( func_mod,
        cepgen::ParametersList()
          .set<string>( "expression", formula )
          .set<vector<string> >( "arguments", args )
      );
    }
    /// Process cloning method
    cepgen::proc::ProcessPtr clone( const cepgen::ParametersList& ) const override {
      return cepgen::proc::ProcessPtr( new TestProcess( *this ) );
    }
    /// Dummy function to be called on phase space definition
    void prepareKinematics() override {
      for ( auto& var : variables_ )
        defineVariable( var, cepgen::proc::Process::Mapping::linear );
    }
    /// Dummy function to be called on events generation
    void fillKinematics( bool ) override {}
    /// Generic formula to compute a weight out of a point in the phase space
    double computeWeight() override {
      return funct_->operator()( variables_ );
    }

  private:
    vector<double> variables_;
    std::shared_ptr<cepgen::utils::Functional> funct_;
};

int
main( int argc, char* argv[] )
{
  bool debug;
  double num_sigma;
  string integrator, func_mod;

  cepgen::initialise();

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "num-sigma,n", "max. number of std.dev.", &num_sigma, 5. )
    .addOptionalArgument( "debug,d", "debugging mode", &debug, false )
    .addOptionalArgument( "integrator,i", "type of integrator used", &integrator, "Vegas" )
    .addOptionalArgument( "functional,f", "type of functional parser user", &func_mod, "ROOT" )
    .parse();

  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;
  else
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  cepgen::Parameters params;
  //--- integrator definition
  auto integr = cepgen::IntegratorFactory::get().build( integrator );

  //--- tests definition
  struct test_t
  {
    cepgen::proc::Process* process;
    double result;
  };
  vector<test_t> tests = {
    { new TestProcess( func_mod, "x^2+y^2", { "x", "y" } ), 2./3 },
    { new TestProcess( func_mod, "x+y^2+z^3", { "x", "y", "z" } ), 13./12. },
    { new TestProcess( func_mod, "1./(1.-cos(x*3.141592654)*cos(y*3.141592654)*cos(z*3.141592654))", { "x", "y", "z" } ), 1.3932039296856768591842462603255 },
  };

  //--- integration part
  size_t i = 0;
  double result, error;
  for ( const auto& test : tests ) {
    params.setProcess( test.process );
    cepgen::Integrand integrand( &params );
    integr->setIntegrand( integrand );
    integr->integrate( result, error );
    if ( fabs( test.result - result ) > num_sigma * error )
      throw CG_FATAL( "main" ) << "Test " << i << ": pull = " << fabs( test.result-result )/error << ".";
    cout << "Test " << i << " passed!" << endl;
    if ( debug )
      CG_INFO( "main" ) << "Test " << i << ": ref.: " << test.result << ", result: " << result << " +/- " << error << ".";
    ++i;
  }

  return 0;
}

