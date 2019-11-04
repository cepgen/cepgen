#include "CepGen/Generator.h"
#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Timer.h"

#include "CepGen/Utils/ArgumentsParser.h"

#include "AbortHandler.h"

#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace cepgen;

int
main( int argc, char* argv[] )
{
  double num_sigma;
  string cfg_filename;
  string integrator;

  ArgumentsParser( argc, argv )
    .addArgument( "cfg", "configuration file", &cfg_filename, 'c' )
    .addOptionalArgument( "num-sigma", "max. number of std.dev.", 3., &num_sigma, 'n' )
    .addOptionalArgument( "integrator", "type of integrator used", "vegas", &integrator, 'i' )
    .parse();

  utils::Logger::get().level = utils::Logger::Level::error;

  utils::Timer tmr;
  Generator mg;
  auto& params = mg.parameters();

  params.integration().type = IntegratorType::Vegas;
  if ( integrator == "plain" )
    params.integration().type = IntegratorType::plain;
  else if ( integrator == "vegas" )
    params.integration().type = IntegratorType::Vegas;
  else if ( integrator == "miser" )
    params.integration().type = IntegratorType::MISER;
  else
    throw CG_FATAL( "main" ) << "Unhandled integrator type: " << integrator << ".";

  CG_LOG( "main" ) << "Testing with " << params.integration().type << " integrator.";

  unsigned short num_tests = 0;
  vector<string> failed_tests, passed_tests;

  CG_INFO( "main" ) << "Initial configuration time: " << tmr.elapsed()*1.e3 << " ms.";
  tmr.reset();

  utils::AbortHandler ctrl_c;

  ifstream cfg( cfg_filename );
  string line;
  try {
    while ( !cfg.eof() ) {
      getline( cfg, line );
      if ( line[0] == '#' || line.empty() )
        continue;
      string config;
      double ref_cs, err_ref_cs;
      stringstream os( line );
      os >> config >> ref_cs >> err_ref_cs;

      mg.setParameters( cepgen::card::Handler::parse( ( "test_processes/"+config+"_cfg.py" ).c_str() ) );
      cout << &params << endl;
      CG_INFO( "main" )
        << "Process: "<< params.processName() << "\n\t"
        << "Configuration time: " << tmr.elapsed()*1.e3 << " ms.";

      tmr.reset();
      mg.clearRun();
      double new_cs, err_new_cs;
      mg.computeXsection( new_cs, err_new_cs );

      const double sigma = ( new_cs-ref_cs ) / hypot( err_new_cs, err_ref_cs );

      CG_INFO( "main" )
        << "Computed cross section:\n\t"
        << "Ref.   = " << ref_cs << " +/- " << err_ref_cs << "\n\t"
        << "CepGen = " << new_cs << " +/- " << err_new_cs << "\n\t"
        << "Pull: " << sigma << ".";

      CG_INFO( "main" ) << "Computation time: " << tmr.elapsed()*1.e3 << " ms.";
      tmr.reset();

      const string test_res = Form( "%-26s\tref=%g\tgot=%g\tpull=%+g", config.c_str(), ref_cs, new_cs, sigma );
      if ( fabs( sigma ) < num_sigma )
        passed_tests.emplace_back( test_res );
      else
        failed_tests.emplace_back( test_res );
      num_tests++;
      cout << "Test " << passed_tests.size() << "/" << num_tests << " finished." << endl;
    }
  } catch ( const Exception& e ) {}
  if ( failed_tests.size() != 0 ) {
    ostringstream os_failed, os_passed;
    for ( const auto& fail : failed_tests )
      os_failed << " (*) " << fail << endl;
    for ( const auto& pass : passed_tests )
      os_passed << " (*) " << pass << endl;
    throw CG_FATAL( "main" )
      << "Some tests failed!\n"
      << os_failed.str() << "\n"
      << "Passed tests:\n"
      << os_passed.str() << ".";
  }

  CG_INFO( "main" ) << "ALL TESTS PASSED!";

  return 0;
}
