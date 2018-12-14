#include "CepGen/Generator.h"
//#include "CepGen/Parameters.h"
#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Timer.h"

#include "AbortHandler.h"

#include <unordered_map>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <math.h>

using namespace std;
using namespace cepgen;

int
main( int argc, char* argv[] )
{
  utils::AbortHandler ctrl_c;

  const double num_sigma = 3.0;

  if ( argc < 2 )
    throw CG_FATAL( "mail" ) << "Usage: " << argv[0] << " [config file] [integrator=vegas]";

    utils::Logger::get().level = utils::Logger::Level::nothing;

  utils::Timer tmr;
  Generator mg;

  mg.parameters->integrator.type = IntegratorType::Vegas;
  if ( argc > 2 ) {
    std::string integrator( argv[2] );
    if ( integrator == "plain" )
      mg.parameters->integrator.type = IntegratorType::plain;
    else if ( integrator == "vegas" )
      mg.parameters->integrator.type = IntegratorType::Vegas;
    else if ( integrator == "miser" )
      mg.parameters->integrator.type = IntegratorType::MISER;
    else
      throw CG_FATAL( "main" ) << "Unhandled integrator type: " << integrator << ".";
  }

  { cout << "Testing with " << mg.parameters->integrator.type << " integrator" << endl; }

  unsigned short num_tests = 0;
  vector<string> failed_tests, passed_tests;

  CG_INFO( "main" ) << "Initial configuration time: " << tmr.elapsed()*1.e3 << " ms.";
  tmr.reset();

  ifstream cfg( argv[1] );
  string line;
  while ( !cfg.eof() ) {
    getline( cfg, line );
    if ( line[0] == '#' || line.empty() )
      continue;
    string config;
    double ref_cs, err_ref_cs;
    stringstream os( line );
    os >> config >> ref_cs >> err_ref_cs;

    try {
      mg.setParameters( cepgen::card::Handler::parse( ( "test_processes/"+config+"_cfg.py" ).c_str() ) );
      //CG_INFO( "main" ) << mg.parameters.get();
      CG_INFO( "main" )
        << "Process: "<< mg.parameters->processName() << "\n\t"
        << "Configuration time: " << tmr.elapsed()*1.e3 << " ms.";

      tmr.reset();
      mg.clearRun();
      double cg_cs, err_cg_cs;
      mg.computeXsection( cg_cs, err_cg_cs );

      const double sigma = fabs( ref_cs-cg_cs ) / std::hypot( err_cg_cs, err_ref_cs );

      CG_INFO( "main" )
        << "Computed cross section:\n\t"
        << "Ref.   = " << ref_cs << " +/- " << err_ref_cs << "\n\t"
        << "CepGen = " << cg_cs << " +/- " << err_cg_cs << "\n\t"
        << "Pull: " << sigma << ".";

      CG_INFO( "main" ) << "Computation time: " << tmr.elapsed()*1.e3 << " ms.";
      tmr.reset();

      const string test_res = Form( "%-16s\tref=%g\tgot=%g\tpull=%+g", config.c_str(), ref_cs, cg_cs, sigma );
      if ( fabs( sigma ) < num_sigma )
        passed_tests.emplace_back( test_res );
      else
        failed_tests.emplace_back( test_res );
      num_tests++;
      cout << "Test " << passed_tests.size() << "/"
                      << num_tests << " passed!" << endl;
    } catch ( Exception& e ) {}
  }
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

