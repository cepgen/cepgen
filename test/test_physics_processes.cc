#include "CepGen/Cards/PythonHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Timer.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/ProgressBar.h"

#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace cepgen;

int main( int argc, char* argv[] )
{
  double num_sigma;
  string cfg_filename;
  string integrator;
  bool debug, quiet;

  ArgumentsParser( argc, argv )
    .addArgument( "cfg", "configuration file", &cfg_filename, 'c' )
    .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
    .addOptionalArgument( "quiet", "quiet mode", false, &quiet, 'q' )
    .addOptionalArgument( "num-sigma", "max. number of std.dev.", 3., &num_sigma, 'n' )
    .addOptionalArgument( "integrator", "type of integrator used", "vegas", &integrator, 'i' )
    .parse();

  if ( debug )
    utils::Logger::get().level = utils::Logger::Level::information;
  else if ( quiet )
    utils::Logger::get().level = utils::Logger::Level::error;

  utils::Timer tmr;
  Generator gen;
  IntegratorType integr;

  if ( integrator == "plain" )
    integr = IntegratorType::plain;
  else if ( integrator == "vegas" )
    integr = IntegratorType::Vegas;
  else if ( integrator == "miser" )
    integr = IntegratorType::MISER;
  else
    throw CG_FATAL( "main" ) << "Unhandled integrator type: " << integrator << ".";

  CG_LOG( "main" ) << "Testing with " << integr << " integrator.";

  vector<string> failed_tests, passed_tests;

  CG_INFO( "main" ) << "Initial configuration time: " << tmr.elapsed()*1.e3 << " ms.";
  tmr.reset();

  utils::AbortHandler ctrl_c;

  struct Test
  {
    string filename;
    double ref_cs, err_ref_cs;
  };
  vector<Test> tests;

  { // parse the various tests to be performed
    ifstream cfg( cfg_filename );
    string line;
    while ( !cfg.eof() ) {
      getline( cfg, line );
      if ( line[0] == '#' || line.empty() )
        continue;
      stringstream os( line );
      Test test;
      os >> test.filename >> test.ref_cs >> test.err_ref_cs;
      tests.emplace_back( test );
    }
  }

  CG_INFO( "main" )
    << "Will run " << utils::s( "test", tests.size() ) << ".";

  std::unique_ptr<utils::ProgressBar> progress;
  if ( debug )
    progress.reset( new utils::ProgressBar( tests.size() ) );

  try {
    unsigned short num_tests = 0;
    for ( const auto& test : tests ) {
      gen.parameters().clearProcess();

      const std::string filename = "test_processes/"+test.filename+"_cfg.py";
      gen.setParameters( cepgen::card::PythonHandler( filename ).parameters() );
      gen.parameters().integrator->setName<int>( (int)integr );
      CG_INFO( "main" )
        << "Process: "<< gen.parameters().processName() << "\n\t"
        << "File: " << filename << "\n\t"
        << "Configuration time: " << tmr.elapsed()*1.e3 << " ms.";

      tmr.reset();

      double new_cs, err_new_cs;
      gen.computeXsection( new_cs, err_new_cs );

      const double ratio = new_cs/test.ref_cs;
      const double err_ratio = ratio * hypot( err_new_cs/new_cs, test.err_ref_cs/test.ref_cs );
      const double pull = ( new_cs-test.ref_cs ) / hypot( err_new_cs, test.err_ref_cs );

      const bool success = fabs( pull ) < num_sigma;

      CG_INFO( "main" )
        << "Computed cross section:\n\t"
        << "Ref.   = " << test.ref_cs << " +/- " << test.err_ref_cs << "\n\t"
        << "CepGen = " << new_cs << " +/- " << err_new_cs << "\n\t"
        << "Ratio: " << ratio << " +/- " << err_ratio << "\n\t"
        << "Pull: " << pull << " (abs(pull) "
        << ( success ? "<" : ">" ) << " " << num_sigma << ").";

      CG_INFO( "main" )
        << "Computation time: " << tmr.elapsed()*1.e3 << " ms.";
      tmr.reset();

      const string test_res = utils::format(
        "%-26s\tref=%g\tgot=%g\tpull=%+g",
        test.filename.c_str(), test.ref_cs, new_cs, pull );
      if ( success )
        passed_tests.emplace_back( test_res );
      else
        failed_tests.emplace_back( test_res );
      num_tests++;
      if ( debug )
        progress->update( num_tests );
      CG_LOG( "main" )
        << "Test " << num_tests << "/" << tests.size() << " finished. "
        << "Success: " << utils::yesno( success ) << ".";
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
