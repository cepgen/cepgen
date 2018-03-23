#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Core/Timer.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "abort.h"

#include <unordered_map>
#include <assert.h>
#include <string.h>

using namespace std;

int
main( int argc, char* argv[] )
{
  typedef vector<pair<const char*,pair<double,double> > > KinematicsMap;
  typedef vector<pair<float, KinematicsMap> > ValuesAtCutMap;

  AbortHandler ctrl_c;

  // values defined at pt(single lepton)>15 GeV, |eta(single lepton)|<2.5, mX<1000 GeV
  // process -> { pt cut -> { kinematics -> ( sigma, delta(sigma) ) } }
  vector<pair<const char*,ValuesAtCutMap> > values_map = {
    //--- LPAIR values at sqrt(s) = 13 TeV
    /*{ "lpair", {
      { 3.0, { // pt cut
        { "elastic",    { 2.0871703e1, 3.542e-2 } },
        { "singlediss", { 1.5042536e1, 3.256e-2 } },
        { "doublediss", { 1.38835e1, 4.03018e-2 } }
      } },
      { 15.0, { // pt cut
        { "elastic",    { 4.1994803e-1, 8.328e-4 } },
        { "singlediss", { 4.8504819e-1, 1.171e-3 } },
        { "doublediss", { 6.35650e-1, 1.93968e-3 } }
      } },
    } },*/
    //--- PPtoWW values
    { "pptoww", {
      { 0.0, { // pt cut
        { "elastic",         { 0.273, 0.01 } },
        { "elastic",         { 0.273, 0.01 } }, // FIXME
        { "singlediss_lux",  { 0.409, 0.01 } },
        { "doublediss_lux",  { 1.090, 0.01 } },
        { "singlediss_allm", { 0.318, 0.01 } },
        { "doublediss_allm", { 0.701, 0.01 } }
      } }
    } },
    //--- PPtoLL values
    { "pptoll", {
      { 3.0, { // pt cut
        { "elastic",       { 2.0936541e1, 1.4096e-2 } },
        { "singlediss_su", { 1.4844881e1, 2.0723e-2 } }, // SU, qt<50
        { "doublediss_su", { 1.3580637e1, 2.2497e-2 } }, // SU, qt<50
      } },
      { 15.0, { // pt cut
        { "elastic",       { 4.2515888e-1, 3.0351e-4 } },
        { "singlediss_su", { 4.4903253e-1, 5.8970e-4 } }, // SU, qt<50
        { "doublediss_su", { 5.1923819e-1, 9.6549e-4 } }, // SU, qt<50
        /*{ "2_singlediss", { 4.6710287e-1, 6.4842e-4 } }, // SU, qt<500
        { "3_doublediss", { 5.6316353e-1, 1.1829e-3 } }, // SU, qt<500*/
      } },
    } },
  };

  const double num_sigma = 3.0;

  if ( argc < 3 || strcmp( argv[2], "debug" ) != 0 ) {
    CepGen::Logger::get().level = CepGen::Logger::Nothing;
  }

  Timer tmr;
  CepGen::Generator mg;

  if ( argc > 1 && strcmp( argv[1], "plain" ) == 0 )
    mg.parameters->integrator.type = CepGen::Integrator::Plain;
  if ( argc > 1 && strcmp( argv[1], "vegas" ) == 0 )
    mg.parameters->integrator.type = CepGen::Integrator::Vegas;
  if ( argc > 1 && strcmp( argv[1], "miser" ) == 0 )
    mg.parameters->integrator.type = CepGen::Integrator::MISER;

  { cout << "Testing with " << mg.parameters->integrator.type << " integrator" << endl; }

  mg.parameters->kinematics.setSqrtS( 13.e3 );
  mg.parameters->kinematics.cuts.central[CepGen::Cuts::eta_single].in( -2.5, 2.5 );
  mg.parameters->kinematics.cuts.remnants[CepGen::Cuts::mass].max() = 1000.;
  //mg.parameters->integrator.ncvg = 50000;

  Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  unsigned short num_tests = 0, num_tests_passed = 0;
  vector<string> failed_tests, passed_tests;

  try {
    for ( const auto& values_vs_generator : values_map ) { // loop over all generators
      const string generator = values_vs_generator.first;
      if ( generator == "lpair"  )
        mg.parameters->setProcess( new CepGen::Process::GamGamLL );
      else if ( generator == "pptoll" ) {
        mg.parameters->setProcess( new CepGen::Process::PPtoLL );
        mg.parameters->kinematics.cuts.initial[CepGen::Cuts::qt] = { 0., 50. };
      }
      else if ( generator == "pptoww" ) {
        mg.parameters->setProcess( new CepGen::Process::PPtoWW );
        mg.parameters->kinematics.setSqrtS( 13.e3 );
        //mg.parameters->kinematics.cuts.initial[CepGen::Cuts::qt] = { 0., 50. };
      }
      else { InError( Form( "Unrecognized generator mode: %s", values_vs_generator.first ) ); break; }

      for ( const auto& values_vs_cut : values_vs_generator.second ) { // loop over the single lepton pT cut
        mg.parameters->kinematics.cuts.central[CepGen::Cuts::pt_single].min() = values_vs_cut.first;
        for ( const auto& values_vs_kin : values_vs_cut.second ) { // loop over all possible kinematics
          const string kin_mode = values_vs_kin.first;

          if ( kin_mode.find( "elastic"    ) != string::npos )
            mg.parameters->kinematics.mode = CepGen::Kinematics::ElasticElastic;
          else if ( kin_mode.find( "singlediss" ) != string::npos )
            mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticElastic;
          else if ( kin_mode.find( "doublediss" ) != string::npos )
            mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticInelastic;
          else {
            InError( Form( "Unrecognized kinematics mode: %s", values_vs_kin.first ) );
            break;
          }

          if ( kin_mode.find( "_su" ) != string::npos )
            mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SzczurekUleshchenko;
          else if ( kin_mode.find( "_lux" ) != string::npos )
            mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::Schaefer;
          else if ( kin_mode.find( "_allm" ) != string::npos )
            mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::ALLM97;
          else
            mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SuriYennie;

          //mg.parameters->dump();
          Information( Form( "Process: %s/%s\n\t"
                             "Configuration time: %.3f ms",
                             values_vs_generator.first,
                             values_vs_kin.first,
                             tmr.elapsed()*1.e3 ) );
          tmr.reset();

          mg.clearRun();
          const double xsec_ref = values_vs_kin.second.first;
          const double err_xsec_ref = values_vs_kin.second.second;
          double xsec_cepgen, err_xsec_cepgen;
          mg.computeXsection( xsec_cepgen, err_xsec_cepgen );

          const double sigma = fabs( xsec_ref-xsec_cepgen ) / std::hypot( err_xsec_cepgen, err_xsec_ref );

          Information( Form( "Computed cross section:\n\t"
                             "Ref.   = %.3e +/- %.3e\n\t"
                             "CepGen = %.3e +/- %.3e\n\t"
                             "Pull: %.6f",
                             xsec_ref, err_xsec_ref,
                             xsec_cepgen, err_xsec_cepgen,
                             sigma ) );

          Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
          tmr.reset();

          ostringstream oss; oss << values_vs_kin.first;
          string test_res = Form( "%-10s", values_vs_generator.first )+"\t"+
                            Form( "pt-gt-%.1f", values_vs_cut.first )+"\t"+
                            Form( "%-16s", oss.str().c_str() )+"\t"
                            "ref="+Form( "%g", xsec_ref )+"\t"
                            "got="+Form( "%g", xsec_cepgen )+"\t"
                            "pull="+Form( "%+g", sigma );
          if ( fabs( sigma ) < num_sigma ) {
            passed_tests.emplace_back( test_res );
            num_tests_passed++;
          }
          else
            failed_tests.emplace_back( test_res );

          num_tests++;
          cout << "Test " << num_tests_passed << "/"
                          << num_tests << " passed!" << endl;
        }
      }
    }
  } catch ( CepGen::Exception& e ) {}
  if ( failed_tests.size() != 0 ) {
    ostringstream os_failed, os_passed;
    for ( const auto& fail : failed_tests )
      os_failed << " (*) " << fail << endl;
    for ( const auto& pass : passed_tests )
      os_passed << " (*) " << pass << endl;
    throw CepGen::Exception( __PRETTY_FUNCTION__,
      Form( "Some tests failed!\n"
            "%s\n"
            "Passed tests:\n"
            "%s",
            os_failed.str().c_str(),
            os_passed.str().c_str() ), CepGen::FatalError );
  }

  Information( "ALL TESTS PASSED!" );

  return 0;
}
