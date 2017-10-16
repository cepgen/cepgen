#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Core/Timer.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "abort.h"

#include <unordered_map>
#include <assert.h>

int
main( int argc, char* argv[] )
{
  typedef std::vector<std::pair<std::string,std::pair<double,double> > > KinematicsMap;
  typedef std::vector<std::pair<float, KinematicsMap> > ValuesAtCutMap;

  AbortHandler ctrl_c;

  // values defined at pt(single lepton)>15 GeV, |eta(single lepton)|<2.5, mX<1000 GeV
  // process -> { pt cut -> { kinematics -> ( sigma, delta(sigma) ) } }
  std::vector<std::pair<std::string,ValuesAtCutMap> > values_map = {
    //--- LPAIR values at sqrt(s) = 13 TeV
    { "lpair", {
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
    } },
    //--- PPTOLL values
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

  if ( argc == 1 || strcmp( argv[1], "debug" ) != 0 ) {
    CepGen::Logger::get().level = CepGen::Logger::Nothing;
  }

  Timer tmr;
  CepGen::Generator mg;

  mg.parameters->kinematics.setSqrtS( 13.e3 );
  mg.parameters->kinematics.central_cuts[CepGen::Cuts::eta_single].in( -2.5, 2.5 );
  mg.parameters->kinematics.remnant_cuts[CepGen::Cuts::mass].max() = 1000.;
  //mg.parameters->vegas.ncvg = 50000;
  mg.parameters->vegas.itvg = 5;

  Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  unsigned short num_tests = 0, num_tests_passed = 0;
  std::vector<std::string> failed_tests, passed_tests;

  try {
    for ( const auto& values_vs_generator : values_map ) { // loop over all generators
      if      ( values_vs_generator.first == "lpair"  ) mg.parameters->setProcess( new CepGen::Process::GamGamLL );
      else if ( values_vs_generator.first == "pptoll" ) {
        mg.parameters->setProcess( new CepGen::Process::PPtoLL );
        mg.parameters->kinematics.initial_cuts[CepGen::Cuts::qt].max() = 50.0;
      }
      else { InError( Form( "Unrecognized generator mode: %s", values_vs_generator.first.c_str() ) ); break; }

      for ( const auto& values_vs_cut : values_vs_generator.second ) { // loop over the single lepton pT cut
        mg.parameters->kinematics.central_cuts[CepGen::Cuts::pt_single].min() = values_vs_cut.first;
        for ( const auto& values_vs_kin : values_vs_cut.second ) { // loop over all possible kinematics
          if      ( values_vs_kin.first.find( "elastic"    ) != std::string::npos ) mg.parameters->kinematics.mode = CepGen::Kinematics::ElasticElastic;
          else if ( values_vs_kin.first.find( "singlediss" ) != std::string::npos ) mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticElastic;
          else if ( values_vs_kin.first.find( "doublediss" ) != std::string::npos ) mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticInelastic;
          else { InError( Form( "Unrecognized kinematics mode: %s", values_vs_kin.first.c_str() ) ); break; }

          if ( values_vs_kin.first.find( "_su" ) != std::string::npos ) mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SzczurekUleshchenko;
          else mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SuriYennie;

          Information( Form( "Process: %s/%s\n\tConfiguration time: %.3f ms", values_vs_generator.first.c_str(), values_vs_kin.first.c_str(), tmr.elapsed()*1.e3 ) );
          tmr.reset();

          mg.clearRun();
          const double xsec_ref = values_vs_kin.second.first, err_xsec_ref = values_vs_kin.second.second;
          double xsec_cepgen, err_xsec_cepgen;
          mg.computeXsection( xsec_cepgen, err_xsec_cepgen );

          const double sigma = ( fabs( xsec_ref-xsec_cepgen ) ) / sqrt( err_xsec_cepgen*err_xsec_cepgen + err_xsec_ref*err_xsec_ref );

          Information( Form( "Computed cross section:\n\tRef.   = %.3e +/- %.3e\n\tCepGen = %.3e +/- %.3e\n\tPull: %.6f", xsec_ref, err_xsec_ref, xsec_cepgen, err_xsec_cepgen, sigma ) );

          Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
          tmr.reset();

          std::string test_res = values_vs_generator.first+":"+
                                 Form( "pt-gt-%.1f", values_vs_cut.first )+":"+
                                 values_vs_kin.first+":"
                                 "ref="+Form( "%.3e", xsec_ref )+":"
                                 "got="+Form( "%.3e", xsec_cepgen )+":"
                                 "pull="+Form( "%.3f", sigma );
          if ( fabs( sigma )<num_sigma ) {
            passed_tests.emplace_back( test_res );
            num_tests_passed++;
          }
          else {
            failed_tests.emplace_back( test_res );
          }
          num_tests++;
          std::cout << "Test " << num_tests_passed << "/" << num_tests << " passed!" << std::endl;
        }
      }
    }
  } catch ( CepGen::Exception& e ) {}
  if ( failed_tests.size() != 0 ) {
    std::ostringstream os_failed, os_passed;
    for ( const auto& fail : failed_tests ) os_failed << " (*) " << fail << std::endl;
    for ( const auto& pass : passed_tests ) os_passed << " (*) " << pass << std::endl;
    throw CepGen::Exception( __PRETTY_FUNCTION__, Form( "Some tests failed!\n%s\nPassed tests:\n%s", os_failed.str().c_str(), os_passed.str().c_str() ), CepGen::FatalError );
  }
  else Information( "ALL TESTS PASSED!" );

  return 0;
}
