#include "CepGen/Generator.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"

#include <assert.h>

int
main( int argc, char* argv[] )
{
  typedef std::map<std::string,std::pair<double,double> > KinematicsMap;
  typedef std::map<float, KinematicsMap> ValuesAtCutMap;

  // values defined at pt(single lepton)>15 GeV, |eta(single lepton)|<2.5, mX<1000 GeV
  // process -> { pt cut -> { kinematics -> ( sigma, delta(sigma) ) } }
  std::map<std::string,ValuesAtCutMap> values_map = {
    //--- LPAIR values at sqrt(s) = 13 TeV
    { "1_lpair", {
      { 15.0, { // pt cut
          { "1_elastic",    { 4.1994803e-1, 8.328e-4 } },
          { "2_singlediss", { 4.8504819e-1, 1.171e-3 } },
          { "3_doublediss", { 6.35650e-1, 1.93968e-3 } }
      } }
    } },
    //--- PPTOLL values
    { "2_pptoll", { {} } },
  };

  const double num_sigma = 3.0;

  CepGen::Logger::get().level = CepGen::Logger::Nothing;

  Timer tmr;
  CepGen::Generator mg;

  mg.parameters->kinematics.setSqrtS( 13.e3 );
  mg.parameters->kinematics.eta_min = -2.5;
  mg.parameters->kinematics.eta_max = 2.5;
  mg.parameters->kinematics.mx_max = 1000.;
  mg.parameters->vegas.ncvg = 50000;
  mg.parameters->vegas.itvg = 5;

  Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  unsigned short num_tests = 0, num_tests_passed = 0;

  for ( const auto& values_vs_generator : values_map ) { // loop over all generators
    if      ( values_vs_generator.first == "1_lpair"  ) mg.parameters->setProcess( new CepGen::Process::GamGamLL );
    else if ( values_vs_generator.first == "2_pptoll" ) mg.parameters->setProcess( new CepGen::Process::PPtoLL );
    else { InError( Form( "Unrecognized generator mode: %s", values_vs_generator.first.c_str() ) ); break; }

    for ( const auto& values_vs_cut : values_vs_generator.second ) { // loop over the single lepton pT cut
      mg.parameters->kinematics.pt_min = values_vs_cut.first;
      for ( const auto& values_vs_kin : values_vs_cut.second ) { // loop over all possible kinematics
        if      ( values_vs_kin.first == "1_elastic"    ) mg.parameters->kinematics.mode = CepGen::Kinematics::ElasticElastic;
        else if ( values_vs_kin.first == "2_singlediss" ) mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticElastic;
        else if ( values_vs_kin.first == "3_doublediss" ) mg.parameters->kinematics.mode = CepGen::Kinematics::InelasticInelastic;
        else { InError( Form( "Unrecognized kinematics mode: %s", values_vs_kin.first.c_str() ) ); break; }

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

        if ( fabs( sigma )<num_sigma ) num_tests_passed++;
        else throw CepGen::Exception( __PRETTY_FUNCTION__, Form( "Test %s/%s failed!", values_vs_generator.first.c_str(), values_vs_kin.first.c_str() ), CepGen::FatalError );
        num_tests++;
        std::cout << "Test " << num_tests_passed << "/" << num_tests << " passed!" << std::endl;
      }
    }
  }

  Information( "ALL TESTS PASSED!" );

  return 0;
}
