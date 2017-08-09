#include "CepGen/Generator.h"
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
    { "lpair",
      {
        {
          15.0, { // pt cut
            { "elastic",    { 4.16891e-1, 7.53702e-4 } },
            { "singlediss", { 4.87144e-1, 1.05158e-3 } },
            { "doublediss", { 6.35650e-1, 1.93968e-3 } }
          }
        }
      }
    },
    //--- PPTOLL values
    //{ "pptoll", { {} } }
  };

  const double num_sigma = 3.0;

  Timer tmr;
  {
    CepGen::Generator mg;

    //Logger::GetInstance()->Level = Logger::Debug;

    mg.parameters->setSqrtS( 13.e3 );
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    mg.parameters->ncvg = 50000;
    mg.parameters->itvg = 5;

    Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
    tmr.reset();
    for ( const auto& values_vs_generator : values_map ) { // loop over all generators
      if      ( values_vs_generator.first.compare( "lpair"  ) == 0 ) mg.parameters->setProcess( new CepGen::Process::GamGamLL );
      else if ( values_vs_generator.first.compare( "pptoll" ) == 0 ) mg.parameters->setProcess( new CepGen::Process::PPtoLL );
      else { InError( Form( "Unrecognized generator mode: %s", values_vs_generator.first.c_str() ) ); break; }

      for ( const auto& values_vs_cut : values_vs_generator.second ) { // loop over the single lepton pT cut
        mg.parameters->minpt = values_vs_cut.first;
        for ( const auto& values_vs_kin : values_vs_cut.second ) { // loop over all possible kinematics
          if      ( values_vs_kin.first.compare( "elastic"    ) == 0 ) mg.parameters->process_mode = CepGen::Kinematics::ElasticElastic;
          else if ( values_vs_kin.first.compare( "singlediss" ) == 0 ) mg.parameters->process_mode = CepGen::Kinematics::ElasticInelastic;
          else if ( values_vs_kin.first.compare( "doublediss" ) == 0 ) mg.parameters->process_mode = CepGen::Kinematics::InelasticInelastic;
          else { InError( Form( "Unrecognized kinematics mode: %s", values_vs_kin.first.c_str() ) ); break; }

          Information( Form( "Process: %s/%s ; configuration time: %.3f ms", values_vs_generator.first.c_str(), values_vs_kin.first.c_str(), tmr.elapsed()*1.e3 ) );

          const float xsec_ref = values_vs_kin.second.first, err_xsec_ref = values_vs_kin.second.second;

          double xsec_cepgen, err_xsec_cepgen;
          mg.computeXsection( xsec_cepgen, err_xsec_cepgen );

          const double sigma = ( fabs( xsec_ref-xsec_cepgen ) ) / sqrt( err_xsec_cepgen*err_xsec_cepgen + err_xsec_ref*err_xsec_ref );
          std::cout << sigma << ":::" << ( xsec_ref-xsec_cepgen ) << std::endl;

          Information( Form( "Computed cross section:\n\tRef.   = %.2e +/- %.2e\n\tCepGen = %.2e +/- %.2e", xsec_ref, err_xsec_ref, xsec_cepgen, err_xsec_cepgen ) );
          Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
          assert( fabs( sigma )<num_sigma );
        }
      }
    }
  }

  return 0;
}
