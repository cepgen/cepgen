#include "CepGen/Generator.h"
#include <assert.h>

int
main( int argc, char* argv[] )
{
  if ( argc<3 ) {
    InError( Form( "Usage: %s [process] [process_type]", argv[0] ) );
    exit( 0 );
  }

  typedef std::map<CepGen::Kinematics::ProcessMode,std::pair<double,double> > KinematicsMap;

  // values defined at pt(single lepton)>15 GeV, |eta(single lepton)|<2.5, mX<1000 GeV
  std::map<std::string,KinematicsMap> values_map = {
    //--- LPAIR values
    { "lpair",
      {
        { CepGen::Kinematics::ElasticElastic,     { 4.16891e-1, 7.53702e-4 } },
        { CepGen::Kinematics::ElasticInelastic,   { 4.87144e-1, 1.05158e-3 } },
        { CepGen::Kinematics::InelasticInelastic, { 6.35650e-1, 1.93968e-3 } }
      }
    }
    //--- PPTOLL values
    //{ "pptoll", {} }
  };

  double xsec, err;

  Timer tmr;
  {
    CepGen::Generator mg;
    if      ( strcmp( argv[1], "lpair"  )==0 ) mg.parameters->process = new CepGen::Process::GamGamLL;
    else if ( strcmp( argv[1], "pptoll" )==0 ) mg.parameters->process = new CepGen::Process::PPtoLL;

    if      ( strcmp( argv[2], "elastic"    )==0 ) mg.parameters->process_mode = CepGen::Kinematics::ElasticElastic;
    else if ( strcmp( argv[2], "singlediss" )==0 ) mg.parameters->process_mode = CepGen::Kinematics::ElasticInelastic;
    else if ( strcmp( argv[2], "doublediss" )==0 ) mg.parameters->process_mode = CepGen::Kinematics::InelasticInelastic;

    mg.parameters->minpt = 15.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;

    Information( Form( "Configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
    {
      tmr.reset();
      const std::pair<double,double> orig_value = values_map[argv[1]][mg.parameters->process_mode];
      mg.computeXsection( xsec, err );

      const double sigma = ( fabs( orig_value.first-xsec )-orig_value.second)/err;
      std::cout << sigma << std::endl;

      Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
      assert( fabs( sigma )<1.0 );
    }
  }

  return 0;
}
