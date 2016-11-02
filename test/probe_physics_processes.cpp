#include "core/MCGen.h"
#include <assert.h>

int
main( int argc, char* argv[] )
{
  if ( argc<3 ) {
    InError( Form( "Usage: %s [process] [process_type]", argv[0] ) );
    exit( 0 );
  }

  typedef std::map< Kinematics::ProcessMode,std::pair<double,double> > KinematicsMap;

  // values defined at pt(single lepton)>15 GeV, |eta(single lepton)|<2.5, mX<1000 GeV
  std::map<std::string,KinematicsMap> values_map;

  //--- LPAIR values
  KinematicsMap lpair_values;
  lpair_values[Kinematics::ElasticElastic] = std::make_pair( 4.16891e-1, 7.53702e-4 );
  lpair_values[Kinematics::ElasticInelastic] = std::make_pair( 4.87144e-1, 1.05158e-3 );
  lpair_values[Kinematics::InelasticInelastic] = std::make_pair( 6.3565e-1, 1.93968e-3 );
  values_map["lpair"] = lpair_values;

  //--- PPTOLL values
  KinematicsMap pptoll_values;
  values_map["pptoll"] = pptoll_values;

  double xsec, err;

  Timer tmr;
  {
    MCGen mg;
    if ( strcmp( argv[1], "lpair" )==0 )       mg.parameters->process = new GamGamLL;
    else if ( strcmp( argv[1], "pptoll" )==0 ) mg.parameters->process = new PPtoLL;

    if ( strcmp( argv[2], "elastic" )==0 )         mg.parameters->process_mode = Kinematics::ElasticElastic;
    else if ( strcmp( argv[2], "singlediss" )==0 ) mg.parameters->process_mode = Kinematics::ElasticInelastic;
    else if ( strcmp( argv[2], "doublediss" )==0 ) mg.parameters->process_mode = Kinematics::InelasticInelastic;

    mg.parameters->minpt = 15.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;

    Information( Form( "Configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
    {
      tmr.reset();
      const std::pair<double,double> orig_value = values_map[argv[1]][mg.parameters->process_mode];
      mg.ComputeXsection( &xsec, &err );

      const double sigma = ( fabs( orig_value.first-xsec )-orig_value.second)/err;
      std::cout << sigma << std::endl;

      Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
      assert( fabs( sigma )<1.0 );
    }
  }

  return 0;
}
