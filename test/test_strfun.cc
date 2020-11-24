#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/FormFactors/Parameterisation.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include <fstream>

int main( int argc, char* argv[] )
{
  int strfun_type, num_points;
  std::vector<double> q2in, xbjin;
  std::string output_file;

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "sf,s", "structure functions modelling", &strfun_type, 301 )
    .addOptionalArgument( "q2,q", "parton virtuality (GeV^2)", &q2in, std::vector<double>{} )
    .addOptionalArgument( "xbj,x", "Bjorken-x", &xbjin, std::vector<double>{} )
    .addOptionalArgument( "npoints,n", "number of x-points to scan", &num_points, 100 )
    .addOptionalArgument( "output,o", "output file name", &output_file, "flux.scan.output.txt" )
    .parse();

  cepgen::initialise();
  std::vector<double> q2vals, xbjvals;

  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build( strfun_type );
  std::ofstream out( output_file );
  out << "# structure functions: " << sf.get() << "\n";
  if ( q2in.empty() )
    throw CG_FATAL( "main" ) << "At least one value of Q^2 is required!";
  else if ( q2in.size() == 2 ) { // min-max
    q2vals.clear();
    for ( int i = 0; i <= num_points; ++i )
      q2vals.emplace_back( q2in[0]+i*( q2in[1]-q2in[0] )/num_points );
  }
  else
    q2vals = q2in;

  if ( xbjin.empty() )
    throw CG_FATAL( "main" ) << "At least one value of x_Bj is required!";
  else if ( xbjin.size() == 2 ) { // min-max
    xbjvals.clear();
    for ( int i = 0; i <= num_points; ++i )
      xbjvals.emplace_back( xbjin[0]+i*( xbjin[1]-xbjin[0] )/num_points );
  }
  else
    xbjvals = xbjin;

  out << "# q2\txbj\tF_2\tF_L\n";

  for ( const auto& xbj : xbjvals )
    for ( const auto& q2 : q2vals ) {
      auto& sfval = (*sf)( xbj, q2 );
      sfval.computeFL( xbj, q2 );
      out << q2 << "\t" << xbj << "\t" << sfval.F2 << "\t" << sfval.FL << "\n";
    }

  CG_INFO( "main" )
    << "Scan written in \"" << output_file << "\".";
  out.close();

  return 0;
}
