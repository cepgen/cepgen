#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include <fstream>

int main( int argc, char* argv[] )
{
  int formfac_type, strfun_type, num_points;
  double kt2, mx;
  std::string output_file;

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "ff,f", "form factors modelling", &formfac_type, 1 )
    .addOptionalArgument( "sf,s", "structure functions modelling", &strfun_type, 301 )
    .addOptionalArgument( "kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 100. )
    .addOptionalArgument( "mx,m", "diffractive state mass (GeV)", &mx, 1.5 )
    .addOptionalArgument( "npoints,n", "number of x-points to scan", &num_points, 100 )
    .addOptionalArgument( "output,o", "output file name", &output_file, "flux.scan.output.txt" )
    .parse();

  cepgen::initialise();
  const double mi = cepgen::PDG::get().mass( cepgen::PDG::proton );
  const double mi2 = mi*mi, mx2 = mx*mx;

  auto ff = cepgen::formfac::FormFactorsFactory::get().build( formfac_type );
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build( strfun_type );
  ff->setStructureFunctions( sf.get() );
  std::ofstream out( output_file );
  out
    << "# form factors: " << ff.get() << "\n"
    << "# structure functions: " << sf.get() << "\n"
    << "# kt2 = " << kt2 << " GeV^2\n"
    << "# mX = " << mx << " GeV\n";
  for ( int i = 0; i < num_points; ++i ) {
    const double x = i*1./num_points;
    out << x
      << "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Photon_Elastic, x, kt2, *ff, mi2, mx2 )
      << "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Photon_Inelastic_Budnev, x, kt2, *ff, mi2, mx2 )
      //<< "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Gluon_KMR, x, kt2, *ff, mi2, mx2 )
      //<< "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Gluon_KMR_alt, x, kt2, *ff, mi2, mx2 )
      << "\n";
  }
  CG_INFO( "main" )
    << "Scan written in \"" << output_file << "\".";
  out.close();

  return 0;
}
