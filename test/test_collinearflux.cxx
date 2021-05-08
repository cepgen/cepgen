#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Modes.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/FormFactors/Parameterisation.h"

#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
  int mode, strfun_type, num_points;
  double mx, xmin, xmax;
  string ffmode, output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("formfac,f", "form factors modelling", &ffmode, ffmode)
      .addArgument("mode,t", "beam modelling", &mode, (int)cepgen::mode::Beam::ProtonElastic)
      .addOptionalArgument("mx,M", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("xmin,x", "minimal fractional loss", &xmin, 0.)
      .addOptionalArgument("xmax,X", "maximal fractional loss", &xmax, 1.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.scan.output.txt")
      .parse();

  cepgen::initialise();

  ofstream out(output_file);

  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  cepgen::KTFlux ktflux;
  switch ((cepgen::mode::Beam)mode) {
    case cepgen::mode::Beam::ProtonElastic:
      ktflux = cepgen::KTFlux::P_Photon_Elastic;
      break;
    case cepgen::mode::Beam::ProtonInelastic:
      ktflux = cepgen::KTFlux::P_Photon_Inelastic;
      break;
    default:
      throw CG_FATAL("main") << "Invalid beam mode: " << mode << "!";
  }
  auto ff = cepgen::formfac::FormFactorsFactory::get().build(ffmode);
  ff->setStructureFunctions(sf.get());

  cepgen::CollinearFlux flux(ktflux, cepgen::Limits(0., 10000.), ff.get());
  for (int i = 0; i < num_points; ++i) {
    const double x = xmin + i * (xmax - xmin) / (num_points - 1);
    out << x << "\t" << flux(x, mx) << "\n";
  }
  out.close();

  return 0;
}
