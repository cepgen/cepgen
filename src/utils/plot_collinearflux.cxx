/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  vector<int> modes;
  int strfun_type, num_points;
  double mx, xmin, xmax;
  string ffmode, output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("formfac,f", "form factors modelling", &ffmode, ffmode)
      .addOptionalArgument("modes,t", "beam modelling(s)", &modes, vector<int>{(int)cepgen::mode::Beam::ProtonElastic})
      .addOptionalArgument("mx,M", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("xmin,x", "minimal fractional loss", &xmin, 0.)
      .addOptionalArgument("xmax,X", "maximal fractional loss", &xmax, 1.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.scan.output.txt")
      .parse();

  cepgen::initialise();

  ofstream out(output_file);

  out << "# struct. functions: " << strfun_type << "\n"
      << "# form factors: " << ffmode << "\n"
      << "# diffractive mass: " << mx << " GeV/c2\n"
      << "# fractional momentum loss: " << cepgen::Limits(xmin, xmax) << "\n"
      << "# fluxes modes:";
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  vector<cepgen::KTFlux> ktfluxes;
  for (const auto& mode : modes) {
    switch ((cepgen::mode::Beam)mode) {
      case cepgen::mode::Beam::ProtonElastic:
        ktfluxes.emplace_back(cepgen::KTFlux::P_Photon_Elastic);
        break;
      case cepgen::mode::Beam::ProtonInelastic:
        ktfluxes.emplace_back(cepgen::KTFlux::P_Photon_Inelastic);
        break;
      default:
        throw CG_FATAL("main") << "Invalid beam mode: " << mode << "!";
    }
    out << "\t" << (cepgen::mode::Beam)mode;
  }
  out << "\n";
  auto ff = cepgen::formfac::FormFactorsFactory::get().build(ffmode);
  ff->setStructureFunctions(sf.get());

  const cepgen::Limits kt2_limits(0., 10000.);

  const cepgen::CollinearFlux flux(ff.get(), kt2_limits);
  for (int i = 0; i < num_points; ++i) {
    const double x = xmin + i * (xmax - xmin) / (num_points - 1);
    out << x;
    for (const auto& ktflux : ktfluxes)
      out << "\t" << flux(x, mx, ktflux);
    out << "\n";
  }
  out.close();

  return 0;
}
