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
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"

int main(int argc, char* argv[]) {
  std::string formfac_type;
  int strfun_type, num_points;
  double kt2, mx;
  std::string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("ff,f", "form factors modelling", &formfac_type, "StandardDipole")
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 100.)
      .addOptionalArgument("mx,m", "diffractive state mass (GeV)", &mx, 1.5)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .parse();

  cepgen::initialise();
  const double mi = cepgen::PDG::get().mass(cepgen::PDG::proton);
  const double mi2 = mi * mi, mx2 = mx * mx;

  auto ff = cepgen::formfac::FormFactorsFactory::get().build(formfac_type);
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  ff->setStructureFunctions(sf.get());
  std::ofstream out(output_file);
  out << "# form factors: " << ff.get() << "\n"
      << "# structure functions: " << sf.get() << "\n"
      << "# kt2 = " << kt2 << " GeV^2\n"
      << "# mX = " << mx << " GeV\n";
  for (int i = 0; i < num_points; ++i) {
    const double x = i * 1. / num_points;
    out << x << "\t" << cepgen::ktFlux(cepgen::KTFlux::P_Photon_Elastic, x, kt2, *ff, mi2, mx2) << "\t"
        << cepgen::ktFlux(cepgen::KTFlux::P_Photon_Inelastic_Budnev, x, kt2, *ff, mi2, mx2)
        //<< "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Gluon_KMR, x, kt2, *ff, mi2, mx2 )
        //<< "\t" << cepgen::ktFlux( cepgen::KTFlux::P_Gluon_KMR_alt, x, kt2, *ff, mi2, mx2 )
        << "\n";
  }
  CG_INFO("main") << "Scan written in \"" << output_file << "\".";
  out.close();

  return 0;
}
