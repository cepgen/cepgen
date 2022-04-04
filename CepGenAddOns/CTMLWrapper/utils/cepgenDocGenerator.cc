/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Version.h"
#include "CepGenAddOns/CTMLWrapper/DocumentationGenerator.h"

namespace cepgen {
  namespace utils {}
}  // namespace cepgen

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::utils::DocumentationGenerator gen{cepgen::ParametersList()};

  //gen.document("Cards handler", cepgen::card::CardsHandlerFactory::get());
  gen.document("Form factors", cepgen::formfac::FormFactorsFactory::get());
  gen.document("Structure functions", cepgen::strfun::StructureFunctionsFactory::get());
  gen.document("Longitudinal/transverse cross section ratio parameterisations",
               cepgen::sigrat::SigmaRatiosFactory::get());

  return 0;
}
