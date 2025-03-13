/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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

#include <cmath>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Utils/Timer.h"

using namespace kmr;

GluonGrid& GluonGrid::get(const cepgen::ParametersList& params) {
  static GluonGrid instance(!params.empty() ? params : description().parameters());
  return instance;
}

GluonGrid::GluonGrid(const cepgen::ParametersList& params)
    : GridHandler<3, 1>(cepgen::GridType::linear /*grid is already logarithmic*/),
      SteeredObject(params),
      grid_path_(steerPath("path")) {
  CG_INFO("GluonGrid") << "Building the KMR grid evaluator.";

  cepgen::utils::Timer tmr;
  {  // file readout part
    std::ifstream file(grid_path_, std::ios::in);
    if (!file.is_open())
      throw CG_FATAL("GluonGrid") << "Failed to load grid file \"" << grid_path_ << "\"!";

    std::string x, kt2, mu2, fg;
    while (file >> x >> kt2 >> mu2 >> fg)
      insert({std::stod(x), std::stod(kt2), std::stod(mu2)}, {std::stod(fg)});
    file.close();
    initialise();  // initialise the grid after filling its nodes
  }
  const auto limits = boundaries();
  CG_INFO("GluonGrid") << "KMR grid evaluator built in " << tmr.elapsed() << " s.\n\t"
                       << " log(x)    in range " << limits.at(0) << ",\t"
                       << "x    in range " << limits.at(0).compute(std::exp) << "\n\t"
                       << " log(kt^2) in range " << limits.at(1) << ",\t"
                       << "kt^2 in range " << limits.at(1).compute(std::exp) << "\n\t"
                       << " log(mu^2) in range " << limits.at(2) << ",\t"
                       << "mu^2 in range " << limits.at(2).compute(std::exp) << ".";
}

double GluonGrid::operator()(double x, double kt2, double mu2) const {
  return eval({std::log10(x), std::log10(kt2), std::log10(mu2)}).at(0);
}

cepgen::ParametersDescription GluonGrid::description() {
  auto desc = cepgen::ParametersDescription();
  desc.add<std::string>("path", DEFAULT_KMR_GRID_PATH);
  return desc;
}
