/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Modules/SurvivalFactorFactory.h"
#include "CepGen/SurvivalFactors/SurvivalFactor.h"
#include "CepGen/Utils/GridHandler.h"

using namespace std::string_literals;

namespace cepgen::survival {
  /// Survival factor interpolated from a one-dimensional grid of values
  /// \note Input file format: two whitespace-separated columns (variable value, survival factor), '#' for comments.
  ///  Outside the grid range, the value at the closest boundary is used.
  class Grid1D final : public SurvivalFactor {
  public:
    explicit Grid1D(const ParametersList& params)
        : SurvivalFactor(params), variable_(steer<std::string>("variable")), grid_(GridType::linear) {
      const auto path = steer<std::string>("path");
      std::ifstream file(path);
      if (!file.is_open())
        throw CG_FATAL("SurvivalFactor:grid1D") << "Failed to open survival factors grid file '" << path << "'.";
      std::string line;
      size_t num_points = 0;
      while (std::getline(file, line)) {
        if (const auto pos = line.find('#'); pos != std::string::npos)
          line.erase(pos);
        std::istringstream iss(line);
        double coord, value;
        if (!(iss >> coord >> value))
          continue;
        grid_.insert({coord}, {value});
        ++num_points;
      }
      if (num_points < 2)
        throw CG_FATAL("SurvivalFactor:grid1D")
            << "At least two points are required in grid file '" << path << "', " << num_points << " found.";
      grid_.initialise();
      limits_ = grid_.boundaries().at(0);
    }

    static ParametersDescription description() {
      auto desc = SurvivalFactor::description();
      desc.setDescription("1D grid-interpolated survival factor");
      desc.add("path", ""s).setDescription("path to the two-column (variable, S^2) grid file");
      desc.add("variable", "m(4)"s).setDescription("event variable (event browser syntax) the grid depends on");
      return desc;
    }

    double operator()(const Event& event) const override {
      return grid_.eval({limits_.trim(browser_.get(event, variable_))}).at(0);
    }

  private:
    const std::string variable_;
    const utils::EventBrowser browser_;
    GridHandler<1, 1> grid_;
    Limits limits_;
  };
}  // namespace cepgen::survival
using cepgen::survival::Grid1D;
REGISTER_SURVIVAL_FACTOR("grid1D", Grid1D);
