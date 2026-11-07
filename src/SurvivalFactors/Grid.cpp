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

#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Modules/SurvivalFactorFactory.h"
#include "CepGen/SurvivalFactors/SurvivalFactor.h"
#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/String.h"

namespace cepgen::survival {
  /// Survival factor interpolated from a D-dimensional grid of values, as a function of D event-level variables
  /// \note Input file format: D+1 whitespace-separated columns (D coordinates, survival factor), '#' for comments.
  ///  Comment lines of the form "#@ key value" are collected as metadata; a "#@ sqrtS <value>" line (e.g. as written
  ///  by the cepgenSurvivalGrid utility) is checked against the beams centre-of-mass energy.
  ///  Variables follow the event browser syntax, e.g. pt(ob1), xi(ob2), acop(ob1,ob2), or m(int).
  ///  1D grids are interpolated with cubic splines, higher-dimensional grids multilinearly (see GridHandler).
  ///  Outside the grid ranges, either the value at the closest boundary is used, or a fallback survival factor
  ///  modelling is called when one is specified.
  template <size_t D>
  class Grid final : public SurvivalFactor {
  public:
    explicit Grid(const ParametersList& params)
        : SurvivalFactor(params),
          grid_(GridType::linear),
          path_(steer<std::string>("path")),
          variables_(steer<std::vector<std::string> >("variables")) {
      if constexpr (D == 1)  // backward-compatible single-variable steering
        if (variables_.empty())
          variables_ = {steer<std::string>("variable")};
      if (variables_.size() != D)
        throw CG_FATAL("SurvivalFactor:grid")
            << D << "D grid requires " << D << " event variable(s), " << variables_.size() << " given: {"
            << utils::merge(variables_, ", ") << "}.";
      readGrid();
      limits_ = grid_.boundaries();
      if (const auto& fallback = steer<ParametersList>("fallback"); fallback.hasName())
        fallback_ = SurvivalFactorFactory::get().build(fallback);
      CG_DEBUG("SurvivalFactor:grid").log([&](auto& log) {
        log << D << "D survival factors grid loaded from '" << path_ << "':";
        for (size_t i = 0; i < D; ++i)
          log << "\n\t" << variables_[i] << " in " << limits_[i];
      });
    }

    static ParametersDescription description() {
      auto desc = SurvivalFactor::description();
      desc.setDescription(std::to_string(D) + "D grid-interpolated survival factor");
      desc.add("path", std::string{}).setDescription("path to the (D+1)-column (variables, S^2) grid file");
      if constexpr (D == 1)
        desc.add("variable", std::string{"m(4)"})
            .setDescription("event variable (event browser syntax) the grid depends on");
      desc.add("variables", std::vector<std::string>{})
          .setDescription("event variables (event browser syntax) the grid depends on");
      desc.add("fallback", ParametersDescription{})
          .setDescription("survival factor modelling used outside the grid ranges (none: closest boundary value)");
      desc.add("energyTolerance", 1.e-3)
          .setDescription("maximum relative difference between the beams and grid (if specified) energies");
      return desc;
    }

    double operator()(const Event& event) const override {
      checkEnergy(event);
      std::vector<double> coord(D);
      bool in_range = true;
      for (size_t i = 0; i < D; ++i) {
        coord[i] = browser_.get(event, variables_[i]);
        in_range = in_range && limits_[i].contains(coord[i]);
      }
      if (!in_range) {
        if (fallback_)
          return (*fallback_)(event);
        for (size_t i = 0; i < D; ++i)
          coord[i] = limits_[i].trim(coord[i]);
      }
      return grid_.eval(coord).at(0);
    }

  private:
    void readGrid() {
      std::ifstream file(path_);
      if (!file.is_open())
        throw CG_FATAL("SurvivalFactor:grid") << "Failed to open survival factors grid file '" << path_ << "'.";
      std::string line;
      size_t num_points = 0;
      while (std::getline(file, line)) {
        if (line.rfind("#@", 0) == 0) {  // metadata line
          std::istringstream iss(line.substr(2));
          std::string key, value;
          if (iss >> key && std::getline(iss >> std::ws, value))
            metadata_[key] = value;
          continue;
        }
        if (const auto pos = line.find('#'); pos != std::string::npos)
          line.erase(pos);
        std::istringstream iss(line);
        std::vector<double> coord(D);
        double value;
        bool valid = true;
        for (auto& c : coord)
          valid = valid && static_cast<bool>(iss >> c);
        if (!valid || !(iss >> value))
          continue;
        grid_.insert(coord, {value});
        ++num_points;
      }
      if (num_points < (1ul << D))
        throw CG_FATAL("SurvivalFactor:grid") << "At least " << (1ul << D) << " points are required in " << D
                                              << "D grid file '" << path_ << "', " << num_points << " found.";
      grid_.initialise();
    }
    void checkEnergy(const Event& event) const {
      if (energy_checked_)
        return;
      if (const auto it = metadata_.find("sqrtS"); it != metadata_.end()) {
        const auto grid_sqrt_s = std::stod(it->second), sqrt_s = event.cmEnergy();
        if (std::fabs(sqrt_s / grid_sqrt_s - 1.) > steer<double>("energyTolerance"))
          throw CG_FATAL("SurvivalFactor:grid") << "Grid '" << path_ << "' was computed for sqrt(s) = " << grid_sqrt_s
                                                << " GeV, but the beams have sqrt(s) = " << sqrt_s << " GeV.";
      }
      energy_checked_ = true;
    }

    GridHandler<D, 1> grid_;
    const std::string path_;
    const utils::EventBrowser browser_;

    std::vector<std::string> variables_;
    std::array<Limits, D> limits_;
    std::map<std::string, std::string> metadata_;
    std::unique_ptr<SurvivalFactor> fallback_;
    mutable bool energy_checked_{false};
  };

  // dimensionalities of interest (each requires the corresponding GridHandler<D, 1> instantiation)
  using Grid1D = Grid<1>;  ///< e.g. S^2 vs. the two-parton system mass
  using Grid2D = Grid<2>;
  using Grid3D = Grid<3>;
  using Grid4D = Grid<4>;
  using Grid5D = Grid<5>;  ///< e.g. eikonal S^2 vs. (pt(ob1), pt(ob2), acop(ob1,ob2), xi(ob1), xi(ob2))
}  // namespace cepgen::survival
using cepgen::survival::Grid1D;
using cepgen::survival::Grid2D;
using cepgen::survival::Grid3D;
using cepgen::survival::Grid4D;
using cepgen::survival::Grid5D;
REGISTER_SURVIVAL_FACTOR("grid1D", Grid1D);
REGISTER_SURVIVAL_FACTOR("grid2D", Grid2D);
REGISTER_SURVIVAL_FACTOR("grid3D", Grid3D);
REGISTER_SURVIVAL_FACTOR("grid4D", Grid4D);
REGISTER_SURVIVAL_FACTOR("grid5D", Grid5D);
