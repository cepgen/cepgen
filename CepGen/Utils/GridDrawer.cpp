/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/GridDrawer.h"
#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/String.h"

namespace cepgen::utils {
  GridDrawer::GridDrawer(const ParametersList& params)
      : SteeredObject(params), drawer_(DrawerFactory::get().build(params_)) {}

  template <size_t N>
  void GridDrawer::draw(const GridHandler<1, N>& grid, const Drawer::Mode& mode) {
    auto gd = GridDrawer(ParametersList());
    // first prepare the array of plots to generate
    std::array<Graph1D, N> plots;
    for (size_t i = 0; i < N; ++i) {
      plots[i].xAxis().setLabel("x");
      plots[i].yAxis().setLabel(utils::format("var%d", i));
    }
    for (const auto& val : grid.values()) {
      for (size_t i = 0; i < N; ++i)
        plots[i].addPoint(val.first[0], val.second[i]);
    }
    for (size_t i = 0; i < N; ++i)
      gd.drawer_->draw(plots[i], mode);
  }

  template <size_t N>
  void GridDrawer::draw(const GridHandler<2, N>& grid, const Drawer::Mode& mode) {
    auto gd = GridDrawer(ParametersList());
    // first prepare the array of plots to generate
    std::array<Graph2D, N> plots;
    for (size_t i = 0; i < N; ++i) {
      plots[i].xAxis().setLabel("x0");
      plots[i].yAxis().setLabel("x1");
      plots[i].zAxis().setLabel(utils::format("var%d", i));
    }
    for (const auto& val : grid.values()) {
      for (size_t i = 0; i < N; ++i)
        plots[i].addPoint(val.first[0], val.first[1], val.second[i]);
    }
    for (size_t i = 0; i < N; ++i)
      gd.drawer_->draw(plots[i], mode);
  }

  template <size_t N>
  void GridDrawer::draw(const GridHandler<3, N>& grid, const Drawer::Mode& mode) {
    auto gd = GridDrawer(ParametersList());
    // first prepare the array of plots to generate
    std::array<std::array<Graph2D, 3>, N> plots;
    for (size_t i = 0; i < N; ++i) {
      plots[i][0].xAxis().setLabel("x0");
      plots[i][0].yAxis().setLabel("x1");
      plots[i][0].zAxis().setLabel(utils::format("var%d", i));
      plots[i][1].xAxis().setLabel("x0");
      plots[i][1].yAxis().setLabel("x2");
      plots[i][1].zAxis().setLabel(utils::format("var%d", i));
      plots[i][2].xAxis().setLabel("x1");
      plots[i][2].yAxis().setLabel("x2");
      plots[i][2].zAxis().setLabel(utils::format("var%d", i));
    }
    for (const auto& val : grid.values())
      for (size_t i = 0; i < N; ++i) {
        plots[i][0].addPoint(val.first[0], val.first[1], val.second[i]);
        plots[i][1].addPoint(val.first[0], val.first[2], val.second[i]);
        plots[i][2].addPoint(val.first[1], val.first[2], val.second[i]);
      }
    DrawableColl coll;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < 3; ++j)
        coll.emplace_back(&plots[i][j]);
    gd.drawer_->draw(coll, "", "", mode);
  }

  ParametersDescription GridDrawer::description() {
    auto desc = ParametersDescription();
    desc.setName("root");
    return desc;
  }

  template void GridDrawer::draw(const GridHandler<1, 1>&, const Drawer::Mode&);
  template void GridDrawer::draw(const GridHandler<1, 2>&, const Drawer::Mode&);
  template void GridDrawer::draw(const GridHandler<2, 2>&, const Drawer::Mode&);
  template void GridDrawer::draw(const GridHandler<3, 1>&, const Drawer::Mode&);
}  // namespace cepgen::utils
