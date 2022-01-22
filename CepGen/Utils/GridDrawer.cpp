/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

namespace cepgen {
  namespace utils {
    GridDrawer::GridDrawer(const ParametersList& params)
        : SteeredObject(params), drawer_(DrawerFactory::get().build(params_)) {}

    template <size_t N>
    void GridDrawer::draw(const GridHandler<1, N>& grid, const Drawer::Mode& mode) {
      auto gd = GridDrawer(ParametersList());
      // first prepare the array of plots to generate
      std::array<Graph1D, N> plots;
      for (size_t i = 0; i < N; ++i)
        plots[i].setXlabel("x").setYlabel(utils::format("var%d", i));
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
      for (size_t i = 0; i < N; ++i)
        plots[i].setXlabel("x0").setYlabel("x1").setZlabel(utils::format("var%d", i));
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
        plots[i][0].setXlabel("x0").setYlabel("x1").setZlabel(utils::format("var%d", i));
        plots[i][1].setXlabel("x0").setYlabel("x2").setZlabel(utils::format("var%d", i));
        plots[i][2].setXlabel("x1").setYlabel("x2").setZlabel(utils::format("var%d", i));
      }
      for (const auto& val : grid.values())
        for (size_t i = 0; i < N; ++i) {
          plots[i][0].addPoint(val.first[0], val.first[1], val.second[i]);
          plots[i][1].addPoint(val.first[0], val.first[2], val.second[i]);
          plots[i][2].addPoint(val.first[1], val.first[2], val.second[i]);
        }
      for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < 3; ++j)
          gd.drawer_->draw(plots[i][j], mode);
    }

    ParametersDescription GridDrawer::description() {
      auto desc = ParametersDescription();
      desc.add<std::string>(ParametersList::MODULE_NAME, "root");
      return desc;
    }

    template void GridDrawer::draw(const GridHandler<1, 1>&, const Drawer::Mode&);
    template void GridDrawer::draw(const GridHandler<1, 2>&, const Drawer::Mode&);
    template void GridDrawer::draw(const GridHandler<2, 2>&, const Drawer::Mode&);
    template void GridDrawer::draw(const GridHandler<3, 1>&, const Drawer::Mode&);
  }  // namespace utils
}  // namespace cepgen
