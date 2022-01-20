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

#include <bitset>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/GridHandler.h"

namespace cepgen {
  namespace utils {
    Drawer::Drawer(const ParametersList& params) : NamedModule(params) {}

    Drawer::Mode operator|(const Drawer::Mode& lhs, const Drawer::Mode& rhs) {
      std::bitset<7> mod1((int)lhs), mod2((int)rhs);
      return (Drawer::Mode)(mod1 | mod2).to_ulong();
    }

    bool operator&(const Drawer::Mode& lhs, const Drawer::Mode& rhs) {
      //return ((int)lhs > (int)rhs) - ((int)lhs < (int)rhs);
      return (int)lhs & (int)rhs;
    }

    template <size_t N>
    const Drawer& Drawer::draw(const GridHandler<1, N>& grid, const Mode& mode) const {
      // first prepare the array of plots to generate
      std::array<Graph1D, N> plots;
      for (const auto& val : grid.values()) {
        for (size_t i = 0; i < N; ++i)
          plots[i].addPoint(val.first[0], val.second[i]);
      }
      for (size_t i = 0; i < N; ++i)
        draw(plots[i], mode);
      return *this;
    }

    template <size_t N>
    const Drawer& Drawer::draw(const GridHandler<2, N>& grid, const Mode& mode) const {
      // first prepare the array of plots to generate
      std::array<Graph2D, N> plots;
      for (const auto& val : grid.values()) {
        for (size_t i = 0; i < N; ++i)
          plots[i].addPoint(val.first[0], val.first[1], val.second[i]);
      }
      for (size_t i = 0; i < N; ++i)
        draw(plots[i], mode);
      return *this;
    }
    template const Drawer& Drawer::draw(const GridHandler<2, 2>& grid, const Mode&) const;
  }  // namespace utils
}  // namespace cepgen
