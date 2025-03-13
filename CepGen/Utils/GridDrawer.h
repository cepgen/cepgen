/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#ifndef CepGen_Utils_GridDrawer_h
#define CepGen_Utils_GridDrawer_h

#include <memory>

#include "CepGen/Utils/Drawer.h"

namespace cepgen {
  template <size_t D, size_t N>
  class GridHandler;
}

namespace cepgen::utils {
  /// Utility object to draw a grid values mapping
  class GridDrawer final : public SteeredObject<GridDrawer> {
  public:
    static ParametersDescription description();

    /// Debugging drawing routing for single-dimensional grids
    template <size_t N>
    static void draw(const GridHandler<1, N>&, const Drawer::Mode& mode = Drawer::Mode::none);
    /// Debugging drawing routing for double-dimensional grids
    template <size_t N>
    static void draw(const GridHandler<2, N>&, const Drawer::Mode& mode = Drawer::Mode::none);
    /// Debugging drawing routing for triple-dimensional grids
    template <size_t N>
    static void draw(const GridHandler<3, N>&, const Drawer::Mode& mode = Drawer::Mode::none);

  private:
    explicit GridDrawer(const ParametersList&);
    const std::unique_ptr<Drawer> drawer_;
  };
}  // namespace cepgen::utils

#endif
