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

#include <cmath>  // pow

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"

namespace cepgen {
  GridParameters::GridParameters(size_t ndim) : ndim_(ndim) {
    //--- build and populate the grid
    coord_t coord(ndim, 0);
    for (size_t i = 0; i < (size_t)pow(M_BIN, ndim_); ++i) {
      size_t jj = i;
      for (size_t j = 0; j < ndim; ++j) {
        size_t tmp = jj * INV_M_BIN;
        coord[j] = jj - tmp * M_BIN;
        jj = tmp;
      }
      coords_.emplace_back(coord);
      num_points_.emplace_back(0ul);
      f_max_.emplace_back(0.);
    }
  }

  size_t GridParameters::size() const { return coords_.size(); }

  const GridParameters::coord_t& GridParameters::n(size_t coord) const { return coords_.at(coord); }

  void GridParameters::setValue(size_t coord, double val) {
    //--- update function local and global maxima if needed
    f_max_.at(coord) = std::max(f_max_.at(coord), val);
    f_max_global_ = std::max(f_max_global_, val);
  }

  double GridParameters::maxValue(size_t coord) const { return f_max_.at(coord); }

  size_t GridParameters::numPoints(size_t coord) const { return num_points_.at(coord); }

  void GridParameters::increment(size_t coord) { num_points_.at(coord)++; }

  void GridParameters::shoot(const Integrator* integr, size_t coord, std::vector<double>& out) const {
    const auto& nv = coords_.at(coord);
    for (size_t i = 0; i < nv.size(); ++i)
      out[i] = (integr->uniform() + nv.at(i)) * INV_M_BIN;
  }

  void GridParameters::dump() const {
    CG_INFO("GridParameters:dump").log([&](auto& info) {
      for (size_t i = 0; i < coords_.size(); ++i)
        info << "\nn[" << i << "]: "
             << "coord=" << coords_.at(i) << ", "
             << "num points: " << num_points_.at(i) << ", "
             << "max=" << f_max_.at(i) << ".";
    });
  }
}  // namespace cepgen
