/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include <functional>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

GridParameters::GridParameters(size_t mbin, size_t ndim) : mbin_(mbin), inv_mbin_(1. / mbin_), num_dimensions_(ndim) {
  // build and populate the grid
  coord_t coord(ndim, 0);
  for (size_t i = 0; i < static_cast<size_t>(std::pow(mbin_, num_dimensions_)); ++i) {
    generateCoordinates(coord, i);
    coords_.emplace_back(coord);
    num_points_.emplace_back(0ul);
    f_max_.emplace_back(0.);
  }
}

void GridParameters::setValue(size_t coordinate, float value) {
  // update function local and global maxima if needed
  f_max_.at(coordinate) = std::max(f_max_.at(coordinate), value);
  f_max_global_ = std::max(f_max_global_, value);
}

void GridParameters::shoot(const Integrator* integrator, size_t coordinate, std::vector<double>& out) const {
  CG_ASSERT(integrator != nullptr);
  const auto& nv = coords_.at(coordinate);
  for (size_t i = 0; i < nv.size(); ++i)
    out[i] = (integrator->uniform() + nv.at(i)) * inv_mbin_;
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

void GridParameters::generateCoordinates(coord_t& coord, size_t i) const {
  size_t jj = i;
  for (size_t j = 0; j < num_dimensions_; ++j) {
    const size_t tmp = jj * inv_mbin_;
    coord[j] = jj - tmp * mbin_;
    jj = tmp;
  }
}

bool GridParameters::correct(size_t bin) {
  if (f_max2_ <= f_max_.at(bin))
    return true;
  f_max_old_ = f_max_.at(bin);
  f_max_diff_ = f_max2_ - f_max_old_;
  correction_ = (num_points_.at(bin) - 1.) * f_max_diff_ / f_max_global_;
  if (f_max2_ >= f_max_global_)
    correction_ *= f_max2_ / f_max_global_;
  setValue(bin, f_max2_);
  correction_ -= correction2_;
  correction2_ = 0.;
  f_max2_ = 0.;
  return false;
}

void GridParameters::rescale(size_t bin, float weight) {
  if (weight <= f_max_.at(bin))
    return;
  f_max2_ = std::max(f_max2_, weight);
  correction_ += 1.;
  correction2_ -= 1.;
}

void GridParameters::initCorrectionCycle(size_t bin, float weight) {
  f_max_old_ = f_max_.at(bin);
  f_max_diff_ = weight - f_max_old_;
  setValue(bin, weight);
  correction_ = (num_points_.at(bin) - 1) * f_max_diff_ / f_max_global_ - 1.;

  CG_DEBUG("GridParameters:initCorrectionCycle")
      << "Correction " << correction_ << " will be applied "
      << "for phase space bin " << bin << " (" << utils::s("point", num_points_.at(bin), true) << "). "
      << "Maxima ratio: " << (f_max_diff_ / f_max_global_) << ".";
}
