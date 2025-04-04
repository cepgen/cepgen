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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <limits>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/GridHandler.h"

//#define GRID_HANDLER_DEBUG 1

using namespace cepgen;

template <size_t D, size_t N>
GridHandler<D, N>::GridHandler(const GridType& grid_type) : grid_type_(grid_type) {
  for (size_t i = 0; i < D; ++i)
    accelerators_.emplace_back(gsl_interp_accel_alloc(), gsl_interp_accel_free);
}

template <size_t D, size_t N>
typename GridHandler<D, N>::values_t GridHandler<D, N>::eval(const coord_t& in_coords) const {
  if (!initialised_)
    throw CG_FATAL("GridHandler") << "Grid extrapolator called but not initialised!";

  values_t out{};
  coord_t coord = in_coords;
  switch (grid_type_) {
    case GridType::logarithmic: {
      std::transform(coord.begin(), coord.end(), coord.begin(), [](const auto& c) { return log10(c); });
    } break;
    case GridType::square: {
      std::transform(coord.begin(), coord.end(), coord.begin(), [](const auto& c) { return c * c; });
    } break;
    default:
      break;
  }
  switch (D) {  // dimension of the vector space coordinate to evaluate
    case 1: {
      for (size_t i = 0; i < N; ++i) {
        if (const auto res =
                gsl_spline_eval_e(splines_1d_.at(i).get(), coord.at(0), accelerators_.at(0).get(), &out[i]);
            res != GSL_SUCCESS) {
          out[i] = 0.;
          CG_WARNING("GridHandler") << "Failed to evaluate the value (N=" << i << ") "
                                    << "for x = " << in_coords.at(0) << " in grid with boundaries " << boundaries()
                                    << ". GSL error: " << gsl_strerror(res);
        }
      }
    } break;
    case 2: {
#ifdef GSL_VERSION_ABOVE_2_1
      const double x = coord.at(0), y = coord.at(1);
      for (size_t i = 0; i < N; ++i) {
        if (const auto res = gsl_spline2d_eval_e(
                splines_2d_.at(i).get(), x, y, accelerators_.at(0).get(), accelerators_.at(1).get(), &out[i]);
            res != GSL_SUCCESS) {
          out[i] = 0.;
          CG_WARNING("GridHandler") << "Failed to evaluate the value (N=" << i << ") "
                                    << "for x = " << x << " / y = " << y << " in grid with boundaries " << boundaries()
                                    << ". GSL error: " << gsl_strerror(res);
        }
      }
#else
      //--- retrieve the indices of the bin in the set
      coord_t before(D), after(D);
      findIndices(coord, before, after);
      //--- find boundaries values
      const gridpoint_t &ext_11 = values_raw_.at({before[0], before[1]}),
                        &ext_12 = values_raw_.at({before[0], after[1]}),
                        &ext_21 = values_raw_.at({after[0], before[1]}), &ext_22 = values_raw_.at({after[0], after[1]});
      //--- now that we have the boundaries, we may interpolate
      coord_t c_d(D);
      for (size_t i = 0; i < D; ++i)
        c_d[i] = (after[i] != before[i]) ? (coord.at(i) - before[i]) / (after[i] - before[i]) : 0.;
      const gridpoint_t ext_1 = ext_11 * (1. - c_d[0]) + ext_21 * c_d[0];
      const gridpoint_t ext_2 = ext_12 * (1. - c_d[0]) + ext_22 * c_d[0];
      out = ext_1 * (1. - c_d[1]) + ext_2 * c_d[1];
#endif
    } break;
    case 3: {
      // retrieve the indices of the bin in the set
      coord_t before(D), after(D);
      findIndices(coord, before, after);
      // find boundaries values
      const grid_point_t &ext_111 = values_raw_.at({before[0], before[1], before[2]}),
                         &ext_112 = values_raw_.at({before[0], before[1], after[2]}),
                         &ext_121 = values_raw_.at({before[0], after[1], before[2]}),
                         &ext_122 = values_raw_.at({before[0], after[1], after[2]}),
                         &ext_211 = values_raw_.at({after[0], before[1], before[2]}),
                         &ext_212 = values_raw_.at({after[0], before[1], after[2]}),
                         &ext_221 = values_raw_.at({after[0], after[1], before[2]}),
                         &ext_222 = values_raw_.at({after[0], after[1], after[2]});
      // now that we have the boundaries, we may interpolate
      coord_t c_d(D);
      for (size_t i = 0; i < D; ++i)
        c_d[i] = (after[i] != before[i]) ? (coord.at(i) - before[i]) / (after[i] - before[i]) : 0.;
      const grid_point_t ext_11 = ext_111 * (1. - c_d[0]) + ext_211 * c_d[0];
      const grid_point_t ext_12 = ext_112 * (1. - c_d[0]) + ext_212 * c_d[0];
      const grid_point_t ext_21 = ext_121 * (1. - c_d[0]) + ext_221 * c_d[0];
      const grid_point_t ext_22 = ext_122 * (1. - c_d[0]) + ext_222 * c_d[0];
      const grid_point_t ext_1 = ext_11 * (1. - c_d[1]) + ext_21 * c_d[1];
      const grid_point_t ext_2 = ext_12 * (1. - c_d[1]) + ext_22 * c_d[1];
      out = ext_1 * (1. - c_d[2]) + ext_2 * c_d[2];
    } break;
    default:
      throw CG_FATAL("GridHandler") << "Unsupported number of dimensions: " << N << ".\n\t"
                                    << "Please contact the developers to add such a new feature.";
  }
  return out;
}

template <size_t D, size_t N>
void GridHandler<D, N>::insert(const coord_t& coord, values_t value) {
  auto modified_coordinate = coord;
  if (grid_type_ != GridType::linear)
    for (auto& c : modified_coordinate)
      switch (grid_type_) {
        case GridType::logarithmic:
          c = std::log10(c);
          break;
        case GridType::square:
          c *= c;
          break;
        default:
          break;
      }
  if (values_raw_.count(modified_coordinate) > 0)
    CG_WARNING("GridHandler") << "Duplicate coordinate detected for x=" << coord << ".";
  values_raw_[modified_coordinate] = value;
  initialised_ = false;
}

template <size_t D, size_t N>
void GridHandler<D, N>::initialise() {
  if (values_raw_.empty())
    throw CG_ERROR("GridHandler") << "Empty grid.";
  gsl_set_error_handler_off();
  //--- start by building grid coordinates from raw values
  for (auto& coordinate : coordinates_)
    coordinate.clear();
  for (const auto& val : values_raw_) {
    unsigned short i = 0;
    for (const auto& c : val.first) {
      if (std::find(coordinates_.at(i).begin(), coordinates_.at(i).end(), c) == coordinates_.at(i).end())
        coordinates_.at(i).emplace_back(c);
      ++i;
    }
  }
  for (auto& c : coordinates_)
    std::sort(c.begin(), c.end());
#ifdef GRID_HANDLER_DEBUG
  CG_DEBUG("GridHandler").log([&](auto& dbg) {
    dbg << "Grid dump:";
    // debugging of the grid coordinates
    unsigned short i = 0;
    for (const auto& cs : coords_) {
      dbg << "\n>> coordinate " << (i++) << " has " << utils::s("member", cs.size(), true) << ":";
      unsigned short j = 0;
      for (const auto& val : cs)
        dbg << (j++ % 20 == 0 ? "\n  " : " ") << val;
    }
  });
#endif
  //--- particularise by dimension
  switch (D) {
    case 1: {  //--- x |-> (f1,...)
      const gsl_interp_type* type = gsl_interp_cspline;
      //const gsl_interp_type* type = gsl_interp_steffen;
#ifdef GSL_VERSION_ABOVE_2_1
      const unsigned short min_size = gsl_interp_type_min_size(type);
#else
      const unsigned short min_size = type->min_size;
#endif
      if (min_size >= values_raw_.size())
        throw CG_FATAL("GridHandler") << "Not enough points for \"" << type->name << "\" type of interpolation.\n\t"
                                      << "Minimum required: " << min_size << ", got " << values_raw_.size() << "!";
      for (size_t i = 0; i < N; ++i) {
        values_[i].reset(new double[values_raw_.size()]);
        splines_1d_.emplace_back(gsl_spline_alloc(type, values_raw_.size()), gsl_spline_free);
      }
      // transform a map(coord -> values) to
      // - 1 vector(coordinates)
      // - N vectors(values)
      std::vector<double> x_vec;
      size_t i = 0;
      for (const auto& vals : values_raw_) {
        x_vec.emplace_back(vals.first.at(0));
        unsigned short j = 0;
        for (const auto& val : vals.second)
          values_[j++].get()[i] = val;
        ++i;
      }
      // initialise spline interpolation objects (one for each value)
      for (size_t j = 0; j < splines_1d_.size(); ++j)
        if (const auto res = gsl_spline_init(splines_1d_.at(j).get(), &x_vec[0], values_[j].get(), values_raw_.size());
            res != GSL_SUCCESS)
          CG_WARNING("GridHandler:initialisation")
              << "Failed to initialise spline for dimension " << j << ". GSL error: " << gsl_strerror(res);
    } break;
    case 2: {  //--- (x,y) |-> (f1,...)
#ifdef GSL_VERSION_ABOVE_2_1
      if (values_raw_.size() < gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
        CG_WARNING("GridHandler") << "The grid size is too small (" << values_raw_.size() << " < "
                                  << gsl_interp2d_type_min_size(gsl_interp2d_bicubic)
                                  << ") for bicubic interpolation. Switching to a bi-linear interpolation mode.";
      //const gsl_interp2d_type* type =
      //    (values_raw_.size() > gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ? gsl_interp2d_bicubic
      //                                                                           : gsl_interp2d_bilinear);
      const gsl_interp2d_type* type = gsl_interp2d_bilinear;
      splines_2d_.clear();
      for (size_t i = 0; i < N; ++i) {
        values_[i].reset(new double[coordinates_.at(0).size() * coordinates_.at(1).size()]);
        splines_2d_.emplace_back(gsl_spline2d_alloc(type, coordinates_.at(0).size(), coordinates_.at(1).size()),
                                 gsl_spline2d_free);
      }

      // second loop over all points to populate the grid
      for (const auto& [coordinate, value] : values_raw_) {
        const auto &coord_x = coordinate.at(0), &coord_y = coordinate.at(1);
        // retrieve the index of the bin in the set
        const auto id_x =
                       std::distance(coordinates_.at(0).begin(),
                                     std::lower_bound(coordinates_.at(0).begin(), coordinates_.at(0).end(), coord_x)),
                   id_y =
                       std::distance(coordinates_.at(1).begin(),
                                     std::lower_bound(coordinates_.at(1).begin(), coordinates_.at(1).end(), coord_y));
        for (size_t i = 0; i < splines_2d_.size(); ++i)
          if (const auto res = gsl_spline2d_set(splines_2d_.at(i).get(), values_[i].get(), id_x, id_y, value[i]);
              res != GSL_SUCCESS)
            CG_WARNING("GridHandler:initialisation")
                << "Failed to set value " << i << " for spline. GSL error: " << gsl_strerror(res);
      }

      // initialise spline interpolation objects (one for each value)
      const coord_t &x_vec = coordinates_.at(0), &y_vec = coordinates_.at(1);
      for (size_t i = 0; i < splines_2d_.size(); ++i)
        gsl_spline2d_init(splines_2d_.at(i).get(), &x_vec[0], &y_vec[0], values_[i].get(), x_vec.size(), y_vec.size());
#else
      CG_WARNING("GridHandler") << "GSL version ≥ 2.1 is required for spline bilinear interpolation.\n\t"
                                << "Version " << GSL_VERSION << " is installed on this system!\n\t"
                                << "Will use a simple bilinear approximation instead.";
#endif
    } break;
    default:
      break;
  }
  initialised_ = true;
  CG_DEBUG("GridHandler").log([&](auto& log) {
    log << "Grid evaluator initialised with boundaries: " << boundaries() << ".";
#ifdef GRID_HANDLER_DEBUG
    log << "\n"
        << "Values handled:\n"
        << values_raw_;
#endif
  });
}

template <size_t D, size_t N>
std::array<Limits, D> GridHandler<D, N>::boundaries() const {
  std::array<Limits, D> out;
  const auto min_val = min(), max_val = max();
  for (size_t i = 0; i < D; ++i)
    out[i] = {min_val[i], max_val[i]};
  return out;
}

template <size_t D, size_t N>
std::array<double, D> GridHandler<D, N>::min() const {
  std::array<double, D> out{};
  size_t i = 0;
  for (const auto& coordinate : coordinates_) {  // loop over all dimensions
    const auto& min = std::min_element(coordinate.begin(), coordinate.end());
    out[i++] = min != coordinate.end() ? *min : std::numeric_limits<double>::infinity();
  }
  return out;
}

template <size_t D, size_t N>
std::array<double, D> GridHandler<D, N>::max() const {
  std::array<double, D> out{};
  size_t i = 0;
  for (const auto& coordinate : coordinates_) {  // loop over all dimensions
    const auto& max = std::max_element(coordinate.begin(), coordinate.end());
    out[i++] = max != coordinate.end() ? *max : std::numeric_limits<double>::infinity();
  }
  return out;
}

template <size_t D, size_t N>
void GridHandler<D, N>::findIndices(const coord_t& coord, coord_t& min, coord_t& max) const {
  for (size_t i = 0; i < D; ++i) {
    if (const auto& coordinate_i = coordinates_.at(i);  // extract all coordinates registered for this dimension
        coord.at(i) < *coordinate_i.begin()) {          // under the range
      CG_DEBUG_LOOP("GridHandler:indices") << "Coordinate " << i << " in underflow range "
                                           << "(" << coord.at(i) << " < " << *coordinate_i.begin() << ").";
      min[i] = max[i] = *coordinate_i.begin();
    } else if (coord.at(i) > *coordinate_i.rbegin()) {  // over the range
      CG_DEBUG_LOOP("GridHandler:indices") << "Coordinate " << i << " in overflow range "
                                           << "(" << coord.at(i) << " > " << *coordinate_i.rbegin() << ").";
      min[i] = max[i] = *coordinate_i.rbegin();
    } else {  // in between two coordinates
      auto it_coord = std::lower_bound(coordinate_i.begin(), coordinate_i.end(), coord.at(i));
      max[i] = *it_coord;
      if (it_coord != coordinate_i.begin())
        --it_coord;
      min[i] = *it_coord;
    }
  }
}

//----------------------------------------------------------------------------
// grid manipulation utility
//----------------------------------------------------------------------------

template <size_t D, size_t N>
typename GridHandler<D, N>::grid_point_t GridHandler<D, N>::grid_point_t::operator*(double c) const {
  grid_point_t out = *this;
  std::transform(out.begin(), out.end(), out.begin(), [&c](const auto& a) { return a * c; });
  return out;
}

template <size_t D, size_t N>
typename GridHandler<D, N>::grid_point_t GridHandler<D, N>::grid_point_t::operator+(const grid_point_t& rhs) const {
  grid_point_t out = *this;
  for (size_t i = 0; i < out.size(); ++i)
    out[i] += rhs[i];
  return out;
}

namespace cepgen {  // template specialisation for the few cases handled
  template class GridHandler<1, 1>;
  template class GridHandler<1, 2>;
  template class GridHandler<2, 2>;
  template class GridHandler<3, 1>;
}  // namespace cepgen
