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

#ifndef CepGen_Utils_GridHandler_h
#define CepGen_Utils_GridHandler_h

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_version.h>
#if defined(GSL_MAJOR_VERSION) && (GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1))
#include <gsl/gsl_spline2d.h>
#define GSL_VERSION_ABOVE_2_1 1
#endif

#include <array>
#include <map>
#include <memory>

#include "CepGen/Utils/Limits.h"

namespace cepgen {
  /// Interpolation type for the grid coordinates
  enum struct GridType { linear, logarithmic, square };
  /// A generic class for \f$\mathbb{R}^D\mapsto\mathbb{R}^N\f$ grid interpolation
  /// \tparam D Number of variables in the grid (dimension)
  /// \tparam N Number of values handled per point
  template <size_t D, size_t N = 1>
  class GridHandler {
  public:
    explicit GridHandler(const GridType& grid_type);  ///< Build a grid interpolator from a grid type
    virtual ~GridHandler() = default;

    typedef std::vector<double> coord_t;     ///< Coordinates container
    typedef std::array<double, N> values_t;  ///< Value(s) at a given coordinate

    values_t eval(const coord_t& in_coords) const;  ///< Interpolate a point to a given coordinate

    void insert(const coord_t& coord, values_t value);                         ///< Insert a new value in the grid
    inline std::map<coord_t, values_t> values() const { return values_raw_; }  ///< List of values in the grid

    void initialise();                         ///< Initialise the grid and all useful interpolators/accelerators
    std::array<Limits, D> boundaries() const;  ///< Grid boundaries (collection of (min,max))
    std::array<double, D> min() const;         ///< Lowest bound of the grid coordinates
    std::array<double, D> max() const;         ///< Highest bound of the grid coordinates

  protected:
    const GridType grid_type_;                ///< Type of interpolation for the grid members
    std::map<coord_t, values_t> values_raw_;  ///< List of coordinates and associated value(s) in the grid
    /// Grid interpolation accelerator
    std::vector<std::unique_ptr<gsl_interp_accel, void (*)(gsl_interp_accel*)> > accel_;
    std::vector<std::unique_ptr<gsl_spline, void (*)(gsl_spline*)> > splines_1d_;  ///< Splines for linear interpolations
#ifdef GSL_VERSION_ABOVE_2_1
    /// Splines for bilinear interpolations
    std::vector<std::unique_ptr<gsl_spline2d, void (*)(gsl_spline2d*)> > splines_2d_;
#endif
    std::array<coord_t, D> coords_;                    ///< Coordinates building up the grid
    std::array<std::unique_ptr<double[]>, N> values_;  ///< Values for all points in the grid

  private:
    void findIndices(const coord_t& coord, coord_t& min, coord_t& max) const;  ///< Lower/upper indices for a coordinate
    /// A single value in grid coordinates
    struct grid_point_t : values_t {
      grid_point_t(const values_t& arr) : values_t(arr) {}
      grid_point_t operator*(double c) const;
      grid_point_t operator+(const grid_point_t& rhs) const;
    };
    bool init_{false};  ///< Has the extrapolator been initialised?
  };
}  // namespace cepgen

#endif
