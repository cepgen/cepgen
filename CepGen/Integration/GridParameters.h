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

#ifndef CepGen_Integration_GridParameters_h
#define CepGen_Integration_GridParameters_h

#include <cstddef>
#include <vector>

namespace cepgen {
  class Integrator;
  /// A parameters placeholder for the grid integration helper
  class GridParameters {
  public:
    /// Build a generation grid for a ndim-dimensional phase space
    explicit GridParameters(size_t ndim);

    /// Coordinates definition
    typedef std::vector<unsigned short> coord_t;

    /// Dump the grid coordinates
    void dump() const;

    /// Grid multiplicity
    size_t size() const;
    /// Number of times a phase space point has been randomly selected
    const coord_t& n(size_t coord) const;
    /// Global function maximum
    float globalMax() const { return f_max_global_; }
    /// Maximal function value for a given grid coordinate
    float maxValue(size_t) const;
    /// Set the function value for a given grid coordinate
    void setValue(size_t, float);
    /// Shoot a phase space point for a grid coordinate
    void shoot(const Integrator* integ, size_t coord, std::vector<double>& out) const;
    /// Specify a new trial has been attempted for bin
    void increment(size_t coord);
    /// Number of points already shot for a given grid coordinate
    size_t numPoints(size_t coord) const;
    /// Has the grid been prepared
    bool prepared() const { return gen_prepared_; }
    /// Mark the grid as prepared
    void setPrepared(bool prepared = true) { gen_prepared_ = prepared; }
    /// Correction to apply on the next phase space point generation
    float correctionValue() const { return correc_; }
    /// Set the correction to apply on the next phase space point generation
    void setCorrectionValue(float correc) { correc_ = correc; }
    /// Apply the correction requested at the previous generation
    bool correct(size_t);
    void rescale(size_t, float);
    void initCorrectionCycle(size_t, float);
    double maxValueDiff() const { return f_max_diff_; }
    double maxHistValue() const { return f_max_old_; }

    /// Integration grid size parameter
    static constexpr unsigned short M_BIN = 3;
    /// Weight of each grid coordinate
    static constexpr double INV_M_BIN = 1. / M_BIN;

  private:
    void generateCoordinates(coord_t&, size_t) const;
    /// Phase space multiplicity
    size_t ndim_{0};
    /// Has the grid been already prepared?
    bool gen_prepared_{false};
    /// Correction to apply on the next phase space point generation
    float correc_{0.};
    float correc2_{0.};
    /// Point coordinates in grid
    std::vector<coord_t> coords_;
    /// Number of functions values evaluated for this point
    std::vector<size_t> num_points_;
    /// Maximal value of the function at one given point
    std::vector<float> f_max_;
    /// Maximal value of the function in the considered integration range
    float f_max_global_{0.};
    float f_max2_{0.};
    float f_max_diff_{0.};
    float f_max_old_{0.};
  };
}  // namespace cepgen

#endif
