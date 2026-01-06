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

#ifndef CepGen_Integration_GridParameters_h
#define CepGen_Integration_GridParameters_h

#include <vector>

namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

namespace cepgen {
  /// A parameters placeholder for the grid integration helper
  class GridParameters {
  public:
    /// Build a generation grid for a ndim-dimensional phase space
    explicit GridParameters(size_t m_bin, size_t num_dimensions);

    using coord_t = std::vector<unsigned short>;  ///< Coordinates definition

    void dump() const;  ///< Dump the grid coordinates

    inline size_t size() const { return coordinates_.size(); }  ///< Grid multiplicity
    /// Number of times a phase space point has been randomly selected
    inline const coord_t& n(size_t coord) const { return coordinates_.at(coord); }
    inline float globalMax() const { return f_max_global_; }  ///< Global function maximum

    /// Maximal function value for a given grid coordinate
    inline float maxValue(size_t coord) const { return f_max_.at(coord); }
    inline double maxValueDiff() const { return f_max_diff_; }
    inline double maxHistValue() const { return f_max_old_; }

    void setValue(size_t, float);  ///< Set the function value for a given grid coordinate
    /// Shoot a phase space point for a grid coordinate
    void shoot(utils::RandomGenerator& random_generator, size_t coordinate, std::vector<double>& out) const;
    /// Number of points already shot for a given grid coordinate
    inline size_t numPoints(size_t coordinate) const { return num_points_.at(coordinate); }
    /// Specify a new trial has been attempted for bin
    inline void increment(size_t coordinate) { num_points_.at(coordinate)++; }

    inline bool prepared() const { return gen_prepared_; }                       ///< Has the grid been prepared?
    inline void setPrepared(bool prepared = true) { gen_prepared_ = prepared; }  ///< Mark the grid as prepared

    /// Correction to apply on the next phase space point generation
    inline float correctionValue() const { return correction_; }
    /// Set the correction to apply on the next phase space point generation
    inline void setCorrectionValue(float correction) { correction_ = correction; }
    bool correct(size_t);  ///< Apply the correction requested at the previous generation

    void rescale(size_t, float);
    void initCorrectionCycle(size_t, float);

  private:
    void generateCoordinates(coord_t&, size_t) const;

    const size_t mbin_;         ///< Integration grid size parameter
    const double inv_mbin_;     ///< Weight of each grid coordinate
    size_t num_dimensions_{0};  ///< Phase space multiplicity
    bool gen_prepared_{false};  ///< Has the grid been already prepared?
    float correction_{0.};      ///< Correction to apply on the next phase space point generation
    float correction2_{0.};
    std::vector<coord_t> coordinates_;  ///< Point coordinates in grid
    std::vector<size_t> num_points_;    ///< Number of functions values evaluated for this point
    std::vector<float> f_max_;          ///< Maximal value of the function at one given point
    float f_max_global_{0.};            ///< Maximal value of the function in the considered integration range
    float f_max2_{0.};
    float f_max_diff_{0.};
    float f_max_old_{0.};
  };
}  // namespace cepgen

#endif
