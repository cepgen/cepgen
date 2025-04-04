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

#ifndef CepGen_Physics_GluonGrid_h
#define CepGen_Physics_GluonGrid_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/GridHandler.h"

/// Kimber-Martin-Ryskin un-integrated gluon densities
namespace kmr {
  /// A KMR un-integrated gluon densities grid interpolator
  class GluonGrid : cepgen::GridHandler<3, 1>, public cepgen::SteeredObject<GluonGrid> {
  public:
    static GluonGrid& get(const cepgen::ParametersList& params = {});  ///< Retrieve the grid interpolator singleton
    GluonGrid(const GluonGrid&) = delete;
    void operator=(const GridHandler&) = delete;

    static cepgen::ParametersDescription description();

    const std::string& path() const { return grid_path_; }      ///< Retrieve the path to the interpolation grid values
    double operator()(double x, double kt2, double mu2) const;  ///< Compute the gluon flux

  private:
    static constexpr const char* DEFAULT_KMR_GRID_PATH = "gluon_mmht2014nlo_Watt.dat";
    explicit GluonGrid(const cepgen::ParametersList&);
    /// Location of the grid to be interpolated
    const std::string grid_path_;
  };
}  // namespace kmr

#endif
