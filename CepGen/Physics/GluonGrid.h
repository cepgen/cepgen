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

#ifndef CepGen_Physics_GluonGrid_h
#define CepGen_Physics_GluonGrid_h

#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/ParametersDescription.h"

#define DEFAULT_KMR_GRID_PATH "gluon_mmht2014nlo_Watt.dat"

/// Kimber-Martin-Ryskin unintegrated gluon densities
namespace kmr {
  /// A KMR unintegrated gluon densities grid interpolator
  class GluonGrid : private cepgen::GridHandler<3, 1> {
  public:
    /// Retrieve the grid interpolator (singleton)
    static GluonGrid& get(const cepgen::ParametersList& params = {});

    GluonGrid(const GluonGrid&) = delete;
    void operator=(const GridHandler&) = delete;

    static cepgen::ParametersDescription parametersDescription();

    /// Retrieve the path to the interpolation grid values
    const std::string& path() const { return grid_path_; }

    /// Compute the gluon flux
    double operator()(double x, double kt2, double mu2) const;

  private:
    explicit GluonGrid(const cepgen::ParametersList&);
    /// Location of the grid to be interpolated
    const std::string grid_path_;
  };
}  // namespace kmr

#endif
