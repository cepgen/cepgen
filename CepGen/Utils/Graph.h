/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Graph_h
#define CepGen_Utils_Graph_h

#include <set>

#include "CepGen/Utils/Drawable.h"

namespace cepgen::utils {
  /// A one-dimensional graph object
  class Graph1D final : public Drawable {
  public:
    explicit Graph1D(const std::string& name = "", const std::string& title = "");

    Graph1D& addPoint(double x, double y);                        ///< Add one value to the graph
    Graph1D& addPoint(double x, double y, double ex, double ey);  ///< Add one value and its uncertainties to the graph
    inline const axis_t& points() const { return values_; }       ///< Retrieve all values in the graph
    double minimum() const;                                       ///< Minimum value registered in this graph
    double maximum() const;                                       ///< Maximum value registered in this graph
    double chi2(const Graph1D&) const;  ///< Compute the \f$\chi^{2}\f$ between this graph and another

    std::set<double> xCoords() const;   ///< List of horizontal axis coordinates
    const Value valueAt(double) const;  ///< Retrieve the value of the graph at a given coordinate

    inline bool isGraph1D() const override { return true; }

  private:
    axis_t values_;
  };

  /// A two-dimensional graph object
  class Graph2D final : public Drawable {
  public:
    explicit Graph2D(const std::string& name = "", const std::string& title = "");

    Graph2D& addPoint(double x, double y, double z);  ///< Add one value to the graph
    /// Add one value and its uncertainties to the graph
    Graph2D& addPoint(double x, double y, double z, double ex, double ey, double ez);
    inline const dualaxis_t& points() const { return values_; }  ///< Retrieve all values in the graph
    void dumpPoints(std::ostream&) const;                        ///< List all values registered in the graph

    std::set<double> xCoords() const;           ///< List of horizontal axis coordinates
    std::set<double> yCoords() const;           ///< List of vertical axis coordinates
    const Value valueAt(double, double) const;  ///< Retrieve the value of the graph at the given coordinates

    inline bool isGraph2D() const override { return true; }

  private:
    dualaxis_t values_;
  };
}  // namespace cepgen::utils

#endif
