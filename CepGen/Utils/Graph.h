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

#ifndef CepGen_Utils_Graph_h
#define CepGen_Utils_Graph_h

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace utils {
    /// A one-dimensional graph object
    class Graph1D : public Drawable {
    public:
      explicit Graph1D(const std::string& name = "", const std::string& title = "");

      /// Add one value to the graph
      Graph1D& addPoint(double x, double y);
      /// Add one value and its uncertainties to the graph
      Graph1D& addPoint(double x, double y, double ex, double ey);
      /// Retrieve all values in the graph
      const axis_t& points() const { return values_; }
      /// Minimum value registered in this graph
      double minimum() const;
      /// Maximum value registered in this graph
      double maximum() const;
      /// Compute the \f$\chi^{2}\f$ between this graph and another
      double chi2(const Graph1D&) const;

      /// List of horizontal axis coordinates
      std::set<double> xCoords() const;
      /// Retrieve the value of the graph at a given coordinate
      const value_t valueAt(double) const;

      bool isGraph1D() const override { return true; }

    private:
      axis_t values_;
    };

    /// A two-dimensional graph object
    class Graph2D : public Drawable {
    public:
      explicit Graph2D(const std::string& name = "", const std::string& title = "");

      /// Add one value to the graph
      Graph2D& addPoint(double x, double y, double z);
      /// Add one value and its uncertainties to the graph
      Graph2D& addPoint(double x, double y, double z, double ex, double ey, double ez);
      /// Retrieve all values in the graph
      const dualaxis_t& points() const { return values_; }
      /// List all values registered in the graph
      void dumpPoints(std::ostream&) const;

      /// List of horizontal axis coordinates
      std::set<double> xCoords() const;
      /// List of vertical axis coordinates
      std::set<double> yCoords() const;
      /// Retrieve the value of the graph at the given coordinates
      const value_t valueAt(double, double) const;

      bool isGraph2D() const override { return true; }

    private:
      dualaxis_t values_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
