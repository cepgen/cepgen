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

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Graph1D::Graph1D(const std::string& name, const std::string& title) : Drawable(name, title) {}

    void Graph1D::addPoint(double x, double y) { values_[coord_t{x}] = value_t{y}; }

    double Graph1D::minimum() const {
      return std::min_element(values_.begin(), values_.end(), CompareAxisByValue())->second.value;
    }

    double Graph1D::maximum() const {
      return std::max_element(values_.begin(), values_.end(), CompareAxisByValue())->second.value;
    }

    Graph2D::Graph2D(const std::string& name, const std::string& title) : Drawable(name, title) {}

    void Graph2D::addPoint(double x, double y, double z) { values_[coord_t{x}][coord_t{y}] = value_t{z}; }

    void Graph2D::dumpPoints(std::ostream& os) const {
      os << "Points registered in the 2D graph:";
      size_t np = 0ul;
      for (const auto& xaxis : values_)
        for (const auto& yaxis : xaxis.second)
          os << utils::format(
              "\n%6zu: (%5g, %5g) = %5g", np++, xaxis.first.value, yaxis.first.value, yaxis.second.value);
    }
  }  // namespace utils
}  // namespace cepgen
