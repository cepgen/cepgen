/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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
#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::utils;

Graph1D::Graph1D(const std::string& name, const std::string& title) : Drawable(name, title) {}

Graph1D& Graph1D::addPoint(double x, double y) {
  values_[coord_t{x}] = Value{y};
  return *this;
}

Graph1D& Graph1D::addPoint(double x, double y, double ex, double ey) {
  values_[coord_t{x, ex}] = Value{y, ey};
  return *this;
}

double Graph1D::minimum() const {
  return std::min_element(values_.begin(), values_.end(), CompareAxisByValue())->second;
}

double Graph1D::maximum() const {
  return std::max_element(values_.begin(), values_.end(), CompareAxisByValue())->second;
}

double Graph1D::chi2(const Graph1D& oth) const {
  double chi2 = 0.;
  if (values_.size() != oth.values_.size())
    throw CG_ERROR("Graph1D:chi2") << "Graphs must have the same number of elements to compute chi^2!";
  for (const auto& [x_coord, x_value] : values_) {
    if (oth.values_.count(x_coord) == 0)
      throw CG_ERROR("Graph1D:chi2") << "Failed to retrieve the value for coordinate=" << x_coord.value << "!\n"
                                     << "Please ensure the two graphs have the same values definition.";
    const auto& val2 = oth.values_.at(x_coord);
    double norm = std::pow(x_value.uncertainty(), 2) + std::pow(val2.uncertainty(), 2);
    if (norm == 0.)
      norm = 1.;
    chi2 += std::pow(x_value - val2, 2) / norm;
  }
  return chi2;
}

std::set<double> Graph1D::xCoords() const {
  std::set<double> coords;
  for (const auto& [x_coord, x_value] : values_)
    coords.insert(x_coord.value);
  return coords;
}

const Value& Graph1D::valueAt(double val) const {
  auto it = std::find_if(values_.begin(), values_.end(), [&val](const auto& xv) { return xv.first.value == val; });
  if (it == values_.end())
    throw CG_ERROR("Graph1D:valueAt") << "Failed to retrieve a point a the coordinate x=" << val << ".";
  return it->second;
}

Graph2D::Graph2D(const std::string& name, const std::string& title) : Drawable(name, title) {}

Graph2D& Graph2D::addPoint(double x, double y, double z) {
  values_[coord_t{x}][coord_t{y}] = Value{z};
  return *this;
}

Graph2D& Graph2D::addPoint(double x, double y, double z, double ex, double ey, double ez) {
  values_[coord_t{x, ex}][coord_t{y, ey}] = Value{z, ez};
  return *this;
}

void Graph2D::dumpPoints(std::ostream& os) const {
  os << "Points registered in the 2D graph:";
  size_t np = 0ul;
  for (const auto& [x_coord, x_value] : values_)
    for (const auto& [y_coord, y_value] : x_value)
      os << format("\n%6zu: (%5g, %5g) = %5g", np++, x_coord.value, y_coord.value, y_value);
}

std::set<double> Graph2D::xCoords() const {
  std::set<double> coords;
  for (const auto& [x_coord, x_value] : values_)
    coords.insert(x_coord.value);
  return coords;
}

std::set<double> Graph2D::yCoords() const {
  std::set<double> coords;
  for (const auto& [x_coord, x_value] : values_)
    for (const auto& [y_coord, y_value] : x_value)
      coords.insert(y_coord.value);
  return coords;
}

const Value& Graph2D::valueAt(double x_index, double y_index) const {
  for (const auto& [x_coord, x_value] : values_)
    if (x_coord.value == x_index)
      for (const auto& [y_coord, y_value] : x_value)
        if (y_coord.value == y_index)
          return y_value;
  throw CG_ERROR("Graph2D:valueAt") << "Failed to retrieve a point a the coordinate (x=" << x_index << ", y=" << y_index
                                    << ").";
}
