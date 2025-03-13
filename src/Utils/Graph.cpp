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
  for (const auto& kv1 : values_) {
    if (oth.values_.count(kv1.first) == 0)
      throw CG_ERROR("Graph1D:chi2") << "Failed to retrieve the value for coordinate=" << kv1.first.value << "!\n"
                                     << "Please ensure the two graphs have the same values definition.";
    const auto& val1 = kv1.second;
    const auto& val2 = oth.values_.at(kv1.first);
    double norm = std::pow(val1.uncertainty(), 2) + std::pow(val2.uncertainty(), 2);
    if (norm == 0.)
      norm = 1.;
    chi2 += std::pow(val1 - val2, 2) / norm;
  }
  return chi2;
}

std::set<double> Graph1D::xCoords() const {
  std::set<double> coords;
  for (const auto& val : values_)
    coords.insert(val.first.value);
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
  for (const auto& x_axis : values_)
    for (const auto& yaxis : x_axis.second)
      os << format("\n%6zu: (%5g, %5g) = %5g", np++, x_axis.first.value, yaxis.first.value, yaxis.second);
}

std::set<double> Graph2D::xCoords() const {
  std::set<double> coords;
  for (const auto& x_value : values_)
    coords.insert(x_value.first.value);
  return coords;
}

std::set<double> Graph2D::yCoords() const {
  std::set<double> coords;
  for (const auto& x_value : values_)
    for (const auto& y_value : x_value.second)
      coords.insert(y_value.first.value);
  return coords;
}

const Value& Graph2D::valueAt(double x_value, double y_value) const {
  for (const auto& xv : values_)
    if (xv.first.value == x_value)
      for (const auto& yv : xv.second)
        if (yv.first.value == y_value)
          return yv.second;
  throw CG_ERROR("Graph2D:valueAt") << "Failed to retrieve a point a the coordinate (x=" << x_value << ", y=" << y_value
                                    << ").";
}
