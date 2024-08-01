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

#ifndef CepGen_Utils_Drawable_h
#define CepGen_Utils_Drawable_h

#include <map>
#include <string>

#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/Value.h"

namespace cepgen::utils {
  /// A generic object which can be drawn in the standard output
  class Drawable {
  public:
    explicit Drawable(const std::string& name = "", const std::string& title = "") : name_(name), title_(title) {}
    virtual ~Drawable() {}

    inline const std::string& name() const { return name_; }        ///< Drawable name
    inline void setName(const std::string& name) { name_ = name; }  ///< Set the drawable name

    inline const std::string& title() const { return title_; }          ///< Drawable name
    inline void setTitle(const std::string& title) { title_ = title; }  ///< Set the drawable title

    /// Metadata for an axis
    class AxisInfo {
    public:
      /// Set the axis title
      inline AxisInfo& setLabel(const std::string& label) {
        label_ = label;
        return *this;
      }
      inline const std::string& label() const { return label_; }  ///< Axis title
      /// Set the minimum range
      inline AxisInfo& setMinimum(double min) {
        lim_.min() = min;
        return *this;
      }
      /// Set the maximum range
      inline AxisInfo& setMaximum(double max) {
        lim_.max() = max;
        return *this;
      }
      /// Set the full axis range
      inline AxisInfo& setRange(const Limits& lim) {
        lim_ = lim;
        return *this;
      }
      inline const Limits& range() const { return lim_; }  ///< Axis range

    private:
      std::string label_;  ///< axis title
      Limits lim_;         ///< axis limits
    };

    AxisInfo& xAxis() { return xaxis_; }
    const AxisInfo& xAxis() const { return xaxis_; }
    AxisInfo& yAxis() { return yaxis_; }
    const AxisInfo& yAxis() const { return yaxis_; }
    AxisInfo& zAxis() { return zaxis_; }
    const AxisInfo& zAxis() const { return zaxis_; }

    /// Generic bin coordinate and its human-readable label
    struct coord_t {
      /// Sorting helper for axis coordinates
      bool operator<(const coord_t& oth) const { return value < oth.value; }
      double value{0.};        ///< Bin central value
      double value_unc = 0.;   ///< Bin uncertainty
      std::string label = "";  ///< Human-readable description of the bin
    };
    /// Metadata for an axis (coordinates and bins value)
    typedef std::map<coord_t, Value> axis_t;
    /// Comparator of an axis by the values it holds
    struct CompareAxisByValue {
      bool operator()(const std::pair<coord_t, Value>& lhs, const std::pair<coord_t, Value>& rhs) {
        return lhs.second < rhs.second;
      }
    };
    /// Metadata for a two-dimensional axis definition (coordinates and bins values)
    typedef std::map<coord_t, axis_t> dualaxis_t;

    inline virtual bool isHist1D() const { return false; }   ///< Is this drawable a one-dimensional histogram?
    inline virtual bool isHist2D() const { return false; }   ///< Is this drawable a two-dimensional histogram?
    inline virtual bool isGraph1D() const { return false; }  ///< Is this drawable a one-dimensional graph?
    inline virtual bool isGraph2D() const { return false; }  ///< Is this drawable a two-dimensional graph?

  protected:
    std::string name_;   ///< Computer-readable name
    std::string title_;  ///< Human-readable title
    AxisInfo xaxis_;     ///< x-axis metadata
    AxisInfo yaxis_;     ///< y-axis metadata
    AxisInfo zaxis_;     ///< z-axis metadata
  };
}  // namespace cepgen::utils

#endif
