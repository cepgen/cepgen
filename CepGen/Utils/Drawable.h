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

#ifndef CepGen_Utils_Drawable_h
#define CepGen_Utils_Drawable_h

#include <map>
#include <string>
#include <vector>

#include "CepGen/Utils/Limits.h"
#include "Drawer.h"

namespace cepgen {
  namespace utils {
    /// A generic object which can be drawn in the standard output
    class Drawable {
    public:
      explicit Drawable(const std::string& name = "", const std::string& title = "") : name_(name), title_(title) {}
      ///< Copy constructor
      Drawable(const Drawable& oth) : xlabel_(oth.xlabel_), ylabel_(oth.ylabel_) {}

      /// Drawable name
      const std::string& name() const { return name_; }
      /// Set the drawable name
      void setName(const std::string& name) { name_ = name; }

      /// Drawable name
      const std::string& title() const { return title_; }
      /// Set the drawable title
      void setTitle(const std::string& title) { title_ = title; }

      /// x-axis label
      const std::string& xLabel() const { return xlabel_; }
      /// Set the x-axis label
      Drawable& setXlabel(const std::string& lab) {
        xlabel_ = lab;
        return *this;
      }
      /// y-axis label
      const std::string& yLabel() const { return ylabel_; }
      /// Set the y-axis label
      Drawable& setYlabel(const std::string& lab) {
        ylabel_ = lab;
        return *this;
      }
      /// z-axis label
      const std::string& zLabel() const { return zlabel_; }
      /// Set the z-axis label
      Drawable& setZlabel(const std::string& lab) {
        zlabel_ = lab;
        return *this;
      }

      /// Generic bin coordinate and its human-readable label
      struct coord_t {
        /// Sorting helper for axis coordinates
        bool operator<(const coord_t& oth) const { return value < oth.value; }
        double value;            ///< Bin central value
        std::string label = "";  ///< Human-readable description of the bin
      };
      /// Helper view of a pair of bin value and its uncertainty
      struct value_t {
        /// Sorting helper for bin values
        bool operator<(const value_t& oth) const { return value < oth.value; }
        double value;           ///< Single bin content
        double value_unc = 0.;  ///< Uncertainty on bin content
      };
      /// Metadata for an axis (coordinates and bins value)
      typedef std::map<coord_t, value_t> axis_t;
      /// Comparator of an axis by the values it holds
      struct CompareAxisByValue {
        bool operator()(const std::pair<coord_t, value_t>& lhs, const std::pair<coord_t, value_t>& rhs) {
          return lhs.second.value < rhs.second.value;
        }
      };
      /// Metadata for a two-dimensional axis definition (coordinates and bins values)
      typedef std::map<coord_t, axis_t> dualaxis_t;

      virtual bool isHist1D() const { return false; }   ///< Is this drawable a one-dimensional histogram?
      virtual bool isHist2D() const { return false; }   ///< Is this drawable a two-dimensional histogram?
      virtual bool isGraph1D() const { return false; }  ///< Is this drawable a one-dimensional graph?
      virtual bool isGraph2D() const { return false; }  ///< Is this drawable a two-dimensional graph?

    protected:
      std::string name_;    ///< Computer-readable name
      std::string title_;   ///< Human-readable title
      std::string xlabel_;  ///< x-axis title
      std::string ylabel_;  ///< y-axis title
      std::string zlabel_;  ///< z-axis title
    };
  }  // namespace utils
}  // namespace cepgen

#endif
