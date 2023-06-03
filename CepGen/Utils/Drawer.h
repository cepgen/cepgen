/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#ifndef CepGen_Utils_Drawer_h
#define CepGen_Utils_Drawer_h

#include <cstdint>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  namespace utils {
    class Graph1D;
    class Graph2D;
    class Hist1D;
    class Hist2D;
    class Drawable;
    /// A collection of drawable objects
    typedef std::vector<const Drawable*> DrawableColl;
    /// A generic drawing utility
    class Drawer : public NamedModule<std::string> {
    public:
      /// Build a drawing utility
      explicit Drawer(const ParametersList& params);

      class Mode {
      public:
        enum value_t : uint16_t {
          none = 0,
          logx = 1 << 0,
          logy = 1 << 1,
          logz = 1 << 2,
          nostack = 1 << 3,
          grid = 1 << 4,
          col = 1 << 5,
          cont = 1 << 6
        };
        Mode() : value_(none) {}
        Mode(int val) : value_((value_t)val) {}
        Mode(const value_t& val) : value_(val) {}

        friend std::ostream& operator<<(std::ostream&, const Mode&);
        friend bool operator&(const Mode&, const Mode::value_t&);

        const value_t& value() const { return value_; }

      private:
        value_t value_;
      };

      /// Draw a one-dimensional graph
      virtual const Drawer& draw(const Graph1D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a two-dimensional graph
      virtual const Drawer& draw(const Graph2D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a one-dimensional histogram
      virtual const Drawer& draw(const Hist1D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a two-dimensional histogram
      virtual const Drawer& draw(const Hist2D&, const Mode& mode = Mode::none) const = 0;

      /// Draw a collection of drawables
      virtual const Drawer& draw(const DrawableColl&,
                                 const std::string& name = "",
                                 const std::string& title = "",
                                 const Mode& mode = Mode::none) const = 0;

      /// Output operator (when necessary)
      virtual std::ostream& operator<<(std::ostream& os) const { return os; }

    protected:
      friend class Drawable;
      friend class Graph1D;
      friend class Graph2D;
      friend class Hist1D;
      friend class Hist2D;
    };
  }  // namespace utils
  utils::Drawer::Mode operator|(const utils::Drawer::Mode&, const utils::Drawer::Mode::value_t&);
  utils::Drawer::Mode operator|(const utils::Drawer::Mode::value_t&, const utils::Drawer::Mode::value_t&);
}  // namespace cepgen
cepgen::utils::Drawer::Mode& operator|=(cepgen::utils::Drawer::Mode&, const cepgen::utils::Drawer::Mode::value_t&);

#endif

