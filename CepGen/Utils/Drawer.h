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

#ifndef CepGen_Utils_Drawer_h
#define CepGen_Utils_Drawer_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  namespace utils {
    class Drawable;
    class Graph1D;
    class Graph2D;
    class Hist1D;
    class Hist2D;
    /// A generic drawing utility
    class Drawer : public NamedModule<std::string> {
    public:
      /// Build a drawing utility
      explicit Drawer(const ParametersList& params);

      enum struct Mode : int16_t { none = 0, logx, logy, logz, nostack };
      friend Mode operator|(const Mode&, const Mode&);
      friend bool operator&(const Mode&, const Mode&);

      /// Draw a one-dimensional graph
      virtual const Drawer& draw(const Graph1D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a two-dimensional graph
      virtual const Drawer& draw(const Graph2D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a one-dimensional histogram
      virtual const Drawer& draw(const Hist1D&, const Mode& mode = Mode::none) const = 0;
      /// Draw a two-dimensional histogram
      virtual const Drawer& draw(const Hist2D&, const Mode& mode = Mode::none) const = 0;

      /// A collection of drawable objects
      typedef std::vector<const Drawable*> DrawableColl;
      /// Draw a collection of drawables
      virtual const Drawer& draw(const DrawableColl&,
                                 const std::string& name = "",
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
}  // namespace cepgen

#endif
