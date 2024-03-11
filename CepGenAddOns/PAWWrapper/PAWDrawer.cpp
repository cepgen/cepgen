/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2022  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGenAddOns/PAWWrapper/PAWCommons.h"

extern "C" {
#include <cfortran/hbook.h>
}

namespace cepgen {
  namespace utils {
    /**
     * Handler for the storage of events in a PAW/HBOOK format
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Apr 2022
     */
    class PAWDrawer : public Drawer {
    public:
      /// Class constructor
      explicit PAWDrawer(const ParametersList& params) : Drawer(params) { HLIMIT(PAWC_SIZE); }

      static ParametersDescription description();

      const PAWDrawer& draw(const Graph1D&, const Mode&) const override;
      const PAWDrawer& draw(const Graph2D&, const Mode&) const override;
      const PAWDrawer& draw(const Hist1D&, const Mode&) const override;
      const PAWDrawer& draw(const Hist2D&, const Mode&) const override;

      const PAWDrawer& draw(const DrawableColl&,
                            const std::string& name = "",
                            const std::string& title = "",
                            const Mode& mode = Mode::none) const override;

    private:
    };

    ParametersDescription PAWDrawer::description() {
      auto desc = Drawer::description();
      return desc;
    }

    const PAWDrawer& PAWDrawer::draw(const Graph1D&, const Mode&) const {
      CG_WARNING("PAWDrawer:draw") << "Not yet implemented.";
      return *this;
    }

    const PAWDrawer& PAWDrawer::draw(const Graph2D&, const Mode&) const {
      CG_WARNING("PAWDrawer:draw") << "Not yet implemented.";
      return *this;
    }

    const PAWDrawer& PAWDrawer::draw(const Hist1D& hist, const Mode&) const {
      int ihist = 1;
      CG_LOG << "haha=" << ((char*)hist.name().data());
      HBOOK1(ihist, (char*)hist.name().data(), hist.nbins(), hist.range().min(), hist.range().max(), 0.);
      CG_LOG << "haha=" << ihist;
      for (size_t i = 0; i < hist.nbins(); ++i)
        HFILL(ihist, hist.binRange(i).x(0.5), 0., hist.value(i));
      HPRINT(ihist);
      return *this;
    }

    const PAWDrawer& PAWDrawer::draw(const Hist2D&, const Mode&) const {
      CG_WARNING("PAWDrawer:draw") << "Not yet implemented.";
      return *this;
    }

    const PAWDrawer& PAWDrawer::draw(const DrawableColl&,
                                     const std::string& name,
                                     const std::string& title,
                                     const Mode& mode) const {
      CG_WARNING("PAWDrawer:draw") << "Not yet implemented.";
      return *this;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("paw", PAWDrawer)
