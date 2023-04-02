/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include <boost/histogram.hpp>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

namespace bh = boost::histogram;

namespace cepgen {
  namespace utils {
    class BoostDrawer : public Drawer {
    public:
      explicit BoostDrawer(const ParametersList& params) : Drawer(params) {}

      static ParametersDescription description() {
        auto desc = Drawer::description();
        return desc;
      }

      const BoostDrawer& draw(const Graph1D&, const Mode&) const override;
      const BoostDrawer& draw(const Graph2D&, const Mode&) const override;
      const BoostDrawer& draw(const Hist1D&, const Mode&) const override;
      const BoostDrawer& draw(const Hist2D&, const Mode&) const override;

      const BoostDrawer& draw(const DrawableColl&,
                              const std::string& name = "",
                              const std::string& title = "",
                              const Mode& = Mode::none) const override;

    private:
    };

    const BoostDrawer& BoostDrawer::draw(const Graph1D&, const Mode&) const {
      CG_WARNING("BoostDrawer") << "Not yet implemented.";
      return *this;
    }

    const BoostDrawer& BoostDrawer::draw(const Graph2D&, const Mode&) const {
      CG_WARNING("BoostDrawer") << "Not yet implemented.";
      return *this;
    }

    const BoostDrawer& BoostDrawer::draw(const Hist1D& hist, const Mode&) const {
      std::vector<double> ax_bins{hist.binRange(0).min()};
      for (size_t i = 0; i < hist.nbins(); ++i)
        ax_bins.emplace_back(hist.binRange(i).max());
      auto ax = bh::axis::variable<double>(ax_bins.begin(), ax_bins.end(), hist.xAxis().label());
      auto h = bh::make_histogram(std::move(ax));
      for (size_t i = 0; i < hist.nbins(); ++i)
        h[i] = hist.value(i);
      return *this;
    }

    const BoostDrawer& BoostDrawer::draw(const Hist2D& hist, const Mode&) const {
      std::vector<double> ax_bins{hist.binRangeX(0).min()}, ay_bins{hist.binRangeY(0).min()};
      for (size_t i = 0; i < hist.nbinsX(); ++i)
        ax_bins.emplace_back(hist.binRangeX(i).max());
      auto ax = bh::axis::variable<double>(ax_bins.begin(), ax_bins.end(), hist.xAxis().label());
      for (size_t i = 0; i < hist.nbinsY(); ++i)
        ay_bins.emplace_back(hist.binRangeY(i).max());
      auto ay = bh::axis::variable<double>(ay_bins.begin(), ay_bins.end(), hist.yAxis().label());
      auto h = bh::make_histogram(std::move(ax), std::move(ay));
      for (size_t i = 0; i < hist.nbinsX(); ++i)
        for (size_t j = 0; j < hist.nbinsY(); ++j)
          h[(i, j)] = hist.value(i, j);
      return *this;
    }

    const BoostDrawer& BoostDrawer::draw(const DrawableColl&,
                                         const std::string&,
                                         const std::string&,
                                         const Mode&) const {
      CG_WARNING("BoostDrawer") << "Not yet implemented.";
      return *this;
    }
  }  // namespace utils
}  // namespace cepgen
using cepgen::utils::BoostDrawer;
REGISTER_DRAWER("boost", BoostDrawer);
