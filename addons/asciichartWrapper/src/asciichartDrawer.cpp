/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <ascii/ascii.h>

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Piper.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

using namespace std::string_literals;

namespace cepgen::utils {
  /// asciichart drawable objects drawing utility
  class asciichartDrawer final : public Drawer {
  public:
    /// Class constructor
    explicit asciichartDrawer(const ParametersList& params)
        : Drawer(params),
          type_(steerAs<int, ascii::Asciichart::Type>("type")),
          height_(steer<int>("height")),
          show_legend_(steer<bool>("showLegend")) {
      if (height_ <= 1)
        throw CG_FATAL("asciichartDrawer") << "Invalid chart height specified: " << height_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Drawer::description();
      desc.setDescription("asciichart drawing utility");
      desc.addAs<int>("type", ascii::Asciichart::Type::LINE);
      desc.add("height", 8);
      desc.add("showLegend", true);
      return desc;
    }

    const asciichartDrawer& draw(const Graph1D& graph, const Mode&) const override {
      std::vector<double> serie;
      for (const auto& [xval, yval] : graph.points())
        serie.emplace_back(yval);
      ascii::Asciichart chart({{graph.title(), serie}});
      CG_LOG << chart.type(type_).height(height_).Plot();
      return *this;
    }
    const asciichartDrawer& draw(const Graph2D&, const Mode&) const override {
      CG_WARNING("asciichartDrawer:draw") << "Unsupported graphical element.";
      return *this;
    }
    const asciichartDrawer& draw(const Hist1D& hist, const Mode&) const override {
      std::vector<double> serie;
      for (const auto& value : hist.values())
        serie.emplace_back(value);
      ascii::Asciichart chart({{hist.title(), serie}});
      CG_LOG << chart.type(type_).height(height_).Plot();
      return *this;
    }
    const asciichartDrawer& draw(const Hist2D&, const Mode&) const override {
      CG_WARNING("asciichartDrawer:draw") << "Unsupported graphical element.";
      return *this;
    }

    const asciichartDrawer& draw(const DrawableColl& coll,
                                 const std::string&,
                                 const std::string&,
                                 const Mode&) const override {
      std::unordered_map<std::string, std::vector<double> > series;
      for (const auto& obj : coll) {
        if (obj->isHist1D())
          if (const auto* hist = dynamic_cast<const Hist1D*>(obj); hist)
            for (const auto& value : hist->values())
              series[hist->name()].emplace_back(value);
        if (obj->isGraph1D())
          if (const auto* graph = dynamic_cast<const Graph1D*>(obj); graph)
            for (const auto& [xval, yval] : graph->points())
              series[graph->name()].emplace_back(yval);
      }
      ascii::Asciichart chart(series);
      CG_LOG << chart.type(type_).height(height_).Plot();
      return *this;
    }

  private:
    const ascii::Asciichart::Type type_;
    const int height_;
    const bool show_legend_;
  };
}  // namespace cepgen::utils
using cepgen::utils::asciichartDrawer;
REGISTER_DRAWER("asciichart", asciichartDrawer);
