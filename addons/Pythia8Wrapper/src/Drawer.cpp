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

#include <Pythia8/Basics.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Version.h"

namespace cepgen::pythia8 {
  class Drawer final : public utils::Drawer {
  public:
    explicit Drawer(const ParametersList& params) : utils::Drawer(params), hist_plot_(steer<bool>("histPlot")) {}

    static ParametersDescription description() {
      auto desc = cepgen::utils::Drawer::description();
      desc.setDescription("Pythia 8 plotter");
      desc.add("histPlot", true).setDescription("write Python code that can generate a PDF file with the spectra?");
      return desc;
    }

    const Drawer& draw(const utils::Graph1D&, const Mode&) const override {
      CG_WARNING("pythia8:Drawer:draw") << "1D graph plotter not (yet) implemented.";
      return *this;
    }
    const Drawer& draw(const utils::Graph2D&, const Mode&) const override {
      CG_WARNING("pythia8:Drawer:draw") << "2D graph plotter not (yet) implemented.";
      return *this;
    }
    const Drawer& draw(const utils::Hist1D& hist, const Mode& mode) const override {
      const auto out = convert(hist, mode);
      CG_LOG << out;
      if (hist_plot_) {
        Pythia8::HistPlot hp(hist.name());
        hp.plotFrame(
            "plot", out, hist.title(), hist.xAxis().label(), hist.yAxis().label(), "h", "void", mode & Mode::logy);
      }
      return *this;
    }
    const Drawer& draw(const utils::Hist2D&, const Mode&) const override {
      CG_WARNING("pythia8:Drawer:draw") << "Not yet implemented.";
      return *this;
    }

    const Drawer& draw(const utils::DrawableColl& objs,
                       const std::string& name = "",
                       const std::string& title = "",
                       const Mode& mode = Mode::none) const override {
      if (!hist_plot_) {
        CG_WARNING("pythia8:Drawer:draw") << "Not yet implemented.";
        return *this;
      }
      std::vector<Pythia8::Hist> histograms;
      const utils::Drawable* first_histogram = nullptr;
      for (const auto* obj : objs)
        if (obj->isHist1D()) {
          auto* hist = dynamic_cast<const utils::Hist1D*>(obj);
          if (!first_histogram)
            first_histogram = hist;
          histograms.emplace_back(convert(*hist, mode));
        } else
          CG_WARNING("pythia8:Drawer:draw") << "Multi-plotter only supports 1D histograms.";
      if (histograms.empty())
        return *this;
      if (!first_histogram)
        throw CG_ERROR("pythia8:Drawer:draw") << "First histogram was not found in list of drawable objects.";
      Pythia8::HistPlot hp(name);
      hp.frame("plot", title, first_histogram->xAxis().label(), first_histogram->yAxis().label());
      for (const auto& hist : histograms)
        hp.add(hist);
      hp.plot(mode & Mode::logy);
      return *this;
    }

  private:
    Pythia8::Hist convert(const utils::Hist1D& hist, const Mode& mode) const {
      Pythia8::Hist out(hist.title(), hist.nbins(), hist.range().min(), hist.range().max(), mode & Mode::logx);
      for (size_t ibin = 0; ibin < hist.nbins(); ++ibin)
        out.fill(hist.binRange(ibin).x(0.5), hist.value(ibin));
      return out;
    }
    const bool hist_plot_;
  };
}  // namespace cepgen::pythia8
using Pythia8Drawer = cepgen::pythia8::Drawer;
REGISTER_DRAWER("pythia8", Pythia8Drawer);
