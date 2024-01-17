/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include <TGraph2DErrors.h>
#include <TH2D.h>
#include <TMultiGraph.h>

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

namespace cepgen {
  namespace utils {
    class ROOTDrawer : public Drawer {
    public:
      explicit ROOTDrawer(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.add<std::string>("filename", "canvas").setDescription("default filename for the output");
        desc.add<std::string>("format", "pdf").setDescription("default extension for the output");
        desc.add<int>("palette", kLightTemperature).setDescription("ROOT colour palette to use");
        return desc;
      }

      const ROOTDrawer& draw(const Graph1D&, const Mode&) const override;
      const ROOTDrawer& draw(const Graph2D&, const Mode&) const override;
      const ROOTDrawer& draw(const Hist1D&, const Mode&) const override;
      const ROOTDrawer& draw(const Hist2D&, const Mode&) const override;

      const ROOTDrawer& draw(const DrawableColl&,
                             const std::string& name = "",
                             const std::string& title = "",
                             const Mode& mode = Mode::none) const override;

    private:
      static void setMode(ROOTCanvas&, const Mode&);
      static void postDraw(TH1*, const Drawable&);
      static TString delatexify(const std::string&);
      static TGraphErrors convert(const Graph1D&);
      static TGraph2DErrors convert(const Graph2D&);
      static TH1D convert(const Hist1D&);
      static TH2D convert(const Hist2D&);

      const std::string def_filename_;
      const std::string def_extension_;
    };

    ROOTDrawer::ROOTDrawer(const ParametersList& params)
        : Drawer(params), def_filename_(steer<std::string>("filename")), def_extension_(steer<std::string>("format")) {
      gStyle->SetPalette(steer<int>("palette"));
    }

    const ROOTDrawer& ROOTDrawer::draw(const Graph1D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle(), mode & Mode::ratio);
      setMode(canv, mode);
      gr.Draw("al");
      gr.GetHistogram()->SetTitle(delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label()));
      canv.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canv.Save(def_extension_);
      return *this;
    }

    const ROOTDrawer& ROOTDrawer::draw(const Graph2D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle(), mode & Mode::ratio);
      setMode(canv, mode);
      if (mode & Mode::col)
        gr.Draw("colz");
      else if (mode & Mode::cont)
        gr.Draw("cont");
      else
        gr.Draw("surf3");
      gr.GetHistogram()->SetTitle(
          delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label() + ";" + graph.zAxis().label()));
      canv.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canv.Save(def_extension_);
      return *this;
    }

    const ROOTDrawer& ROOTDrawer::draw(const Hist1D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), h.GetTitle(), mode & Mode::ratio);
      setMode(canv, mode);
      h.Draw();
      canv.Prettify(&h);
      postDraw(&h, hist);
      canv.Save(def_extension_);
      return *this;
    }

    const ROOTDrawer& ROOTDrawer::draw(const Hist2D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), h.GetTitle(), mode & Mode::ratio);
      setMode(canv, mode);
      h.Draw("colz");
      canv.Prettify(&h);
      postDraw(&h, hist);
      canv.Save(def_extension_);
      return *this;
    }

    const ROOTDrawer& ROOTDrawer::draw(const DrawableColl& objs,
                                       const std::string& name,
                                       const std::string& title,
                                       const Mode& mode) const {
      ROOTCanvas canv(name.empty() ? def_filename_ : name, delatexify(title).Data(), mode & Mode::ratio);
      auto* mg = canv.Make<TMultiGraph>();
      auto* hs = canv.Make<THStack>();
      setMode(canv, mode);
      Drawable* first = nullptr;
      DrawableColl plots_2d;
      for (size_t i = 0; i < objs.size(); ++i) {
        const auto* obj = objs.at(i);
        auto colour = ROOTCanvas::colours.at(i % ROOTCanvas::colours.size());
        auto style = i + 1;
        if (obj->isHist1D()) {
          auto* hist = new TH1D(convert(*dynamic_cast<const Hist1D*>(obj)));
          hist->SetLineColor(colour);
          hist->SetLineStyle(style);
          hs->Add(hist);
          canv.AddLegendEntry(hist, hist->GetTitle(), "l");
        } else if (obj->isGraph1D()) {
          auto* gr = new TGraphErrors(convert(*dynamic_cast<const Graph1D*>(obj)));
          gr->SetLineColor(colour);
          gr->SetLineStyle(style);
          mg->Add(gr);
          canv.AddLegendEntry(gr, gr->GetTitle(), "l");
        } else {
          plots_2d.emplace_back(obj);
          CG_DEBUG("ROOTDrawer:draw") << "Adding a 2-dimensional drawable '" << obj->name() << "' to the stack.";
          continue;
        }
        if (!first)
          first = const_cast<Drawable*>(obj);
      }
      const bool has_hists = hs->GetHists() && !hs->GetHists()->IsEmpty();
      const bool has_graphs = mg->GetListOfGraphs() && !mg->GetListOfGraphs()->IsEmpty();
      if (has_hists || has_graphs) {
        if (has_hists)
          hs->Draw(mode & Mode::nostack ? "nostack" : "");
        if (has_graphs)
          mg->Draw((std::string("l") + (!has_hists ? "a" : "")).c_str());
        if (has_hists) {
          postDraw(hs->GetHistogram(), *first);
          canv.Prettify(hs);
        } else if (has_graphs) {
          postDraw(mg->GetHistogram(), *first);
          canv.Prettify(mg);
        }
        canv.Save(def_extension_);
      }
      for (size_t i = 0; i < plots_2d.size(); ++i) {
        const auto* obj = plots_2d.at(i);
        const std::string postfix = i == 0 ? "(" : i == plots_2d.size() - 1 ? ")" : "";
        if (obj->isHist2D()) {
          const auto* hist = dynamic_cast<const Hist2D*>(obj);
          auto* h = new TH2D(convert(*hist));
          setMode(canv, mode);
          h->Draw("colz");
          canv.Prettify(h);
          postDraw(h, *hist);
        } else if (obj->isGraph2D()) {
          const auto* graph = dynamic_cast<const Graph2D*>(obj);
          auto* gr = new TGraph2D(convert(*graph));
          setMode(canv, mode);
          if (mode & Mode::col)
            gr->Draw("colz");
          else if (mode & Mode::cont)
            gr->Draw("cont");
          else
            gr->Draw("surf3");
          gr->GetHistogram()->SetTitle(
              delatexify(";" + graph->xAxis().label() + ";" + graph->yAxis().label() + ";" + graph->zAxis().label()));
          canv.Prettify(gr->GetHistogram());
          postDraw(gr->GetHistogram(), *graph);
        }
        canv.Print(utils::format("%s_multi.%s%s", canv.GetName(), def_extension_.data(), postfix.data()).data());
      }
      return *this;
    }

    void ROOTDrawer::setMode(ROOTCanvas& canv, const Mode& mode) {
      canv.SetLegendX1(0.175);
      if (mode & Mode::logx)
        canv.SetLogx();
      if (mode & Mode::logy)
        canv.SetLogy();
      if (mode & Mode::logz)
        canv.SetLogz();
      if (mode & Mode::grid)
        canv.SetGrid();
    }

    void ROOTDrawer::postDraw(TH1* obj, const Drawable& dr) {
      const auto &xrng = dr.xAxis().range(), &yrng = dr.yAxis().range();
      obj->GetXaxis()->SetTitle(delatexify(dr.xAxis().label()));
      obj->GetYaxis()->SetTitle(delatexify(dr.yAxis().label()));
      if (xrng.valid())
        obj->GetXaxis()->SetLimits(xrng.min(), xrng.max());
      if (yrng.valid()) {
        if (yrng.hasMin())
          obj->SetMinimum(yrng.min());
        if (yrng.hasMax())
          obj->SetMaximum(yrng.max());
      }
    }

    TString ROOTDrawer::delatexify(const std::string& tok) {
      auto out = utils::replace_all(tok, {{"$", ""}});
      return TString(out);
    }

    TGraphErrors ROOTDrawer::convert(const Graph1D& graph) {
      TGraphErrors gr;
      gr.SetTitle(delatexify(graph.title()));
      int i = 0;
      for (const auto& it : graph.points()) {
        const auto& val = it.second;
        gr.SetPoint(i, it.first.value, val);
        gr.SetPointError(i, it.first.value_unc, val.uncertainty());
        ++i;
      }
      gr.SetLineWidth(3);
      return gr;
    }

    TGraph2DErrors ROOTDrawer::convert(const Graph2D& graph) {
      TGraph2DErrors gr;
      gr.SetTitle(delatexify(graph.title()));
      int i = 0;
      for (const auto& it_x : graph.points()) {
        const auto& ax_x = it_x.first.value;
        for (const auto& it_y : it_x.second) {
          const auto& ax_y = it_y.first.value;
          const auto& val = it_y.second;
          gr.SetPoint(i, ax_x, ax_y, val);
          gr.SetPointError(i, 0., 0., val.uncertainty());
          ++i;
        }
      }
      return gr;
    }

    TH1D ROOTDrawer::convert(const Hist1D& hist) {
      const auto bins = hist.bins(Histogram::BinMode::both);
      TH1D h(hist.name().c_str(), delatexify(hist.title()), bins.size() - 1, bins.data());
      h.SetBinContent(0, hist.underflow());
      for (size_t i = 0; i < hist.nbins(); ++i) {
        const auto val = hist.value(i);
        h.SetBinContent(i + 1, val);
        h.SetBinError(i + 1, val.uncertainty());
      }
      h.SetBinContent(hist.nbins() + 1, hist.overflow());
      h.GetXaxis()->SetTitle(delatexify(hist.xAxis().label()));
      h.GetYaxis()->SetTitle(delatexify(hist.yAxis().label()));
      return h;
    }

    TH2D ROOTDrawer::convert(const Hist2D& hist) {
      const auto bins_x = hist.binsX(Histogram::BinMode::both), bins_y = hist.binsY(Histogram::BinMode::both);
      TH2D h(hist.name().c_str(),
             delatexify(hist.title()),
             bins_x.size() - 1,
             bins_x.data(),
             bins_y.size() - 1,
             bins_y.data());
      for (size_t ix = 0; ix < hist.nbinsX(); ++ix)
        for (size_t iy = 0; iy < hist.nbinsY(); ++iy) {
          const auto val = hist.value(ix, iy);
          h.SetBinContent(ix + 1, iy + 1, val);
          h.SetBinError(ix + 1, iy + 1, val.uncertainty());
        }
      h.SetBinContent(0, 0, hist.outOfRange().at(Hist2D::contents_t::LT_LT));
      h.SetBinContent(0, 1, hist.outOfRange().at(Hist2D::contents_t::LT_IN));
      h.SetBinContent(0, hist.nbinsY() + 1, hist.outOfRange().at(Hist2D::contents_t::LT_GT));
      h.SetBinContent(1, 0, hist.outOfRange().at(Hist2D::contents_t::IN_LT));
      h.SetBinContent(1, hist.nbinsY() + 1, hist.outOfRange().at(Hist2D::contents_t::IN_GT));
      h.SetBinContent(hist.nbinsX() + 1, 0, hist.outOfRange().at(Hist2D::contents_t::GT_LT));
      h.SetBinContent(hist.nbinsX() + 1, 1, hist.outOfRange().at(Hist2D::contents_t::GT_IN));
      h.SetBinContent(hist.nbinsX() + 1, hist.nbinsY() + 1, hist.outOfRange().at(Hist2D::contents_t::GT_GT));
      h.GetXaxis()->SetTitle(delatexify(hist.xAxis().label()));
      h.GetYaxis()->SetTitle(delatexify(hist.yAxis().label()));
      h.GetZaxis()->SetTitle(delatexify(hist.zAxis().label()));
      return h;
    }
  }  // namespace utils
}  // namespace cepgen
typedef cepgen::utils::ROOTDrawer RootDrawer;
REGISTER_DRAWER("root", RootDrawer);
