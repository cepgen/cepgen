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

#include <TGraph2DErrors.h>
#include <TH2D.h>
#include <TMultiGraph.h>

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGenRoot/ROOTCanvas.h"

using namespace std::string_literals;

namespace cepgen::root {
  class Drawer final : public utils::Drawer {
  public:
    explicit Drawer(const ParametersList& params)
        : utils::Drawer(params),
          def_filename_(steer<std::string>("filename")),
          def_extension_(steer<std::string>("format")) {
      gStyle->SetPalette(steer<int>("palette"));
    }

    static ParametersDescription description() {
      auto desc = utils::Drawer::description();
      desc.add("filename", "canvas"s).setDescription("default filename for the output");
      desc.add("format", "pdf"s).setDescription("default extension for the output");
      desc.add<int>("palette", kLightTemperature).setDescription("ROOT colour palette to use");
      return desc;
    }

    const Drawer& draw(const utils::Graph1D& graph, const Mode& mode) const override {
      auto gr = convert(graph);
      ROOTCanvas canvas(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle(), mode & Mode::ratio);
      setMode(canvas, mode);
      gr.Draw("al");
      gr.GetHistogram()->SetTitle(delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label()));
      canvas.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canvas.Save(def_extension_);
      return *this;
    }
    const Drawer& draw(const utils::Graph2D& graph, const Mode& mode) const override {
      auto gr = convert(graph);
      ROOTCanvas canvas(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle(), mode & Mode::ratio);
      setMode(canvas, mode);
      if (mode & Mode::col)
        gr.Draw("colz");
      else if (mode & Mode::cont)
        gr.Draw("cont");
      else
        gr.Draw("surf3");
      gr.GetHistogram()->SetTitle(
          delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label() + ";" + graph.zAxis().label()));
      canvas.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canvas.Save(def_extension_);
      return *this;
    }
    const Drawer& draw(const utils::Hist1D& histogram, const Mode& mode) const override {
      auto h = convert(histogram);
      ROOTCanvas canvas(histogram.name().empty() ? def_filename_ : histogram.name(), h.GetTitle(), mode & Mode::ratio);
      setMode(canvas, mode);
      h.Draw();
      canvas.Prettify(&h);
      postDraw(&h, histogram);
      canvas.Save(def_extension_);
      return *this;
    }
    const Drawer& draw(const utils::Hist2D& histogram, const Mode& mode) const override {
      auto h = convert(histogram);
      ROOTCanvas canvas(histogram.name().empty() ? def_filename_ : histogram.name(), h.GetTitle(), mode & Mode::ratio);
      setMode(canvas, mode);
      h.Draw("colz");
      canvas.Prettify(&h);
      postDraw(&h, histogram);
      canvas.Save(def_extension_);
      return *this;
    }

    const Drawer& draw(const utils::DrawableColl&,
                       const std::string& name = "",
                       const std::string& title = "",
                       const Mode& mode = Mode::none) const override;

  private:
    static void setMode(ROOTCanvas& canvas, const Mode& mode) {
      canvas.SetLegendX1(0.175);
      if (mode & Mode::logx)
        canvas.SetLogx();
      if (mode & Mode::logy)
        canvas.SetLogy();
      if (mode & Mode::logz)
        canvas.SetLogz();
      if (mode & Mode::grid)
        canvas.SetGrid();
    }

    static void postDraw(TH1* histogram, const utils::Drawable& drawable) {
      const auto &x_range = drawable.xAxis().range(), &y_range = drawable.yAxis().range();
      histogram->GetXaxis()->SetTitle(delatexify(drawable.xAxis().label()));
      histogram->GetYaxis()->SetTitle(delatexify(drawable.yAxis().label()));
      histogram->SetLineWidth(std::max((short)3, histogram->GetLineWidth()));
      if (x_range.valid())
        histogram->GetXaxis()->SetLimits(x_range.min(), x_range.max());
      if (y_range.valid()) {
        if (y_range.hasMin())
          histogram->SetMinimum(y_range.min());
        if (y_range.hasMax())
          histogram->SetMaximum(y_range.max());
      }
    }
    static TString delatexify(const std::string& token) { return {utils::replaceAll(token, {{"$", ""}, {"\\", "#"}})}; }
    static TGraphErrors convert(const utils::Graph1D&);
    static TGraph2DErrors convert(const utils::Graph2D&);
    static TH1D convert(const utils::Hist1D&);
    static TH2D convert(const utils::Hist2D&);

    const std::string def_filename_;
    const std::string def_extension_;
  };

  const Drawer& Drawer::draw(const utils::DrawableColl& objects,
                             const std::string& name,
                             const std::string& title,
                             const Mode& mode) const {
    ROOTCanvas canvas(name.empty() ? def_filename_ : name, delatexify(title).Data(), mode & Mode::ratio);
    auto* mg = canvas.Make<TMultiGraph>();
    auto* hs = canvas.Make<THStack>();
    setMode(canvas, mode);
    utils::Drawable* first = nullptr;
    utils::DrawableColl plots_2d;
    for (size_t i = 0; i < objects.size(); ++i) {
      const auto* obj = objects.at(i);
      auto colour = ROOTCanvas::colours.at(i % ROOTCanvas::colours.size());
      auto style = i + 1;
      if (obj->isHist1D()) {
        if (auto* hist = new TH1D(convert(*dynamic_cast<const utils::Hist1D*>(obj))); hist) {
          hist->SetLineColor(colour);
          hist->SetLineStyle(style);
          hs->Add(hist);
          canvas.AddLegendEntry(hist, hist->GetTitle(), "l");
        }
      } else if (obj->isGraph1D()) {
        if (auto* gr = new TGraphErrors(convert(*dynamic_cast<const utils::Graph1D*>(obj))); gr) {
          gr->SetLineColor(colour);
          gr->SetLineStyle(style);
          mg->Add(gr);
          canvas.AddLegendEntry(gr, gr->GetTitle(), "l");
        }
      } else {
        plots_2d.emplace_back(obj);
        CG_DEBUG("root:Drawer:draw") << "Adding a 2-dimensional drawable '" << obj->name() << "' to the stack.";
        continue;
      }
      if (!first)
        first = const_cast<utils::Drawable*>(obj);
    }
    if (const bool has_hists = hs->GetHists() && !hs->GetHists()->IsEmpty(),
        has_graphs = mg->GetListOfGraphs() && !mg->GetListOfGraphs()->IsEmpty();
        has_hists || has_graphs) {
      if (has_hists)
        hs->Draw(((mode & Mode::bar ? "hist"s : ""s) + (mode & Mode::nostack ? "nostack"s : ""s)).data());
      if (has_graphs)
        mg->Draw(("l"s + (!has_hists ? "a" : "")).c_str());
      if (has_hists) {
        postDraw(hs->GetHistogram(), *first);
        canvas.Prettify(hs);
      } else if (has_graphs) {
        postDraw(mg->GetHistogram(), *first);
        canvas.Prettify(mg);
      }
      canvas.Save(def_extension_);
    }
    for (size_t i = 0; i < plots_2d.size(); ++i) {
      const auto* obj = plots_2d.at(i);
      const std::string postfix = i == 0 ? "(" : i == plots_2d.size() - 1 ? ")" : "";
      if (obj->isHist2D()) {
        if (const auto* hist = dynamic_cast<const utils::Hist2D*>(obj); hist) {
          auto* h = new TH2D(convert(*hist));
          setMode(canvas, mode);
          h->Draw("colz");
          canvas.Prettify(h);
          postDraw(h, *hist);
        }
      } else if (obj->isGraph2D()) {
        if (const auto* graph = dynamic_cast<const utils::Graph2D*>(obj); graph) {
          auto* gr = new TGraph2D(convert(*graph));
          setMode(canvas, mode);
          if (mode & Mode::col)
            gr->Draw("colz");
          else if (mode & Mode::cont)
            gr->Draw("cont");
          else
            gr->Draw("surf3");
          gr->GetHistogram()->SetTitle(
              delatexify(";" + graph->xAxis().label() + ";" + graph->yAxis().label() + ";" + graph->zAxis().label()));
          canvas.Prettify(gr->GetHistogram());
          postDraw(gr->GetHistogram(), *graph);
        }
      }
      canvas.Print(utils::format("%s_multi.%s%s", canvas.GetName(), def_extension_.data(), postfix.data()).data());
    }
    return *this;
  }

  TGraphErrors Drawer::convert(const utils::Graph1D& graph) {
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

  TGraph2DErrors Drawer::convert(const utils::Graph2D& graph) {
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

  TH1D Drawer::convert(const utils::Hist1D& histogram) {
    const auto bins = histogram.bins(utils::Histogram::BinMode::both);
    TH1D h(histogram.name().c_str(), delatexify(histogram.title()), bins.size() - 1, bins.data());
    h.SetBinContent(0, histogram.underflow());
    for (size_t i = 0; i < histogram.nbins(); ++i) {
      const auto val = histogram.value(i);
      h.SetBinContent(i + 1, val);
      h.SetBinError(i + 1, val.uncertainty());
    }
    h.SetBinContent(histogram.nbins() + 1, histogram.overflow());
    h.GetXaxis()->SetTitle(delatexify(histogram.xAxis().label()));
    h.GetYaxis()->SetTitle(delatexify(histogram.yAxis().label()));
    h.SetLineWidth(3);
    return h;
  }

  TH2D Drawer::convert(const utils::Hist2D& histogram) {
    const auto bins_x = histogram.binsX(utils::Histogram::BinMode::both),
               bins_y = histogram.binsY(utils::Histogram::BinMode::both);
    TH2D h(histogram.name().c_str(),
           delatexify(histogram.title()),
           bins_x.size() - 1,
           bins_x.data(),
           bins_y.size() - 1,
           bins_y.data());
    for (size_t ix = 0; ix < histogram.nbinsX(); ++ix)
      for (size_t iy = 0; iy < histogram.nbinsY(); ++iy) {
        const auto val = histogram.value(ix, iy);
        h.SetBinContent(ix + 1, iy + 1, val);
        h.SetBinError(ix + 1, iy + 1, val.uncertainty());
      }
    h.SetBinContent(0, 0, histogram.outOfRange().at(utils::Hist2D::contents_t::LT_LT));
    h.SetBinContent(0, 1, histogram.outOfRange().at(utils::Hist2D::contents_t::LT_IN));
    h.SetBinContent(0, histogram.nbinsY() + 1, histogram.outOfRange().at(utils::Hist2D::contents_t::LT_GT));
    h.SetBinContent(1, 0, histogram.outOfRange().at(utils::Hist2D::contents_t::IN_LT));
    h.SetBinContent(1, histogram.nbinsY() + 1, histogram.outOfRange().at(utils::Hist2D::contents_t::IN_GT));
    h.SetBinContent(histogram.nbinsX() + 1, 0, histogram.outOfRange().at(utils::Hist2D::contents_t::GT_LT));
    h.SetBinContent(histogram.nbinsX() + 1, 1, histogram.outOfRange().at(utils::Hist2D::contents_t::GT_IN));
    h.SetBinContent(
        histogram.nbinsX() + 1, histogram.nbinsY() + 1, histogram.outOfRange().at(utils::Hist2D::contents_t::GT_GT));
    h.GetXaxis()->SetTitle(delatexify(histogram.xAxis().label()));
    h.GetYaxis()->SetTitle(delatexify(histogram.yAxis().label()));
    h.GetZaxis()->SetTitle(delatexify(histogram.zAxis().label()));
    return h;
  }
}  // namespace cepgen::root
using cepgen::root::Drawer;
REGISTER_DRAWER("root", Drawer);
