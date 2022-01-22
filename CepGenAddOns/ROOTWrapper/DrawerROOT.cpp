#include <TGraph2DErrors.h>
#include <TH2D.h>
#include <TMultiGraph.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

namespace cepgen {
  namespace utils {
    class DrawerROOT : public Drawer {
    public:
      explicit DrawerROOT(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.add<std::string>("filename", "canvas").setDescription("default filename for the output");
        desc.add<std::string>("format", "pdf").setDescription("default extension for the output");
        desc.add<int>("palette", kLightTemperature).setDescription("ROOT colour palette to use");
        return desc;
      }

      const DrawerROOT& draw(const Graph1D&, const Mode&) const override;
      const DrawerROOT& draw(const Graph2D&, const Mode&) const override;
      const DrawerROOT& draw(const Hist1D&, const Mode&) const override;
      const DrawerROOT& draw(const Hist2D&, const Mode&) const override;

      const DrawerROOT& draw(const DrawableColl&,
                             const std::string& name = "",
                             const Mode& mode = Mode::none) const override;

    private:
      static void setMode(ROOTCanvas&, const Mode&);
      static TGraphErrors convert(const Graph1D&);
      static TGraph2DErrors convert(const Graph2D&);
      static TH1D convert(const Hist1D&);
      static TH2D convert(const Hist2D&);

      const std::string def_filename_;
      const std::string def_extension_;
    };

    DrawerROOT::DrawerROOT(const ParametersList& params)
        : Drawer(params), def_filename_(steer<std::string>("filename")), def_extension_(steer<std::string>("format")) {
      gStyle->SetPalette(steer<int>("palette"));
    }

    const DrawerROOT& DrawerROOT::draw(const Graph1D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), graph.title());
      setMode(canv, mode);
      gr.Draw("al");
      gr.GetHistogram()->SetTitle((";" + graph.xLabel() + ";" + graph.yLabel()).c_str());
      canv.Prettify(gr.GetHistogram());
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Graph2D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), graph.title());
      setMode(canv, mode);
      gr.Draw("surf3");
      gr.GetHistogram()->SetTitle((";" + graph.xLabel() + ";" + graph.yLabel() + ";" + graph.zLabel()).c_str());
      canv.Prettify(gr.GetHistogram());
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist1D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), hist.title());
      setMode(canv, mode);
      h.Draw();
      canv.Prettify(&h);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist2D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), hist.title());
      setMode(canv, mode);
      h.Draw("colz");
      canv.Prettify(&h);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const DrawableColl& objs, const std::string& name, const Mode& mode) const {
      TMultiGraph mg;
      THStack hs;
      ROOTCanvas canv(name.empty() ? def_filename_ : name, "");
      setMode(canv, mode);
      size_t i = 0;
      for (const auto* obj : objs) {
        if (obj->isHist1D()) {
          auto* hist = new TH1D(convert(*dynamic_cast<const Hist1D*>(obj)));
          hist->SetLineColor(ROOTCanvas::colours.at(i++));
          hs.Add(hist);
          canv.AddLegendEntry(hist, obj->title(), "l");
        } else if (obj->isGraph1D()) {
          auto* gr = new TGraphErrors(convert(*dynamic_cast<const Graph1D*>(obj)));
          gr->SetLineColor(ROOTCanvas::colours.at(i++));
          mg.Add(gr);
          canv.AddLegendEntry(gr, obj->title(), "l");
        } else {
          CG_WARNING("DrawerROOT:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
          continue;
        }
      }
      const bool has_hists = hs.GetHists() && !hs.GetHists()->IsEmpty();
      const bool has_graphs = mg.GetListOfGraphs() && !mg.GetListOfGraphs()->IsEmpty();
      if (has_hists)
        hs.Draw(mode & Mode::nostack ? "nostack" : "");
      if (has_graphs)
        mg.Draw((std::string("l") + (!has_hists ? "a" : "")).c_str());
      if (has_hists)
        canv.Prettify(hs.GetHistogram());
      else if (has_graphs)
        canv.Prettify(mg.GetHistogram());
      canv.Save(def_extension_);
      return *this;
    }

    void DrawerROOT::setMode(ROOTCanvas& canv, const Mode& mode) {
      if (mode & Mode::logx)
        canv.SetLogx();
      if (mode & Mode::logy)
        canv.SetLogy();
      if (mode & Mode::logz)
        canv.SetLogz();
    }

    TGraphErrors DrawerROOT::convert(const Graph1D& graph) {
      TGraphErrors gr;
      gr.SetTitle(graph.title().c_str());
      int i = 0;
      for (const auto& it : graph.points()) {
        gr.SetPoint(i, it.first.value, it.second.value);
        gr.SetPointError(i, 0. /* FIXME not yet supported */, it.second.value_unc);
        ++i;
      }
      return gr;
    }

    TGraph2DErrors DrawerROOT::convert(const Graph2D& graph) {
      TGraph2DErrors gr;
      gr.SetTitle(graph.title().c_str());
      int i = 0;
      for (const auto& it_x : graph.points()) {
        const auto& ax_x = it_x.first.value;
        for (const auto& it_y : it_x.second) {
          const auto& ax_y = it_y.first.value;
          gr.SetPoint(i, ax_x, ax_y, it_y.second.value);
          gr.SetPointError(i, 0., 0., it_y.second.value_unc);
          ++i;
        }
      }
      return gr;
    }

    TH1D DrawerROOT::convert(const Hist1D& hist) {
      const auto& rng = hist.range();
      TH1D h(hist.name().c_str(), hist.title().c_str(), hist.nbins(), rng.min(), rng.max());
      for (size_t i = 0; i < hist.nbins(); ++i) {
        h.SetBinContent(i + 1, hist.value(i));
        h.SetBinError(i + 1, hist.valueUnc(i));
      }
      h.GetXaxis()->SetTitle(hist.xLabel().c_str());
      h.GetYaxis()->SetTitle(hist.yLabel().c_str());
      return h;
    }

    TH2D DrawerROOT::convert(const Hist2D& hist) {
      const auto &rng_x = hist.rangeX(), &rng_y = hist.rangeY();
      TH2D h(hist.name().c_str(),
             hist.title().c_str(),
             hist.nbinsX(),
             rng_x.min(),
             rng_x.max(),
             hist.nbinsY(),
             rng_y.min(),
             rng_y.max());
      for (size_t ix = 0; ix < hist.nbinsX(); ++ix)
        for (size_t iy = 0; iy < hist.nbinsY(); ++iy) {
          h.SetBinContent(ix + 1, iy + 1, hist.value(ix, iy));
          h.SetBinError(ix + 1, iy + 1, hist.valueUnc(ix, iy));
        }
      h.GetXaxis()->SetTitle(hist.xLabel().c_str());
      h.GetYaxis()->SetTitle(hist.yLabel().c_str());
      h.GetZaxis()->SetTitle(hist.zLabel().c_str());
      return h;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("root", DrawerROOT)
