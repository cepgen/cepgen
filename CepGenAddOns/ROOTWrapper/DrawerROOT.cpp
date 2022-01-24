#include <TGraph2DErrors.h>
#include <TH2D.h>
#include <TMultiGraph.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"
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

    DrawerROOT::DrawerROOT(const ParametersList& params)
        : Drawer(params), def_filename_(steer<std::string>("filename")), def_extension_(steer<std::string>("format")) {
      gStyle->SetPalette(steer<int>("palette"));
    }

    const DrawerROOT& DrawerROOT::draw(const Graph1D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle());
      setMode(canv, mode);
      gr.Draw("al");
      gr.GetHistogram()->SetTitle(delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label()));
      canv.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Graph2D& graph, const Mode& mode) const {
      auto gr = convert(graph);
      ROOTCanvas canv(graph.name().empty() ? def_filename_ : graph.name(), gr.GetTitle());
      setMode(canv, mode);
      gr.Draw("surf3");
      gr.GetHistogram()->SetTitle(
          delatexify(";" + graph.xAxis().label() + ";" + graph.yAxis().label() + ";" + graph.zAxis().label()));
      canv.Prettify(gr.GetHistogram());
      postDraw(gr.GetHistogram(), graph);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist1D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), h.GetTitle());
      setMode(canv, mode);
      h.Draw();
      canv.Prettify(&h);
      postDraw(&h, hist);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist2D& hist, const Mode& mode) const {
      auto h = convert(hist);
      ROOTCanvas canv(hist.name().empty() ? def_filename_ : hist.name(), h.GetTitle());
      setMode(canv, mode);
      h.Draw("colz");
      canv.Prettify(&h);
      postDraw(&h, hist);
      canv.Save(def_extension_);
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const DrawableColl& objs,
                                       const std::string& name,
                                       const std::string& title,
                                       const Mode& mode) const {
      TMultiGraph mg;
      THStack hs;
      ROOTCanvas canv(name.empty() ? def_filename_ : name, delatexify(title).Data());
      setMode(canv, mode);
      size_t i = 0;
      Drawable* first = nullptr;
      for (const auto* obj : objs) {
        if (obj->isHist1D()) {
          auto* hist = new TH1D(convert(*dynamic_cast<const Hist1D*>(obj)));
          hist->SetLineColor(ROOTCanvas::colours.at(i++));
          hs.Add(hist);
          canv.AddLegendEntry(hist, hist->GetTitle(), "l");
        } else if (obj->isGraph1D()) {
          auto* gr = new TGraphErrors(convert(*dynamic_cast<const Graph1D*>(obj)));
          gr->SetLineColor(ROOTCanvas::colours.at(i++));
          mg.Add(gr);
          canv.AddLegendEntry(gr, gr->GetTitle(), "l");
        } else {
          CG_WARNING("DrawerROOT:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
          continue;
        }
        if (!first)
          first = const_cast<Drawable*>(obj);
      }
      const bool has_hists = hs.GetHists() && !hs.GetHists()->IsEmpty();
      const bool has_graphs = mg.GetListOfGraphs() && !mg.GetListOfGraphs()->IsEmpty();
      if (has_hists)
        hs.Draw(mode & Mode::nostack ? "nostack" : "");
      if (has_graphs)
        mg.Draw((std::string("l") + (!has_hists ? "a" : "")).c_str());
      if (has_hists) {
        canv.Prettify(hs.GetHistogram());
        postDraw(hs.GetHistogram(), *first);
      } else if (has_graphs) {
        canv.Prettify(mg.GetHistogram());
        postDraw(mg.GetHistogram(), *first);
      }
      canv.Save(def_extension_);
      return *this;
    }

    void DrawerROOT::setMode(ROOTCanvas& canv, const Mode& mode) {
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

    void DrawerROOT::postDraw(TH1* obj, const Drawable& dr) {
      const auto &xrng = dr.xAxis().range(), &yrng = dr.yAxis().range();
      obj->GetXaxis()->SetTitle(delatexify(dr.xAxis().label()));
      obj->GetYaxis()->SetTitle(delatexify(dr.yAxis().label()));
      if (xrng.valid())
        obj->GetXaxis()->SetRangeUser(xrng.min(), xrng.max());
      if (yrng.valid()) {
        if (yrng.hasMin())
          obj->SetMinimum(yrng.min());
        if (yrng.hasMax())
          obj->SetMaximum(yrng.max());
      }
    }

    TString DrawerROOT::delatexify(const std::string& tok) {
      auto out = utils::replace_all(tok, {{"$", ""}});
      return TString(out);
    }

    TGraphErrors DrawerROOT::convert(const Graph1D& graph) {
      TGraphErrors gr;
      gr.SetTitle(delatexify(graph.title()));
      int i = 0;
      for (const auto& it : graph.points()) {
        gr.SetPoint(i, it.first.value, it.second.value);
        gr.SetPointError(i, it.first.value_unc, it.second.value_unc);
        ++i;
      }
      gr.SetLineWidth(2);
      return gr;
    }

    TGraph2DErrors DrawerROOT::convert(const Graph2D& graph) {
      TGraph2DErrors gr;
      gr.SetTitle(delatexify(graph.title()));
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
      TH1D h(hist.name().c_str(), delatexify(hist.title()), hist.nbins(), rng.min(), rng.max());
      for (size_t i = 0; i < hist.nbins(); ++i) {
        h.SetBinContent(i + 1, hist.value(i));
        h.SetBinError(i + 1, hist.valueUnc(i));
      }
      h.GetXaxis()->SetTitle(delatexify(hist.xAxis().label()));
      h.GetYaxis()->SetTitle(delatexify(hist.yAxis().label()));
      return h;
    }

    TH2D DrawerROOT::convert(const Hist2D& hist) {
      const auto &rng_x = hist.rangeX(), &rng_y = hist.rangeY();
      TH2D h(hist.name().c_str(),
             delatexify(hist.title()),
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
      h.GetXaxis()->SetTitle(delatexify(hist.xAxis().label()));
      h.GetYaxis()->SetTitle(delatexify(hist.yAxis().label()));
      h.GetZaxis()->SetTitle(delatexify(hist.zAxis().label()));
      return h;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("root", DrawerROOT)
