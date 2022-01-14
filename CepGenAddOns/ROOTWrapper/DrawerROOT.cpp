#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawable.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "TGraph2DErrors.h"
#include "TH2D.h"

namespace cepgen {
  namespace utils {
    class DrawerROOT : public Drawer {
    public:
      explicit DrawerROOT(const ParametersList&);

      const DrawerROOT& draw(const Graph1D&, const Mode&) const override;
      const DrawerROOT& draw(const Graph2D&, const Mode&) const override;
      const DrawerROOT& draw(const Hist1D&, const Mode&) const override;
      const DrawerROOT& draw(const Hist2D&, const Mode&) const override;

    private:
      static TGraphErrors convert(const Graph1D&);
      static TGraph2DErrors convert(const Graph2D&);
      static TH1D convert(const Hist1D&);
      static TH2D convert(const Hist2D&);
    };

    DrawerROOT::DrawerROOT(const ParametersList& params) : Drawer(params) {}

    const DrawerROOT& DrawerROOT::draw(const Graph1D& graph, const Mode&) const {
      auto gr = convert(graph);
      cepgen::ROOTCanvas canv(graph.name(), graph.title());
      gr.Draw("al");
      canv.Prettify(gr.GetHistogram());
      canv.Save("pdf");
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Graph2D& graph, const Mode&) const {
      auto gr = convert(graph);
      cepgen::ROOTCanvas canv(graph.name());
      gr.Draw("colz");
      canv.Prettify(gr.GetHistogram());
      canv.Save("pdf");
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist1D& hist, const Mode&) const {
      auto h = convert(hist);
      cepgen::ROOTCanvas canv(hist.name(), hist.title());
      h.Draw();
      canv.Prettify(&h);
      canv.Save("pdf");
      return *this;
    }

    const DrawerROOT& DrawerROOT::draw(const Hist2D& hist, const Mode&) const {
      auto h = convert(hist);
      cepgen::ROOTCanvas canv(hist.name(), hist.title());
      h.Draw("colz");
      canv.Prettify(&h);
      canv.Save("pdf");
      return *this;
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
      return h;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("root", DrawerROOT)
