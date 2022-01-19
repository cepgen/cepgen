#include <YODA/Histo1D.h>
#include <YODA/Histo2D.h>
#include <YODA/Scatter2D.h>
#include <YODA/Scatter3D.h>
#include <YODA/WriterAIDA.h>
#include <YODA/WriterFLAT.h>
#include <YODA/WriterYODA.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    template <typename T>
    class DrawerYODA : public Drawer {
    public:
      explicit DrawerYODA(const ParametersList& params)
          : Drawer(params), file_(steer<std::string>("filename")), writer_(&T::create()) {
        if (steer<bool>("compress"))
          writer_->useCompression(true);
        writer_->setPrecision(steer<int>("precision"));
      }

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("YODA/AIDA plotting utility");
        desc.add<std::string>("filename", "plots.yoda");
        desc.add<bool>("compress", false).setDescription("use libz compression?");
        desc.add<int>("precision", 6).setDescription("precision of numerical quantities in output");
        return desc;
      }

      const DrawerYODA& draw(const Graph1D& graph, const Mode&) const override {
        writer_->write(file_, convert(graph));
        return *this;
      }
      const DrawerYODA& draw(const Graph2D& graph, const Mode&) const override {
        writer_->write(file_, convert(graph));
        return *this;
      }
      const DrawerYODA& draw(const Hist1D& hist, const Mode&) const override {
        writer_->write(file_, convert(hist));
        return *this;
      }
      const DrawerYODA& draw(const Hist2D& hist, const Mode&) const override {
        writer_->write(file_, convert(hist));
        return *this;
      }

      const DrawerYODA& draw(const DrawableColl&,
                             const std::string& name = "",
                             const Mode& mode = Mode::none) const override;

    private:
      static YODA::Scatter2D convert(const Graph1D&);
      static YODA::Scatter3D convert(const Graph2D&);
      static YODA::Histo1D convert(const Hist1D&);
      static YODA::Histo2D convert(const Hist2D&);
      static std::string path(const std::string& name) {
        return "/" + utils::replace_all(utils::replace_all(name, ")", "_"), "(", "_");
      }
      mutable std::ofstream file_;
      mutable YODA::Writer* writer_;
    };

    template <typename T>
    const DrawerYODA<T>& DrawerYODA<T>::draw(const DrawableColl& objs, const std::string&, const Mode&) const {
      std::vector<const YODA::AnalysisObject*> objs_coll;
      for (const auto* obj : objs) {
        if (obj->isHist1D()) {
          objs_coll.emplace_back(convert(*dynamic_cast<const Hist1D*>(obj)).newclone());
        } else if (obj->isGraph1D()) {
          objs_coll.emplace_back(convert(*dynamic_cast<const Graph1D*>(obj)).newclone());
        } else {
          CG_WARNING("DrawerYODA:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
          continue;
        }
      }
      writer_->write(file_, objs_coll);
      return *this;
    }

    template <typename T>
    YODA::Scatter2D DrawerYODA<T>::convert(const Graph1D& graph) {
      YODA::Scatter2D gr(path(graph.name()), graph.title());
      for (const auto& it : graph.points())
        gr.addPoint(it.first.value, it.second.value, 0. /* FIXME not yet supported */, it.second.value_unc);
      //gr.setAnnotation("xlabel", graph.xAxis().label());
      //gr.setAnnotation("ylabel", graph.yAxis().label());
      return gr;
    }

    template <typename T>
    YODA::Scatter3D DrawerYODA<T>::convert(const Graph2D& graph) {
      YODA::Scatter3D gr(path(graph.name()), graph.title());
      for (const auto& it_x : graph.points()) {
        const auto& ax_x = it_x.first.value;
        for (const auto& it_y : it_x.second) {
          const auto& ax_y = it_y.first.value;
          gr.addPoint(ax_x, ax_y, it_y.second.value, 0., 0., it_y.second.value_unc);
        }
      }
      //gr.setAnnotation("xlabel", graph.xAxis().label());
      //gr.setAnnotation("ylabel", graph.yAxis().label());
      return gr;
    }

    template <typename T>
    YODA::Histo1D DrawerYODA<T>::convert(const Hist1D& hist) {
      const auto& rng = hist.range();
      YODA::Histo1D h(hist.nbins(), rng.min(), rng.max(), path(hist.name()), hist.title());
      for (size_t i = 0; i < hist.nbins(); ++i)
        h.fillBin(i, hist.value(i), std::pow(hist.valueUnc(i), 2));
      //h.setAnnotation("xlabel", hist.xAxis().label());
      //h.setAnnotation("ylabel", hist.yAxis().label());
      return h;
    }

    template <typename T>
    YODA::Histo2D DrawerYODA<T>::convert(const Hist2D& hist) {
      const auto &rng_x = hist.rangeX(), &rng_y = hist.rangeY();
      YODA::Histo2D h(hist.nbinsX(),
                      rng_x.min(),
                      rng_x.max(),
                      hist.nbinsY(),
                      rng_y.min(),
                      rng_y.max(),
                      path(hist.name()),
                      hist.title());
      for (size_t ix = 0; ix < hist.nbinsX(); ++ix)
        for (size_t iy = 0; iy < hist.nbinsY(); ++iy)
          h.fillBin((ix + 1) * (iy + 1), hist.value(ix, iy), std::pow(hist.valueUnc(ix, iy), 2));
      //h.setAnnotation("xlabel", hist.xAxis().label());
      //h.setAnnotation("ylabel", hist.yAxis().label());
      return h;
    }
  }  // namespace utils
}  // namespace cepgen

typedef cepgen::utils::DrawerYODA<YODA::WriterYODA> DrawerYoda;
typedef cepgen::utils::DrawerYODA<YODA::WriterAIDA> DrawerYodaAida;
typedef cepgen::utils::DrawerYODA<YODA::WriterFLAT> DrawerYodaFlat;
REGISTER_DRAWER("yoda", DrawerYoda)
REGISTER_DRAWER("yoda_aida", DrawerYodaAida)
REGISTER_DRAWER("yoda_flat", DrawerYodaFlat)
