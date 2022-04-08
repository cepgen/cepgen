#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Piper.h"
#include "CepGen/Utils/String.h"

#ifndef GNUPLOT_BIN
#error "Gnuplot executable must be specified using GNUPLOT_BIN!"
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define GNUPLOT TOSTRING(GNUPLOT_BIN)
#endif

namespace cepgen {
  namespace utils {
    /// Gnuplot drawable objects drawing utility
    class GnuplotDrawer final : public Drawer {
    public:
      /// Class constructor
      explicit GnuplotDrawer(const ParametersList& params)
          : Drawer(params),
            extension_(steer<std::string>("extension")),
            persist_(steer<bool>("persist")),
            size_(steer<std::vector<int> >("size")),
            font_(steer<std::string>("font")) {
        if (size_.size() != 2)
          throw CG_FATAL("GnuplotDrawer") << "Invalid canvas size specified: " << size_ << ".";
      }

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("Gnuplot drawing utility");
        desc.add<std::string>("extension", "png");
        desc.add<bool>("persist", false);
        desc.add<std::vector<int> >("size", {640, 480});
        desc.add<std::string>("font", "Helvetica,15");
        return desc;
      }

      const GnuplotDrawer& draw(const Graph1D&, const Mode&) const override;
      const GnuplotDrawer& draw(const Graph2D&, const Mode&) const override;
      const GnuplotDrawer& draw(const Hist1D&, const Mode&) const override;
      const GnuplotDrawer& draw(const Hist2D&, const Mode&) const override;

      const GnuplotDrawer& draw(const DrawableColl&,
                                const std::string& name = "",
                                const std::string& title = "",
                                const Mode& mode = Mode::none) const override;

    private:
      void execute(const Piper::Commands& cmds, const std::string& name) const {
        std::string term;
        if (extension_ == "pdf")
          term = "pdfcairo enhanced";
        else if (extension_ == "png")
          term = "pngcairo transparent enhanced";
        else if (extension_ == "tex")
          term = "epslatex";
        else if (extension_ == "ps")
          term = "postscript nobackground enhanced";
        else
          throw CG_FATAL("GnuplotDrawer:execute") << "Invalid extension set: '" + extension_ + "'";
        if (!font_.empty())
          term += " font '" + font_ + "'";
        term += " size " + merge(size_, ",");
        Piper::Commands full_cmds{"set term " + term, "set output '" + name + "." + extension_ + "'"};
        full_cmds += cmds;
        full_cmds += Piper::Commands{//"set key left bottom",
                                     //"set key right top",
                                     "exit"};
        Piper(std::string(GNUPLOT) + (persist_ ? " -persist" : "")).execute(full_cmds);
        CG_DEBUG("GnuplotDrawer:execute") << "Gnuplot just plotted:\n" << full_cmds;
      }

      static Piper::Commands preDraw(const Drawable& dr, const Mode& mode) {
        Piper::Commands cmds;
        if (!dr.title().empty())
          cmds += "set title '" + dr.title() + "'";
        //if (filling_)
        if (mode & Mode::grid)
          cmds += "set grid x y mx my";
        if (mode & Mode::logx)
          cmds += "set logscale x";
        if (mode & Mode::logy)
          cmds += "set logscale y";
        if (mode & Mode::logz)
          cmds += "set logscale z";
        for (const auto& ai : std::unordered_map<std::string, const Drawable::AxisInfo&>{
                 {"x", dr.xAxis()}, {"y", dr.yAxis()}, {"z", dr.zAxis()}}) {
          if (!ai.second.label().empty())
            cmds += "set " + ai.first + "label '" + ai.second.label() + "'";
          const auto& rng = ai.second.range();
          if (rng.valid())
            cmds += "set " + ai.first + "range [" + std::to_string(rng.min()) + ":" + std::to_string(rng.max()) + "]";
        }
        return cmds;
      }
      const std::string extension_;
      const bool persist_;
      const std::vector<int> size_;
      const std::string font_;
    };

    const GnuplotDrawer& GnuplotDrawer::draw(const Graph1D& graph, const Mode& mode) const {
      auto cmds = preDraw(graph, mode);
      cmds += "$DATA << EOD";
      for (const auto& pt : graph.points())
        cmds +=
            merge(std::vector<double>{pt.first.value, pt.first.value_unc, pt.second.value, pt.second.value_unc}, "\t");
      cmds += "EOD";
      cmds += "plot '$DATA' u 1:3 notitle";
      execute(cmds, graph.name());
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const Graph2D& graph, const Mode& mode) const {
      auto cmds = preDraw(graph, mode);
      cmds += "$DATA << EOD";
      const auto xs = graph.xCoords(), ys = graph.yCoords();
      const std::vector<double> xvec(xs.begin(), xs.end()), yvec(ys.begin(), ys.end());
      cmds += std::to_string(yvec.size()) + "\t" + merge(xvec, "\t");
      for (const auto& y : yvec) {
        std::ostringstream os;
        os << y;
        for (const auto& x : xvec)
          os << "\t" << graph.valueAt(x, y).value;
        cmds += os.str();
      }
      cmds += "EOD";
      cmds += "set autoscale xfix";
      cmds += "set autoscale yfix";
      cmds += "set autoscale cbfix";
      cmds += "plot '$DATA' matrix nonuniform with image notitle";
      execute(cmds, graph.name());
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const Hist1D& hist, const Mode& mode) const {
      auto cmds = preDraw(hist, mode);
      cmds += "set style data histograms";
      cmds += "set style histogram gap 0.";
      //cmds += "set style fill solid 1.0 border -1";
      cmds += "set style fill transparent pattern 2 bo";

      cmds += "$DATA << EOH";
      for (size_t ibin = 0; ibin < hist.nbins(); ++ibin)
        cmds += merge(std::vector<double>{hist.binRange(ibin).x(0.5), hist.value(ibin)}, "\t");
      cmds += "EOH";
      cmds += "set style data lines";
      cmds += "set yrange [" + std::to_string(hist.minimum()) + ":" + std::to_string(hist.maximum()) + "]";
      cmds += "set xtics 1 norangelimit nomirror";
      cmds += "set style fill solid 0.5 noborder";
      cmds += "set jitter spread 0.5";
      cmds += "plot '$DATA' using 1:2 bins=" + std::to_string(hist.nbins()) + " with boxes notitle";
      execute(cmds, hist.name());
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const Hist2D&, const Mode&) const {
      CG_WARNING("GnuplotDrawer:draw") << "Not yet implemented";
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const DrawableColl&,
                                             const std::string& name,
                                             const std::string& title,
                                             const Mode& mode) const {
      CG_WARNING("GnuplotDrawer:draw") << "Not yet implemented";
      return *this;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("gnuplot", GnuplotDrawer)
