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
    /**
     * This object allows to invoke gnuplot, the portable command-line driven
     * graphing utility for Linux, OS/2, MS Windows, OSX, VMS, and many other
     * platforms.
     * @brief Plotting utility used in control plots generation
     */
    class GnuplotDrawer final : public Drawer {
    public:
      /// Class constructor
      explicit GnuplotDrawer(const ParametersList& params)
          : Drawer(params), extension_(steer<std::string>("extension")), persist_(steer<bool>("persist")) {}

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.add<std::string>("extension", "png");
        desc.add<bool>("persist", false);
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
        Piper piper(std::string(GNUPLOT) + (persist_ ? " -persist" : ""));
        if (extension_ == "png")
          piper.execute({{"set term pngcairo transparent enhanced font 'arial,10' fontscale 1.0 size 800, 600"},
                         {"set output '" + name + "." + extension_ + "'"}});
        else if (extension_ == "tex")
          piper.execute({"set term epslatex"});
        else
          throw CG_FATAL("GnuplotDrawer:execute") << "Invalid extension set: '" + extension_ + "'";
        piper.execute({//{"set key left bottom"},
                       {"set key right top"}});
        CG_DEBUG("GnuplotDrawer:execute") << "Gnuplot just plotted:\n" << cmds;
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
          cmds += "set " + ai.first + "label '" + ai.second.label() + "'";
          const auto& rng = ai.second.range();
          if (rng.valid())
            cmds += "set " + ai.first + "range [" + std::to_string(rng.min()) + ":" + std::to_string(rng.max()) + "]";
        }
        return cmds;
      }
      const std::string extension_;
      const bool persist_;
    };

    const GnuplotDrawer& GnuplotDrawer::draw(const Graph1D&, const Mode&) const {
      CG_WARNING("GnuplotDrawer:draw") << "Not yet implemented";
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const Graph2D&, const Mode&) const {
      CG_WARNING("GnuplotDrawer:draw") << "Not yet implemented";
      return *this;
    }

    const GnuplotDrawer& GnuplotDrawer::draw(const Hist1D& hist, const Mode& mode) const {
      Piper::Commands cmds;
      cmds += "set style data histograms";
      cmds += "set style histogram gap 0.";
      //cmds += "set style fill solid 1.0 border -1";
      cmds += "set style fill transparent pattern 2 bo";

      const std::string tmp_filename = "/tmp/hist.gp.dat";  // + randomString(10);
      {
        CG_DEBUG("GnuplotDrawer:draw") << "Temporary file for histogram data: '" << tmp_filename << "'.";
        std::ofstream tmp(tmp_filename);
        for (size_t ibin = 0; ibin < hist.nbins(); ++ibin) {
          const auto& rng = hist.binRange(ibin);
          tmp << rng.min() << "\t" << rng.max() << "\t" << hist.value(ibin) << "\t" << hist.valueUnc(ibin) << "\n";
        }
      }
      cmds += "set xtics auto";
      cmds += "everyNth(lab,N) =((int(column(0)) % N == 0) ? stringcolumn(lab) : \"\"); ";
      cmds += "plot '" + tmp_filename + "' using 2:xtic(everyNth(1, " + std::to_string(hist.nbins()) +
              ")) w histeps t \"" + hist.title() + "\"";
      execute(cmds, hist.name());
      //remove(tmp_filename);
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
