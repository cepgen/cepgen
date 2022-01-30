#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"

#ifndef TD_BIN
#error "Topdrawer executable must be specified using TD_BIN!"
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define TD TOSTRING(TD_BIN)
#endif

namespace cepgen {
  namespace utils {
    class DrawerTopdrawer : public Drawer {
    public:
      explicit DrawerTopdrawer(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("Topdrawer plotter");
        return desc;
      }

      const DrawerTopdrawer& draw(const Graph1D&, const Mode&) const override;
      const DrawerTopdrawer& draw(const Graph2D&, const Mode&) const override;
      const DrawerTopdrawer& draw(const Hist1D&, const Mode&) const override;
      const DrawerTopdrawer& draw(const Hist2D&, const Mode&) const override;

      const DrawerTopdrawer& draw(const DrawableColl&,
                                  const std::string& name = "",
                                  const std::string& title = "",
                                  const Mode& mode = Mode::none) const override;

    private:
      struct Commands : std::vector<std::string> {
        Commands& operator+=(const Commands& oth) {
          std::copy(oth.begin(), oth.end(), std::back_inserter(*this));
          return *this;
        }
        Commands& operator+=(const std::string& str) {
          emplace_back(str);
          return *this;
        }
        friend std::ostream& operator<<(std::ostream& os, const Commands& cmds) {
          std::string sep;
          os << "{";
          for (const auto& cmd : cmds)
            os << sep << cmd, sep = "\n";
          return os << "}";
        }
      };
      static void execute(const Commands&, const std::string&);
      static Commands plot(const Graph1D&);
      static Commands plot(const Hist1D&);
      static Commands preDraw(const Drawable&, const Mode&);
      static Commands postDraw(const Drawable&, const Mode&);
    };

    DrawerTopdrawer::DrawerTopdrawer(const ParametersList& params) : Drawer(params) {}

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Graph1D& graph, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(graph, mode);
      cmds += plot(graph);
      cmds += postDraw(graph, mode);
      execute(cmds, graph.name());
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Graph2D& graph, const Mode& mode) const {
      CG_WARNING("DrawerTopdrawer:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Hist1D& hist, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(hist, mode);
      cmds += plot(hist);
      cmds += postDraw(hist, mode);
      execute(cmds, hist.name());
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Hist2D&, const Mode&) const {
      CG_WARNING("DrawerTopdrawer:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const DrawableColl& objs,
                                                 const std::string& name,
                                                 const std::string& title,
                                                 const Mode& mode) const {
      std::vector<std::string> line_styles = {
          "SOLID", "DOTS", "DASHES", "DAASHES", "DOTDASH", "SPACE", "PATTERNED", "FUNNY", "PERMANENT"};
      size_t plot_id = 0;
      Commands cmds;
      const Drawable* first{nullptr};
      Commands cmds_plots;
      for (const auto* obj : objs) {
        if (obj->isGraph1D()) {
          const auto* gr = dynamic_cast<const Graph1D*>(obj);
          cmds_plots.emplace_back("SET TEXTURE " + line_styles.at(plot_id++));
          cmds_plots += plot(*gr);
          if (!first)
            first = gr;
        } else if (obj->isHist1D()) {
          const auto* hist = dynamic_cast<const Hist1D*>(obj);
          cmds_plots.emplace_back("SET TEXTURE " + line_styles.at(plot_id++));
          cmds_plots += plot(*hist);
          if (!first)
            first = hist;
        }
      }
      cmds += preDraw(*first, mode);
      cmds += cmds_plots;
      cmds += postDraw(*first, mode);
      execute(cmds, name);
      return *this;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::plot(const Graph1D& graph) {
      Commands cmds;
      for (const auto& pt : graph.points())
        cmds += std::to_string(pt.first.value) + "," + std::to_string(pt.second.value) + "," +
                std::to_string(pt.first.value_unc) + "," + std::to_string(pt.second.value_unc);
      cmds += "JOIN";
      return cmds;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::plot(const Hist1D& hist) {
      Commands cmds;
      for (size_t i = 0; i < hist.nbins(); ++i) {
        const auto& bin = hist.binRange(i);
        cmds += std::to_string(bin.x(0.5)) + "," + std::to_string(hist.value(i)) + "," +
                std::to_string(0.5 * bin.range()) + "," + std::to_string(hist.valueUnc(i));
      }
      cmds += "HIST";
      return cmds;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::preDraw(const Drawable& dr, const Mode& mode) {
      Commands cmds;
      cmds += "SET DEVICE POSTSCR ORIENTATION 3";
      cmds += "SET FONT DUPLEX";
      cmds += "SET FILL FULL";
      if (mode & Mode::grid)
        cmds += "SET GRID ON WIDTH=1 DOTS";
      if (mode & Mode::logx)
        cmds += "SET SCALE X LOG";
      if (mode & Mode::logy)
        cmds += "SET SCALE Y LOG";
      const auto& xrng = dr.xAxis().range();
      if (xrng.valid())
        cmds += "SET LIMITS X " + std::to_string(xrng.min()) + " TO " + std::to_string(xrng.max());
      const auto& yrng = dr.yAxis().range();
      if (yrng.valid())
        cmds += "SET LIMITS Y " + std::to_string(yrng.min()) + " TO " + std::to_string(yrng.max());
      return cmds;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::postDraw(const Drawable& dr, const Mode& mode) {
      Commands cmds;
      cmds += "TITLE TOP '" + dr.title() + "'";
      cmds += "TITLE BOTTOM '" + dr.xAxis().label() + "'";
      cmds += "TITLE LEFT '" + dr.yAxis().label() + "'";
      cmds += "EXIT";
      return cmds;
    }

    void DrawerTopdrawer::execute(const Commands& cmds, const std::string& name) {
      const auto cmd = std::string("TOPDRAWER_OUTPUT=") + name + ".ps " + TD;
      std::unique_ptr<FILE, decltype(&pclose)> file(popen(cmd.c_str(), "w"), pclose);
      for (const auto& cmd : cmds)
        fputs((cmd + "\n").c_str(), file.get());
      CG_DEBUG("DrawerTopdrawer:execute") << "Topdrawer just plotted:\n" << cmds;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("topdrawer", DrawerTopdrawer)
