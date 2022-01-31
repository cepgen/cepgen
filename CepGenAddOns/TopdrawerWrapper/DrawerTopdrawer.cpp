#include <cmath>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

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
      static Commands preDraw(const Drawable&, const Mode&);
      static Commands plot(const Graph1D&);
      static Commands plot(const Graph2D&, const Mode&);
      static Commands plot(const Hist1D&);
      static Commands plot(const Hist2D&, const Mode&);
      static Commands postDraw(const Drawable&, const Mode&);
      static Commands stringify(const std::string&, const std::string&);
      static const std::map<std::string, std::pair<char, char> > kSpecChars;
    };

    const std::map<std::string, std::pair<char, char> > DrawerTopdrawer::kSpecChars = {
        {"Alpha", {'A', 'F'}},      {"Beta", {'B', 'F'}},
        {"Chi", {'C', 'F'}},        {"Delta", {'D', 'F'}},
        {"Epsilon", {'E', 'F'}},    {"Phi", {'F', 'F'}},
        {"Gamma", {'G', 'F'}},      {"Eta", {'H', 'F'}},
        {"Iota", {'I', 'F'}},       {"Kappa", {'K', 'F'}},
        {"Lambda", {'L', 'F'}},     {"Mu", {'M', 'F'}},
        {"Nu", {'N', 'F'}},         {"Omicron", {'O', 'F'}},
        {"Pi", {'P', 'F'}},         {"Theta", {'Q', 'F'}},
        {"Rho", {'R', 'F'}},        {"Sigma", {'S', 'F'}},
        {"Tau", {'T', 'F'}},        {"Upsilon", {'U', 'F'}},
        {"Omega", {'W', 'F'}},      {"Xi", {'X', 'F'}},
        {"Psi", {'Y', 'F'}},        {"Zeta", {'Z', 'F'}},
        {"alpha", {'A', 'G'}},      {"beta", {'B', 'G'}},
        {"chi", {'C', 'G'}},        {"delta", {'D', 'G'}},
        {"epsilon", {'E', 'G'}},    {"phi", {'G', 'G'}},
        {"gamma", {'G', 'G'}},      {"eta", {'H', 'G'}},
        {"iota", {'I', 'G'}},       {"kappa", {'K', 'G'}},
        {"lambda", {'L', 'G'}},     {"mu", {'M', 'G'}},
        {"nu", {'N', 'G'}},         {"omicron", {'O', 'G'}},
        {"pi", {'P', 'G'}},         {"theta", {'Q', 'G'}},
        {"rho", {'R', 'G'}},        {"sigma", {'S', 'G'}},
        {"tau", {'T', 'G'}},        {"upsilon", {'U', 'G'}},
        {"omega", {'W', 'G'}},      {"xi", {'X', 'G'}},
        {"psi", {'Y', 'G'}},        {"zeta", {'Z', 'G'}},
        {"simeq", {'C', 'M'}},      {"gt", {'G', 'M'}},
        {"ge", {'H', 'M'}},         {"int", {'I', 'M'}},
        {"icirc", {'J', 'M'}},      {"lt", {'L', 'M'}},
        {"le", {'M', 'M'}},         {"neq", {'N', 'M'}},
        {"sim", {'S', 'M'}},        {"perp", {'T', 'M'}},
        {"dpar", {'Y', 'M'}},       {"infty", {'0', 'M'}},
        {"sqrt", {'2', 'M'}},       {"pm", {'+', 'M'}},
        {"mp", {'-', 'M'}},         {"otimes", {'*', 'M'}},
        {"equiv", {'=', 'M'}},      {"cdot", {'.', 'M'}},
        {"times", {'1', 'O'}},      {"leftarrow", {'L', 'W'}},
        {"rightarrow", {'R', 'W'}}, {"leftrightarrow", {'B', 'W'}},
        {"langle", {'B', 'S'}},     {"rangle", {'E', 'S'}},
        {"hbar", {'H', 'K'}},       {"lambdabar", {'L', 'K'}}};

    DrawerTopdrawer::DrawerTopdrawer(const ParametersList& params) : Drawer(params) {}

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Graph1D& graph, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(graph, mode);
      cmds += plot(graph);
      cmds += stringify("TITLE TOP", graph.title());
      cmds += postDraw(graph, mode);
      execute(cmds, graph.name());
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Graph2D& graph, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(graph, mode);
      cmds += plot(graph, mode);
      cmds += stringify("TITLE TOP", graph.title());
      cmds += postDraw(graph, mode);
      execute(cmds, graph.name());
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Hist1D& hist, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(hist, mode);
      cmds += plot(hist);
      cmds += stringify("TITLE TOP", hist.title());
      cmds += postDraw(hist, mode);
      execute(cmds, hist.name());
      return *this;
    }

    const DrawerTopdrawer& DrawerTopdrawer::draw(const Hist2D& hist, const Mode& mode) const {
      Commands cmds;
      cmds += preDraw(hist, mode);
      cmds += plot(hist, mode);
      cmds += stringify("TITLE TOP", hist.title());
      cmds += postDraw(hist, mode);
      execute(cmds, hist.name());
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
      cmds += stringify("TITLE TOP", title);
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

    DrawerTopdrawer::Commands DrawerTopdrawer::plot(const Graph2D& graph, const Mode& mode) {
      Commands cmds;
      auto to_fortran_float = [](double val) -> std::string {
        return utils::replace_all(utils::format("%g", val), {{"e", "D"}});
      };
      cmds += "READ MESH";
      std::ostringstream osl;
      for (const auto& yval : graph.points().begin()->second)
        osl << " " << to_fortran_float(fabs(yval.first.value) < 1.e-14 ? 0. : yval.first.value);
      cmds += "Y" + osl.str();
      for (const auto& xval : graph.points()) {
        osl.str("");
        osl << "X " << to_fortran_float(xval.first.value) << " Z";
        for (const auto& yval : xval.second)
          osl << " " << (std::isfinite(yval.second.value) ? to_fortran_float(yval.second.value) : "0.");
        cmds += osl.str();
      }
      if (mode & Mode::col)
        cmds += "JOIN";
      else if (mode & Mode::cont)
        cmds += "CONTOUR";
      else {
        cmds += "SET THREE OFF";
        cmds += "PLOT";
      }
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

    DrawerTopdrawer::Commands DrawerTopdrawer::plot(const Hist2D& hist, const Mode& mode) {
      Commands cmds;
      cmds += "READ MESH BINS";
      std::ostringstream osl;
      std::string sep;
      for (size_t iy = 0; iy < hist.nbinsY(); ++iy)
        osl << sep << hist.binRangeY(iy).min(), sep = " ";
      osl << " " << hist.binRangeY(hist.nbinsY() - 1).max();
      cmds += "FOR Y=" + osl.str();
      for (size_t ix = 0; ix < hist.nbinsX(); ++ix) {
        osl.str("");
        osl << "X=" << hist.binRangeX(ix).x(0.5) << " Z=";
        for (size_t iy = 0; iy < hist.nbinsY(); ++iy) {
          osl << " " << hist.value(ix, iy);
        }
        cmds += osl.str();
      }
      if (mode & Mode::col)
        cmds += "JOIN";
      else if (mode & Mode::cont)
        cmds += "CONTOUR";
      else {
        cmds += "SET THREE OFF";
        cmds += "PLOT";
      }
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
      if (mode & Mode::logz)
        cmds += "SET SCALE Z LOG";
      const auto& xrng = dr.xAxis().range();
      if (xrng.valid())
        cmds += "SET LIMITS X " + std::to_string(xrng.min()) + " TO " + std::to_string(xrng.max());
      const auto& yrng = dr.yAxis().range();
      if (yrng.valid())
        cmds += "SET LIMITS Y " + std::to_string(yrng.min()) + " TO " + std::to_string(yrng.max());
      const auto& zrng = dr.zAxis().range();
      if (zrng.valid())
        cmds += "SET LIMITS Z " + std::to_string(zrng.min()) + " TO " + std::to_string(zrng.max());
      return cmds;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::postDraw(const Drawable& dr, const Mode&) {
      Commands cmds;
      cmds += stringify("TITLE BOTTOM", dr.xAxis().label());
      cmds += stringify("TITLE LEFT", dr.yAxis().label());
      return cmds;
    }

    void DrawerTopdrawer::execute(const Commands& cmds, const std::string& name) {
      const auto cmd = std::string("TOPDRAWER_OUTPUT=") + name + ".ps " + TD;
      std::unique_ptr<FILE, decltype(&pclose)> file(popen(cmd.c_str(), "w"), pclose);
      for (const auto& cmd : cmds)
        fputs((cmd + "\n").c_str(), file.get());
      fputs("EXIT", file.get());
      CG_DEBUG("DrawerTopdrawer:execute") << "Topdrawer just plotted:\n" << cmds;
    }

    DrawerTopdrawer::Commands DrawerTopdrawer::stringify(const std::string& label, const std::string& str) {
      bool in_math{false}, in_bs{false};
      std::map<int, std::string> m_spec_char;
      std::string lab, mod;
      for (size_t i = 0; i < str.size(); ++i) {
        const auto ch = str[i];
        if (ch == '$' && (i == 0 || str[i - 1] != '\\')) {
          in_math = !in_math;
          continue;
        }
        if (ch == '\\') {
          in_bs = true;
          m_spec_char[lab.size()] = "";
          lab.push_back('*');
          continue;
        }
        if (in_bs) {
          if (ch == ' ' || ch == '_' || ch == '(' || ch == ')' || ch == '{' || ch == '}' || ch == '[' || ch == ']')
            in_bs = false;
          else if (ch == '\\') {
            m_spec_char[lab.size()] = "";
            lab.push_back('*');
            continue;
          } else {
            m_spec_char.rbegin()->second.push_back(ch);
            continue;
          }
        }
        if (ch != '_' && ch != '{' && ch != '}')
          lab.push_back(ch);
      }
      mod = std::string(lab.size(), ' ');
      CG_LOG << lab << "\n\n" << m_spec_char;
      for (const auto& ch : m_spec_char) {
        const auto& tok = kSpecChars.at(ch.second);
        lab[ch.first] = tok.first;
        mod[ch.first] = tok.second;
      }
      Commands out;
      out += label + " '" + lab + "'";
      out += "CASE" + std::string(label.size() - 4, ' ') + " '" + mod + "'";
      return out;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("topdrawer", DrawerTopdrawer)
