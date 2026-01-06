/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2026  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Piper.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

#ifndef TD_BIN
#error "Topdrawer executable must be specified using TD_BIN!"
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define TD TOSTRING(TD_BIN)
#endif

using namespace std::string_literals;

namespace cepgen::topdrawer {
  class Drawer : public utils::Drawer {
  public:
    explicit Drawer(const ParametersList& params)
        : utils::Drawer(params), font_(utils::toUpper(steer<std::string>("font"))), filling_(steer<bool>("filling")) {}

    static ParametersDescription description() {
      auto desc = utils::Drawer::description();
      desc.setDescription("Topdrawer plotter");
      desc.add("font", "duplex"s).setDescription("Topdrawer font to use");
      desc.add("filling", true).setDescription("allow to fill the whole available space?");
      return desc;
    }

    const Drawer& draw(const utils::Graph1D& graph, const Mode& mode) const override {
      utils::Piper::Commands cmds;
      cmds += preDraw(graph, mode);
      cmds += plot(graph);
      cmds += stringify("TITLE TOP", graph.title());
      cmds += postDraw(graph, mode);
      execute(cmds, graph.name());
      return *this;
    }
    const Drawer& draw(const utils::Graph2D& graph, const Mode& mode) const override {
      utils::Piper::Commands cmds;
      cmds += preDraw(graph, mode);
      cmds += plot(graph, mode);
      cmds += stringify("TITLE TOP", graph.title());
      cmds += postDraw(graph, mode);
      execute(cmds, graph.name());
      return *this;
    }
    const Drawer& draw(const utils::Hist1D& hist, const Mode& mode) const override {
      utils::Piper::Commands cmds;
      cmds += preDraw(hist, mode);
      cmds += plot(hist);
      cmds += stringify("TITLE TOP", hist.title());
      cmds += postDraw(hist, mode);
      execute(cmds, hist.name());
      return *this;
    }
    const Drawer& draw(const utils::Hist2D& hist, const Mode& mode) const override {
      utils::Piper::Commands cmds;
      cmds += preDraw(hist, mode);
      cmds += plot(hist, mode);
      cmds += stringify("TITLE TOP", hist.title());
      cmds += postDraw(hist, mode);
      execute(cmds, hist.name());
      return *this;
    }

    const Drawer& draw(const utils::DrawableColl&,
                       const std::string& name = "",
                       const std::string& title = "",
                       const Mode& mode = Mode::none) const override;

  private:
    static void execute(const utils::Piper::Commands&, const std::string&);
    static utils::Piper::Commands plot(const utils::Graph1D&);
    static utils::Piper::Commands plot(const utils::Graph2D&, const Mode&);
    static utils::Piper::Commands plot(const utils::Hist1D&);
    static utils::Piper::Commands plot(const utils::Hist2D&, const Mode&);
    static utils::Piper::Commands postDraw(const utils::Drawable&, const Mode&);
    static utils::Piper::Commands stringify(const std::string&, const std::string&);
    static const std::unordered_map<std::string, std::pair<char, char> > kSpecialCharacters;

    utils::Piper::Commands preDraw(const utils::Drawable&, const Mode&) const;

    const std::string font_;
    const bool filling_;
  };

  const std::unordered_map<std::string, std::pair<char, char> > Drawer::kSpecialCharacters = {
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

  const Drawer& Drawer::draw(const utils::DrawableColl& objs,
                             const std::string& name,
                             const std::string& title,
                             const Mode& mode) const {
    std::vector<std::string> line_styles = {
        "SOLID", "DOTS", "DASHES", "DAASHES", "DOTDASH", "SPACE", "PATTERNED", "FUNNY", "PERMANENT"};
    size_t plot_id = 0;
    const utils::Drawable* first{nullptr};
    utils::Piper::Commands cmds_plots;
    for (const auto* obj : objs) {
      const auto line_style = (plot_id++) % line_styles.size();
      if (obj->isGraph1D()) {
        const auto* gr = dynamic_cast<const utils::Graph1D*>(obj);
        cmds_plots.emplace_back("SET TEXTURE " + line_styles.at(line_style));
        cmds_plots += plot(*gr);
        if (!first)
          first = gr;
      } else if (obj->isHist1D()) {
        const auto* hist = dynamic_cast<const utils::Hist1D*>(obj);
        cmds_plots.emplace_back("SET TEXTURE " + line_styles.at(line_style));
        cmds_plots += plot(*hist);
        if (!first)
          first = hist;
      } else
        throw CG_FATAL("topdrawer:Drawer:draw") << "Invalid object type to be plotted in multigraph!";
    }
    if (!first)
      throw CG_FATAL("topdrawer:Drawer:draw") << "No object defined as the first drawable in the canvas.";
    utils::Piper::Commands cmds;
    cmds += preDraw(*first, mode);
    cmds += cmds_plots;
    cmds += postDraw(*first, mode);
    cmds += stringify("TITLE TOP", title);
    execute(cmds, name);
    return *this;
  }

  utils::Piper::Commands Drawer::plot(const utils::Graph1D& graph) {
    utils::Piper::Commands cmds;
    for (const auto& pt : graph.points())
      cmds += utils::format("%g,%g,%g,%g", pt.first.value, pt.second, pt.first.value_unc, pt.second.uncertainty());
    cmds += "JOIN";
    return cmds;
  }

  utils::Piper::Commands Drawer::plot(const utils::Graph2D& graph, const Mode& mode) {
    utils::Piper::Commands cmds;
    auto to_fortran_float = [](double val) -> std::string {
      return utils::replaceAll(utils::format("%g", val), {{"e", "D"}});
    };
    cmds += "READ MESH";
    std::ostringstream osl;
    for (const auto& yval : graph.yCoords())
      osl << " " << to_fortran_float(fabs(yval) < 1.e-14 ? 0. : yval);
    cmds += "Y" + osl.str();
    for (const auto& xval : graph.points()) {
      osl.str("");
      osl << "X " << to_fortran_float(xval.first.value) << " Z";
      for (const auto& yval : xval.second)
        osl << " " << (std::isfinite(yval.second) ? to_fortran_float(yval.second) : "0.");
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

  utils::Piper::Commands Drawer::plot(const utils::Hist1D& hist) {
    utils::Piper::Commands cmds;
    for (size_t i = 0; i < hist.nbins(); ++i) {
      const auto& bin = hist.binRange(i);
      const auto& val = hist.value(i);
      cmds += utils::format("%g,%g,%g,%g", bin.x(0.5), static_cast<double>(val), 0.5 * bin.range(), val.uncertainty());
    }
    cmds += "HISTOGRAM";
    return cmds;
  }

  utils::Piper::Commands Drawer::plot(const utils::Hist2D& hist, const Mode& mode) {
    utils::Piper::Commands cmds;
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
        osl << " " << static_cast<double>(hist.value(ix, iy));
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

  utils::Piper::Commands Drawer::preDraw(const utils::Drawable& dr, const Mode& mode) const {
    utils::Piper::Commands cmds;
    cmds += "SET DEVICE POSTSCR ORIENTATION 3";
    cmds += "SET FONT " + font_;
    if (filling_)
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
      cmds += utils::format("SET LIMITS X %g TO %g", xrng.min(), xrng.max());
    const auto& yrng = dr.yAxis().range();
    if (yrng.valid())
      cmds += utils::format("SET LIMITS Y %g TO %g", yrng.min(), yrng.max());
    const auto& zrng = dr.zAxis().range();
    if (zrng.valid())
      cmds += utils::format("SET LIMITS Z %g TO %g", zrng.min(), zrng.max());
    return cmds;
  }

  utils::Piper::Commands Drawer::postDraw(const utils::Drawable& dr, const Mode&) {
    utils::Piper::Commands cmds;
    cmds += stringify("TITLE BOTTOM", dr.xAxis().label());
    cmds += stringify("TITLE LEFT", dr.yAxis().label());
    cmds += stringify("TITLE CENTER 10.8 9.25", "CepGen v" + version::tag);
    return cmds;
  }

  void Drawer::execute(const utils::Piper::Commands& cmds, const std::string& name) {
    utils::Piper("TOPDRAWER_OUTPUT=" + name + ".ps " + TD).execute(cmds).execute({"EXIT"});
    CG_DEBUG("topdrawer:Drawer:execute") << "Topdrawer just plotted:\n" << cmds;
  }

  utils::Piper::Commands Drawer::stringify(const std::string& label, const std::string& str) {
    bool in_math{false}, in_bs{false}, in_sub{false}, in_sup{false};
    std::map<int, std::string> m_spec_char, m_sub_char;
    std::string lab;
    auto str_parsed = utils::parseSpecialChars(str);
    for (size_t i = 0; i < str_parsed.size(); ++i) {
      const auto ch = str_parsed[i];
      if (ch == '$' && (i == 0 || str_parsed[i - 1] != '\\')) {
        in_math = !in_math;
        continue;
      }
      if (ch == '_') {  // check if we are in subscript mode
        in_sub = true;
        m_sub_char[lab.size()] = "";
        continue;
      }
      if (ch == '^') {  // check if we are in superscript mode
        in_sup = true;
        m_sub_char[lab.size()] = "";
        continue;
      }
      if (in_sub || in_sup) {
        if (ch == '{') {
          lab.push_back(in_sup ? '0' : '2');
          continue;
        }
        if (ch == '}') {
          lab.push_back(in_sup ? '1' : '3');
          if (in_sub)
            in_sub = false;
          if (in_sup)
            in_sup = false;
          continue;
        }
        m_sub_char.rbegin()->second.push_back(ch);
        lab.push_back(ch);
        continue;
      }
      if (ch == '\\') {  // check if we have a special character
        in_bs = true;
        m_spec_char[lab.size()] = "";
        lab.push_back('*');
        continue;
      }
      if (in_bs) {
        if (ch == ' ' || ch == '_' || ch == '/' || ch == '(' || ch == ')' || ch == '{' || ch == '}' || ch == '[' ||
            ch == ']')
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
      lab.push_back(ch);  // otherwise assume we are just pushing into the characters buffer
    }
    std::string mod(lab.size(), ' ');
    for (const auto& ch : m_spec_char)
      if (kSpecialCharacters.count(ch.second) > 0)
        std::tie(lab[ch.first], mod[ch.first]) = kSpecialCharacters.at(ch.second);
      else
        CG_WARNING("topdrawer:Drawer:stringify")
            << "Special character '" << ch.second << "' is not defined. Please either define it or use another one.";
    for (const auto& ch : m_sub_char) {
      mod[ch.first] = 'C';
      mod[ch.first + ch.second.size() + 1] = 'C';
    }
    utils::Piper::Commands out;
    out += label + " '" + lab + "'";
    out += "CASE" + std::string(label.size() - 4, ' ') + " '" + mod + "'";
    return out;
  }
}  // namespace cepgen::topdrawer
using TopdrawerDrawer = cepgen::topdrawer::Drawer;
REGISTER_DRAWER("topdrawer", TopdrawerDrawer);
