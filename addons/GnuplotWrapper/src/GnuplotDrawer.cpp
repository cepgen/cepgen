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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Piper.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

#ifndef GNUPLOT_BIN
#error "Gnuplot executable must be specified using GNUPLOT_BIN!"
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define GNUPLOT TOSTRING(GNUPLOT_BIN)
#endif

using namespace cepgen;
using namespace std::string_literals;

/// Gnuplot drawable objects drawing utility
class GnuplotDrawer final : public utils::Drawer {
public:
  /// Class constructor
  explicit GnuplotDrawer(const ParametersList& params)
      : Drawer(params),
        extension_(steer<std::string>("extension")),
        persist_(steer<bool>("persist")),
        size_(steer<std::vector<std::string> >("size")),
        font_(steer<std::string>("font")),
        plot_style_(steer<std::string>("plotStyle")) {
    if (size_.size() != 2)
      throw CG_FATAL("GnuplotDrawer") << "Invalid canvas size specified: " << size_ << ".";
  }

  static ParametersDescription description() {
    auto desc = Drawer::description();
    desc.setDescription("Gnuplot drawing utility");
    desc.add("extension", "png"s);
    desc.add("persist", false);
    desc.add("size", std::vector{"30cm"s, "20cm"s});
    desc.add("font", ""s);
    desc.add("plotStyle", "lp"s);
    return desc;
  }

  const GnuplotDrawer& draw(const utils::Graph1D&, const Mode&) const override;
  const GnuplotDrawer& draw(const utils::Graph2D&, const Mode&) const override;
  const GnuplotDrawer& draw(const utils::Hist1D&, const Mode&) const override;
  const GnuplotDrawer& draw(const utils::Hist2D&, const Mode&) const override;

  const GnuplotDrawer& draw(const utils::DrawableColl&,
                            const std::string& name = "",
                            const std::string& title = "",
                            const Mode& mode = Mode::none) const override;

private:
  void execute(const utils::Piper::Commands& cmds, const std::string& name) const {
    std::string term;
    if (extension_ == "pdf")
      term = "pdfcairo enhanced";
    else if (extension_ == "png")
      term = "pngcairo transparent enhanced";
    else if (extension_ == "tex")
      term = "epslatex";
    else if (extension_ == "ps")
      term = "postscript nobackground enhanced";
    else if (extension_ == "fig")
      term = "fig";
    else
      throw CG_FATAL("GnuplotDrawer:execute") << "Invalid extension set: '" + extension_ + "'";
    if (!font_.empty())
      term += " font '" + font_ + "'";
    term += " size " + utils::merge(size_, ",");
    utils::Piper::Commands full_cmds{"set term " + term, "set output '" + name + "." + extension_ + "'"};
    full_cmds += cmds;
    full_cmds += utils::Piper::Commands{//"set key left bottom",
                                        //"set key right top",
                                        "exit"};
    utils::Piper(std::string{GNUPLOT} + (persist_ ? " -persist" : "")).execute(full_cmds);
    CG_DEBUG("GnuplotDrawer:execute") << "Gnuplot just plotted:\n" << full_cmds;
  }

  static utils::Piper::Commands preDraw(const utils::Drawable& dr, const Mode& mode) {
    utils::Piper::Commands cmds;
    //if (filling_)
    if (mode & Mode::grid)
      cmds += "set grid x y mx my";
    if (mode & Mode::logx)
      cmds += "set logscale x";
    if (mode & Mode::logy)
      cmds += "set logscale y";
    if (mode & Mode::logz)
      cmds += "set logscale z";
    if (!dr.title().empty())
      cmds += "set title " + delatexify(dr.title());
    for (const auto& ai : std::unordered_map<std::string, const utils::Drawable::AxisInfo&>{
             {"x", dr.xAxis()}, {"y", dr.yAxis()}, {"z", dr.zAxis()}}) {
      if (!ai.second.label().empty())
        cmds += "set " + ai.first + "label " + delatexify(ai.second.label());
      const auto& rng = ai.second.range();
      if (rng.valid())
        cmds += "set " + ai.first + "range [" + std::to_string(rng.min()) + ":" + std::to_string(rng.max()) + "]";
    }
    cmds += "set label 'CepGen v" + version::tag + "' at graph 1,1.025 right";
    return cmds;
  }
  static utils::Piper::Commands drawGraph1D(const utils::Graph1D&, const Mode&, const std::string&);
  static utils::Piper::Commands drawHist1D(const utils::Hist1D&, const Mode&);
  static std::string delatexify(const std::string& tok) {
    //return "\""s + utils::replaceAll(tok, {{"\\", "\\\\"}}) + "\"";
    return "'"s + utils::replaceAll(tok, {{"'", "\\'"}}) + "'";
  }

  const std::string extension_;
  const bool persist_;
  const std::vector<std::string> size_;
  const std::string font_, plot_style_;
};

const GnuplotDrawer& GnuplotDrawer::draw(const utils::Graph1D& graph, const Mode& mode) const {
  auto cmds = preDraw(graph, mode);
  cmds += drawGraph1D(graph, mode, plot_style_);
  execute(cmds, graph.name());
  return *this;
}

const GnuplotDrawer& GnuplotDrawer::draw(const utils::Graph2D& graph, const Mode& mode) const {
  auto cmds = preDraw(graph, mode);
  cmds += "$DATA << EOD";
  const auto xs = graph.xCoords(), ys = graph.yCoords();
  const std::vector<double> x_vector(xs.begin(), xs.end()), y_vector(ys.begin(), ys.end());
  cmds += std::to_string(y_vector.size()) + "\t" + utils::merge(x_vector, "\t");
  for (const auto& y : y_vector) {
    std::ostringstream os;
    os << y;
    for (const auto& x : x_vector)
      os << "\t" << static_cast<double>(graph.valueAt(x, y));
    cmds += os.str();
  }
  cmds += "EOD";
  cmds += "set autoscale xfix";
  cmds += "set autoscale yfix";
  cmds += "set autoscale cbfix";
  if (mode & Mode::col) {
    cmds += "set hidden3d";
    cmds += "plot '$DATA' matrix nonuniform with image notitle";
  } else if (mode & Mode::cont) {
    cmds += "set view map";
    cmds += "set contour";
    cmds += "unset surface";
    cmds += "set isosamples 500,100";
    cmds += "set cntrlabel start 25 interval -1 font \",7\"";
    cmds += "splot '$DATA' matrix nonuniform with lines notitle";
  } else {
    cmds += "set hidden3d";
    cmds += "set style data lines";
    cmds += "unset contour";
    cmds += "splot '$DATA' matrix nonuniform notitle";
  }
  execute(cmds, graph.name());
  return *this;
}

const GnuplotDrawer& GnuplotDrawer::draw(const utils::Hist1D& hist, const Mode& mode) const {
  auto cmds = preDraw(hist, mode);
  cmds += drawHist1D(hist, mode);
  execute(cmds, hist.name());
  return *this;
}

const GnuplotDrawer& GnuplotDrawer::draw(const utils::Hist2D& hist, const Mode& mode) const {
  auto cmds = preDraw(hist, mode);
  cmds += "$DATA << EOD";
  {
    std::ostringstream os;
    os << hist.nbinsX();
    for (size_t ix = 0; ix < hist.nbinsX(); ++ix)
      os << "\t" << hist.binRangeX(ix).x(0.5);
    cmds += os.str();
  }
  for (size_t iy = 0; iy < hist.nbinsY(); ++iy) {
    std::ostringstream os;
    os << hist.binRangeY(iy).x(0.5);
    for (size_t ix = 0; ix < hist.nbinsX(); ++ix)
      os << "\t" << static_cast<double>(hist.value(ix, iy));
    cmds += os.str();
  }
  cmds += "EOD";
  if (mode & Mode::col) {
    cmds += "set hidden3d";
    cmds += "plot '$DATA' matrix nonuniform with image notitle";
  } else if (mode & Mode::cont) {
    cmds += "set view map";
    cmds += "set contour";
    cmds += "unset surface";
    cmds += "set isosamples 500,100";
    cmds += "splot '$DATA' matrix nonuniform with lines notitle";
  } else {
    cmds += "set hidden3d";
    cmds += "set style data lines";
    cmds += "unset contour";
    cmds += "splot '$DATA' matrix nonuniform notitle";
  }
  execute(cmds, hist.name());
  return *this;
}

const GnuplotDrawer& GnuplotDrawer::draw(const utils::DrawableColl& objs,
                                         const std::string& name,
                                         const std::string& title,
                                         const Mode& mode) const {
  if (objs.empty())
    return *this;
  auto cmds = preDraw(*objs.at(0), mode);
  cmds += "set title " + delatexify(title);
  std::vector<std::string> plot_cmds, splot_cmds;
  for (const auto* obj : objs) {
    if (obj->isGraph1D()) {
      if (const auto* gr = dynamic_cast<const utils::Graph1D*>(obj); gr) {
        auto gr_cmds = drawGraph1D(*gr, mode, plot_style_);
        auto it = gr_cmds.begin();
        while (it != gr_cmds.end())
          if (utils::startsWith(*it, "plot")) {
            plot_cmds.emplace_back(utils::replaceAll(it->substr(5), " notitle", " title " + delatexify(obj->title())));
            it = gr_cmds.erase(it);
          } else if (utils::startsWith(*it, "splot")) {
            splot_cmds.emplace_back(utils::replaceAll(it->substr(6), " notitle", " title " + delatexify(obj->title())));
            it = gr_cmds.erase(it);
          } else
            ++it;
        if (plot_cmds.empty() && splot_cmds.empty())
          throw CG_FATAL("GnuplotDrawer:draw")
              << "No drawing command found for graph with name \"" << obj->name() << "\"!";
        cmds += gr_cmds;
      }
    } else if (obj->isHist1D()) {
      if (const auto* hist = dynamic_cast<const utils::Hist1D*>(obj); hist) {
        auto h_cmds = drawHist1D(*hist, mode);
        auto it = h_cmds.begin();
        while (it != h_cmds.end())
          if (utils::startsWith(*it, "plot")) {
            plot_cmds.emplace_back(utils::replaceAll(it->substr(5), " notitle", " title " + delatexify(obj->title())));
            it = h_cmds.erase(it);
          } else if (utils::startsWith(*it, "splot")) {
            splot_cmds.emplace_back(utils::replaceAll(it->substr(6), " notitle", " title " + delatexify(obj->title())));
            it = h_cmds.erase(it);
          } else
            ++it;
        if (plot_cmds.empty() && splot_cmds.empty())
          throw CG_FATAL("GnuplotDrawer:draw")
              << "No drawing command found for histogram with name \"" << obj->name() << "\"!";
        cmds += h_cmds;
      }
    }
  }
  if (plot_cmds.empty() && splot_cmds.empty())
    throw CG_FATAL("GnuplotDrawer:draw") << "No drawing command found!";
  if (!plot_cmds.empty() && !splot_cmds.empty())
    throw CG_FATAL("GnuplotDrawer:draw") << "Cannot combine 'flat', and surface-like drawing commands!";
  if (!plot_cmds.empty())
    cmds += "plot " + utils::merge(plot_cmds, ", ");
  else if (!splot_cmds.empty())
    cmds += "splot " + utils::merge(splot_cmds, ", ");
  execute(cmds, name);
  return *this;
}

utils::Piper::Commands GnuplotDrawer::drawGraph1D(const utils::Graph1D& graph, const Mode&, const std::string& style) {
  utils::Piper::Commands cmds;
  auto random_filename = utils::randomString(5);
  cmds += "$DATA_" + random_filename + " << EOD";
  for (const auto& pt : graph.points())
    cmds +=
        utils::merge(std::vector<double>{pt.first.value, pt.first.value_unc, pt.second, pt.second.uncertainty()}, "\t");
  cmds += "EOD";
  cmds += "plot '$DATA_" + random_filename + "' u 1:3 w " + style + " notitle";
  return cmds;
}

utils::Piper::Commands GnuplotDrawer::drawHist1D(const utils::Hist1D& hist, const Mode&) {
  utils::Piper::Commands cmds;
  auto random_filename = utils::randomString(5);
  cmds += "set style data histograms";
  cmds += "set style histogram gap 0.";
  //cmds += "set style fill solid 1.0 border -1";
  cmds += "set style fill transparent pattern 2 bo";

  cmds += "$DATA_" + random_filename + " << EOH";
  for (size_t ibin = 0; ibin < hist.nbins(); ++ibin)
    cmds += utils::merge(std::vector{hist.binRange(ibin).x(0.5), static_cast<double>(hist.value(ibin))}, "\t");
  cmds += "EOH";
  cmds += "set style data lines";
  cmds += "set yrange [" + std::to_string(hist.minimum()) + ":" + std::to_string(hist.maximum()) + "]";
  cmds += "set xtics 1 norangelimit nomirror";
  cmds += "set style fill solid 0.5 noborder";
  cmds += "set jitter spread 0.5";
  cmds += "plot '$DATA_" + random_filename + "' using 1:2 bins=" + std::to_string(hist.nbins()) + " with boxes notitle";
  return cmds;
}
REGISTER_DRAWER("gnuplot", GnuplotDrawer);
