/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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
#include <sciplot/sciplot.hpp>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Piper.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

using namespace std::string_literals;

/// Collection of Sciplot utilities
namespace cepgen::sciplot {
  /// Sciplot drawable objects drawing utility
  class Drawer final : public utils::Drawer {
  public:
    /// Class constructor
    explicit Drawer(const ParametersList& params)
        : utils::Drawer(params),
          palette_name_(steer<std::string>("paletteName")),
          font_name_(steer<std::string>("fontName")),
          width_(steer<int>("width")),
          height_(steer<int>("height")),
          font_size_(steer<int>("fontSize")),
          line_width_(steer<int>("lineWidth")) {}

    static ParametersDescription description() {
      auto desc = utils::Drawer::description();
      desc.setDescription("Sciplot drawing utility");
      desc.add("paletteName", "set1"s);
      desc.add("fontName", "Palatino"s);
      desc.add("width", 360).setDescription("plot width, in points (1cm>~28pt)");
      desc.add("height", 200).setDescription("plot width, in points (1cm>~28pt)");
      desc.add("fontSize", 12);
      desc.add("lineWidth", 3);
      return desc;
    }

    const Drawer& draw(const utils::Graph1D& graph, const Mode& mode) const override {
      ::sciplot::Plot2D plot;
      build(plot, graph);
      style(plot, mode);
      plotAndSave({plot}, graph.name(), graph.title());
      return *this;
    }
    const Drawer& draw(const utils::Graph2D& graph, const Mode&) const override {
      ::sciplot::Plot3D plot;
      build(plot, graph);
      ::sciplot::Figure fig{{plot}};
      plotAndSave({plot}, graph.name(), graph.title());
      return *this;
    }
    const Drawer& draw(const utils::Hist1D& hist, const Mode& mode) const override {
      ::sciplot::Plot2D plot;
      build(plot, hist);
      style(plot, mode);
      plotAndSave({plot}, hist.name(), hist.title());
      return *this;
    }
    const Drawer& draw(const utils::Hist2D& hist, const Mode& mode) const override {
      ::sciplot::Plot3D plot;
      build(plot, hist);
      style(plot, mode);
      plotAndSave({plot}, hist.name(), hist.title());
      return *this;
    }

    const Drawer& draw(const utils::DrawableColl& coll,
                       const std::string& name,
                       const std::string& title,
                       const Mode& mode) const override {
      ::sciplot::Plot2D plot;
      int graph_line_style = 1, hist_line_style = 1;
      for (const auto& obj : coll) {
        if (obj->isHist1D())
          if (const auto* hist = dynamic_cast<const utils::Hist1D*>(obj); hist)
            build(plot, *hist, hist_line_style++);
        if (obj->isGraph1D())
          if (const auto* graph = dynamic_cast<const utils::Graph1D*>(obj); graph)
            build(plot, *graph, graph_line_style++);
      }
      style(plot, mode);
      if (!(mode & Mode::nostack))
        CG_WARNING("sciplot::Drawer::draw") << "Stacked plots are not yet available for this drawer. By default, all "
                                               "distributions will be drawn unstacked.";
      plotAndSave({plot}, name, title);
      return *this;
    }

  private:
    void plotAndSave(const std::vector<::sciplot::PlotVariant>& plots,
                     const std::string& name,
                     const std::string& title) const {
      ::sciplot::Figure fig{{plots}};
      fig.title("CepGen v" + version::tag + (!title.empty() ? " - "s + title : ""s));
      fig.palette(palette_name_);
      ::sciplot::Canvas canvas{{fig}};
      canvas.size(width_, height_);
      canvas.save(name + ".pdf");
    }
    void style(::sciplot::Plot& plot, const Mode& mode) const {
      plot.xtics().fontName(font_name_).fontSize(font_size_);
      plot.ytics().fontName(font_name_).fontSize(font_size_);
      plot.ztics().fontName(font_name_).fontSize(font_size_);
      if (mode & Mode::logx)
        plot.xtics().logscale(10);
      if (mode & Mode::logy)
        plot.ytics().logscale(10);
      if (mode & Mode::logz)
        plot.ztics().logscale(10);
      if (mode & Mode::grid) {
        auto& grid = plot.grid();
        grid.xtics().dashType(3);
        grid.ytics().dashType(3);
        grid.ztics().dashType(3);
      }
    }
    void build(::sciplot::Plot2D& plot, const utils::Graph1D& graph, int line_style = 1) const {
      plot.fontName(font_name_).fontSize(font_size_);
      plot.xlabel(graph.xAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.ylabel(graph.yAxis().label()).fontName(font_name_).fontSize(font_size_);
      if (!graph.empty()) {
        std::vector<double> x_values, y_values;
        for (const auto& [x_value, y_value] : graph.points()) {
          x_values.emplace_back(x_value.value);
          y_values.emplace_back(y_value);
        }
        if (const auto& x_range = graph.xAxis().range(); x_range.valid())
          plot.xrange(x_range.min(), x_range.max());
        if (const auto& y_range = graph.yAxis().range(); y_range.valid())
          plot.yrange(y_range.min(), y_range.max());
        plot.drawCurve(x_values, y_values).label(graph.title()).dashType(line_style).lineWidth(line_width_);
      }
    }
    void build(::sciplot::Plot3D& plot, const utils::Graph2D& graph) const {
      plot.fontName(font_name_).fontSize(font_size_);
      plot.xlabel(graph.xAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.ylabel(graph.yAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.zlabel(graph.zAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.border().clear();
      plot.border().bottomLeftFront();
      plot.border().bottomRightFront();
      plot.border().leftVertical();
      if (!graph.empty()) {
        std::vector<double> x_values, y_values, z_values;
        for (const auto& [x_value, yz_values] : graph.points())
          for (const auto& [y_value, z_value] : yz_values) {
            x_values.emplace_back(x_value.value);
            y_values.emplace_back(y_value.value);
            z_values.emplace_back(z_value);
          }
        if (const auto& x_range = graph.xAxis().range(); x_range.valid())
          plot.xrange(x_range.min(), x_range.max());
        if (const auto& y_range = graph.yAxis().range(); y_range.valid())
          plot.yrange(y_range.min(), y_range.max());
        if (const auto& z_range = graph.zAxis().range(); z_range.valid())
          plot.zrange(z_range.min(), z_range.max());
        plot.drawDots(x_values, y_values, z_values).label(graph.title()).lineWidth(line_width_);
      }
    }
    void build(::sciplot::Plot2D& plot, const utils::Hist1D& hist, int line_style = 1) const {
      plot.fontName(font_name_).fontSize(font_size_);
      plot.xlabel(hist.xAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.ylabel(hist.yAxis().label()).fontName(font_name_).fontSize(font_size_);
      if (!hist.empty()) {
        const auto bins = hist.bins(utils::Histogram::BinMode::both);
        std::vector<double> entries_per_bin, unc_per_bin;
        for (const auto& value : hist.values()) {
          entries_per_bin.emplace_back(value);
          unc_per_bin.emplace_back(value.uncertainty());
        }
        if (const auto& x_range = hist.xAxis().range(); x_range.valid())
          plot.xrange(x_range.min(), x_range.max());
        if (const auto& y_range = hist.yAxis().range(); y_range.valid())
          plot.yrange(y_range.min(), y_range.max());
        // explicitly set template arguments for now
        // (fix for issue highlighted and fixed in https://github.com/sciplot/sciplot/pull/118)
        plot.drawBoxesWithErrorBarsY<std::vector<double>, std::vector<double>, std::vector<double>>(
                bins, entries_per_bin, unc_per_bin)
            .label(hist.title())
            .dashType(line_style)
            .fillIntensity(0.33)
            .borderShow();
        plot.boxWidthRelative(1.);
      }
    }
    void build(::sciplot::Plot3D& plot, const utils::Hist2D& hist) const {
      plot.fontName(font_name_).fontSize(font_size_);
      plot.xlabel(hist.xAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.ylabel(hist.yAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.zlabel(hist.zAxis().label()).fontName(font_name_).fontSize(font_size_);
      plot.border().clear();
      plot.border().bottomLeftFront();
      plot.border().bottomRightFront();
      plot.border().leftVertical();
      if (!hist.empty()) {
        const auto bins_x = hist.binsX(utils::Histogram::BinMode::both),
                   bins_y = hist.binsY(utils::Histogram::BinMode::both);
        std::vector<double> all_bins_x, all_bins_y, entries_per_bin;
        for (size_t ix = 0; ix < hist.nbinsX(); ++ix) {
          const std::vector unitary_bins_x(bins_y.size(), bins_x.at(ix));
          all_bins_x.insert(all_bins_x.end(), unitary_bins_x.begin(), unitary_bins_x.end());
          all_bins_y.insert(all_bins_y.end(), bins_y.begin(), bins_y.end());
          for (size_t iy = 0; iy < hist.nbinsY(); ++iy)
            entries_per_bin.emplace_back(hist.value(ix, iy));
        }
        plot.drawDots(all_bins_x, all_bins_y, entries_per_bin)
            .label(hist.title())
            .lineWidth(line_width_)
            .fillIntensity(0.5)
            .borderShow();
        plot.boxWidthRelative(1.);
      }
    }

    const std::string palette_name_;
    const std::string font_name_;
    const int width_;
    const int height_;
    const int font_size_;
    const int line_width_;
  };
}  // namespace cepgen::sciplot
using SciplotDrawer = cepgen::sciplot::Drawer;
REGISTER_DRAWER("sciplot", SciplotDrawer);
