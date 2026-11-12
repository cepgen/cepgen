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

#include <array>
#include <cmath>
#include <iomanip>

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::utils;
using namespace std::string_literals;

class TextDrawer final : public Drawer {
public:
  explicit TextDrawer(const ParametersList& params)
      : Drawer(params),
        CHAR('*'),
        ERR_CHAR('-'),
        NEG_CHAR('-'),
        MARKERS_CHAR("o+*#@."),
        VALUES_CHAR(" .:oO0@%#"),
        width_(static_cast<size_t>(std::max(5, steer<int>("width")))),
        colourise_(steer<bool>("colourise")) {}

  static ParametersDescription description() {
    auto desc = Drawer::description();
    desc.setDescription("Text-based drawing module");
    desc.add("width", 50);
    desc.add("colourise", true).setDescription("colourise the output (for TTY-compatible displays)");
    return desc;
  }

  const TextDrawer& draw(const Graph1D& graph, const Mode& mode) const override {
    CG_LOG.log([this, &graph, &mode](auto& log) {
      if (!graph.name().empty())
        log << "plot of \"" << graph.name() << "\"\n";
      drawValues(log.stream(), graph, graph.points(), mode, colourise_);
    });
    return *this;
  }
  const TextDrawer& draw(const Graph2D& graph, const Mode& mode) const override {
    CG_LOG.log([this, &graph, &mode](auto& log) {
      if (!graph.name().empty())
        log << "plot of \"" << graph.name() << "\"\n";
      drawValues(log.stream(), graph, graph.points(), mode, colourise_);
    });
    return *this;
  }
  const TextDrawer& draw(const Hist1D& hist, const Mode& mode) const override {
    CG_LOG.log([this, &hist, &mode](auto& log) {
      if (!hist.name().empty())
        log << "plot of \"" << hist.name() << "\"\n";
      drawValues(log.stream(), hist, hist.axis(), mode, colourise_);
      const double bin_width = hist.range().range() / hist.nbins();
      log << "\t"
          << "bin width=" << s("unit", bin_width, true) << ", "
          << "mean=" << hist.mean() << ", "
          << "std.dev.=" << hist.rms() << "\n\t"
          << "integral.=" << hist.integral();
      if (hist.underflow() > 0ull)
        log << ", underflow: " << hist.underflow();
      if (hist.overflow() > 0ull)
        log << ", overflow: " << hist.overflow();
    });
    return *this;
  }
  const TextDrawer& draw(const Hist2D& hist, const Mode& mode) const override {
    CG_LOG.log([this, &hist, &mode](auto& log) {
      if (!hist.name().empty())
        log << "plot of \"" << hist.name() << "\"\n";
      Drawable::dual_axis_t axes;
      for (size_t bin_x = 0; bin_x < hist.nbinsX(); ++bin_x) {
        const auto& range_x = hist.binRangeX(bin_x);
        auto& axis_x = axes[Drawable::coord_t{
            range_x.x(0.5), 0.5 * range_x.range(), format("[%7.2f,%7.2f)", range_x.min(), range_x.max())}];
        for (size_t bin_y = 0; bin_y < hist.nbinsY(); ++bin_y) {
          const auto& range_y = hist.binRangeY(bin_y);
          axis_x[Drawable::coord_t{range_y.x(0.5), 0.5 * range_y.range(), format("%+g", range_y.min())}] =
              hist.value(bin_x, bin_y);
        }
      }
      drawValues(log.stream(), hist, axes, mode, colourise_);
      log << "\t"
          << " x-axis: "
          << "bin width=" << s("unit", hist.rangeX().range() / hist.nbinsX(), true) << ", "
          << "mean=" << hist.meanX() << ","
          << "st.dev.=" << hist.rmsX() << "\n\t"
          << " y-axis: "
          << "bin width=" << s("unit", hist.rangeY().range() / hist.nbinsY(), true) << ", "
          << "mean=" << hist.meanY() << ","
          << "st.dev.=" << hist.rmsY() << ",\n\t"
          << " integral=" << hist.integral();
      if (const auto& content = hist.outOfRange(); content.total() > 0)
        log << ", outside range (in/overflow):\n" << content;
    });
    return *this;
  }

  const TextDrawer& draw(const DrawableColl&,
                         const std::string& name = "",
                         const std::string& title = "",
                         const Mode& mode = Mode::none) const override;

private:
  friend class Drawable;
  static const std::array<Colour, 7> kColours;
  static const std::string kEmptyLabel;

  void drawValues(std::ostream&, const Drawable&, const Drawable::axis_t&, const Mode&, bool effects) const;
  void drawValues(std::ostream&, const Drawable&, const Drawable::dual_axis_t&, const Mode&, bool effects) const;

  const char CHAR, ERR_CHAR, NEG_CHAR;
  const std::string MARKERS_CHAR, VALUES_CHAR;

  static std::string delatexify(const std::string& tok) { return replaceAll(tok, {{"$", ""}}); }
  /// Build a centred "empty plot" banner of exactly the requested width
  static std::string emptyBanner(size_t width) {
    if (width <= kEmptyLabel.size())
      return kEmptyLabel.substr(0, width);
    const auto left = (width - kEmptyLabel.size()) / 2, right = width - kEmptyLabel.size() - left;
    return std::string(left, ' ') + kEmptyLabel + std::string(right, ' ');
  }

  /// Sorting helper for the axis metadata container
  struct map_elements {
    /// Sort two values in an axis
    bool operator()(const std::pair<Drawable::coord_t, Value>& lhs,
                    const std::pair<Drawable::coord_t, Value>& rhs) const {
      return lhs.second < rhs.second;
    }
  };

  const size_t width_{0};
  const bool colourise_{true};
};

const std::array<Colour, 7> TextDrawer::kColours = {
    Colour::red, Colour::cyan, Colour::blue, Colour::magenta, Colour::green, Colour::yellow, Colour::reset};

const std::string TextDrawer::kEmptyLabel = "E M P T Y ";

const TextDrawer& TextDrawer::draw(const DrawableColl& objects,
                                   const std::string& name,
                                   const std::string&,
                                   const Mode& mode) const {
  CG_WARNING("TextDrawer:draw") << "Multi-plots is now only partially supported (no axes rescaling).";
  auto inside_plot = [](const std::string& str) -> std::string {
    std::istringstream ss(str);
    std::ostringstream out;
    for (std::string line; std::getline(ss, line);)
      if (const auto tok = split(line, ':'); tok.size() == 3)
        out << tok.at(1) << "\n";
    return out.str();
  };
  auto replace_plot = [](const std::string& orig, const std::string& new_plot) -> std::string {
    std::istringstream ss(orig), ss_new(new_plot);
    std::ostringstream out;
    for (std::string line; std::getline(ss, line);) {
      if (auto token = split(line, ':'); token.size() == 3) {
        std::getline(ss_new, token[1]);
        token[2].clear();
        out << merge(token, ":") << "\n";
      } else
        out << line << "\n";
    }
    return out.str();
  };
  std::stringstream buf, os_base;
  size_t num_plots = 0;
  auto add_plot = [this, &buf, &num_plots](const std::string& plot) {
    ++num_plots;
    if (plot.empty())
      return;
    std::istringstream ss(plot);
    std::ostringstream out;
    for (std::string line; std::getline(ss, line);) {
      auto base = line.empty() ? ""s : std::string(line.size(), ' ');
      if (!buf.str().empty() && !std::getline(buf, base)) {
        CG_WARNING("TextDrawer:draw") << "Invalid plot to be produced... Aborting the multiplot.";
        return;
      }
      for (size_t j = 0; j < line.size(); ++j) {
        if (line[j] == CHAR)
          base[j] = (num_plots > 1 ? MARKERS_CHAR[num_plots - 2] : CHAR);
        else if (line[j] == ERR_CHAR)
          base[j] = ERR_CHAR;
      }
      out << base << "\n";
    }
    buf.str(out.str());
  };
  std::vector<std::string> plt_names;
  for (const auto* object : objects)
    if (object->isHist1D()) {
      if (const auto* hist = dynamic_cast<const Hist1D*>(object); hist) {
        if (os_base.str().empty()) {
          drawValues(os_base, *hist, hist->axis(), mode, false);
          add_plot(inside_plot(os_base.str()));
        } else {
          std::ostringstream os;
          drawValues(os, *hist, hist->axis(), mode, false);
          add_plot(inside_plot(os.str()));
        }
        plt_names.emplace_back(hist->name());
      }
    } else if (object->isGraph1D()) {
      if (const auto* gr = dynamic_cast<const Graph1D*>(object); gr) {
        if (os_base.str().empty()) {
          drawValues(os_base, *gr, gr->points(), mode, false);
          add_plot(inside_plot(os_base.str()));
        } else {
          std::ostringstream os;
          drawValues(os, *gr, gr->points(), mode, false);
          add_plot(inside_plot(os.str()));
        }
        plt_names.emplace_back(gr->name());
      }
    } else
      CG_WARNING("TextDrawer:draw") << "Cannot add drawable '" << object->name() << "' to the stack.";

  CG_LOG.log([this, &name, &replace_plot, &os_base, &buf, &num_plots, &plt_names](auto& log) {
    if (!name.empty())
      log << "plot of \"" << name << "\"\n";
    log << replace_plot(os_base.str(), buf.str());
    if (num_plots > 1)
      log << "\tLegend:\n\t  " << CHAR << ": " << plt_names.at(0);
    for (size_t i = 1; i < num_plots; ++i)
      log << "\n\t  " << MARKERS_CHAR[i - 1] << ": " << plt_names.at(i);
  });
  return *this;
}

void TextDrawer::drawValues(
    std::ostream& os, const Drawable& dr, const Drawable::axis_t& axis, const Mode& mode, bool effects) const {
  const std::string sep(17, ' ');
  const auto values_range = Limits{std::min_element(axis.begin(), axis.end(), map_elements())->second,
                                   std::max_element(axis.begin(), axis.end(), map_elements())->second *
                                       (mode & Mode::logy ? 5. : 1.2)},
             log_values_range = values_range.compute(std::log).truncate({-10., 10.});
  if (!dr.yAxis().label().empty()) {
    const auto y_label = delatexify(dr.yAxis().label());
    os << sep;
    if (const auto length = 2. + width_ - y_label.size(); length > 0)
      os << std::string(length, ' ');
    os << y_label << "\n";
  }
  os << sep << format("%-5.2f ", mode & Mode::logy ? std::exp(log_values_range.min()) : values_range.min())
     << std::setw(width_ > 11 ? width_ - 11 : 0) << std::left
     << (mode & Mode::logy ? "logarithmic scale" : "linear scale")
     << format("%5.2e", mode & Mode::logy ? std::exp(log_values_range.max()) : values_range.max()) << "\n"
     << sep << std::string(width_ + 2, '.');  // abscissa axis
  size_t idx = 0;
  const auto empty_plot = dr.empty();
  for (const auto& [coordinates, value] : axis) {
    const auto left_label = coordinates.label.empty() ? format("%17g", coordinates.value) : coordinates.label;
    if (empty_plot) {  // empty plot; display the "EMPTY" banner
      os << "\n"
         << left_label << ":" << (idx == axis.size() / 2 ? emptyBanner(width_) : std::string(width_, ' ')) << ":";
    } else {  // value is valid
      const double val = value, unc = value.uncertainty();
      // column (in [0, width-1]) at which a given ordinate is drawn
      const auto column = [this, &mode, &values_range, &log_values_range](double y) -> size_t {
        double frac = 0.;
        if (mode & Mode::logy) {  // log scale
          if (y > 0. && values_range.max() > 0.)
            frac = std::log(y) / log_values_range;
        } else if (values_range.max() > 0.)  // linear scale
          frac = y / values_range;
        if (!std::isfinite(frac) || frac <= 0.)
          return 0;
        return std::min(width_ - 1, static_cast<size_t>(std::ceil(width_ * std::min(frac, 1.))));
      };
      const auto i_value = column(val);
      size_t i_unc_low = 0, i_unc_high = 0;  // error bar extension on each side of the marker
      if (unc > 0.) {
        if (mode & Mode::logy) {  // asymmetric in log scale: span [log(y-dy), log(y+dy)]
          i_unc_low = i_value - std::min(i_value, column(val - unc));
          i_unc_high = std::max(i_value, column(val + unc)) - i_value;
        } else if (const auto frac = unc / values_range.range(); std::isfinite(frac) && frac > 0.)
          i_unc_low = i_unc_high = static_cast<size_t>(std::ceil(width_ * frac));
      }
      i_unc_low = std::min(i_unc_low, i_value);
      i_unc_high = std::min(i_unc_high, width_ - 1 - i_value);
      os << "\n"
         << left_label << ":" << std::string(i_value - i_unc_low, ' ') << std::string(i_unc_low, ERR_CHAR)
         << (effects ? boldify(std::string{CHAR}) : std::string{CHAR}) << std::string(i_unc_high, ERR_CHAR)
         << std::string(width_ - 1 - i_value - i_unc_high, ' ') << ": "
         << (mode & Mode::logy ? format("% 6.2e +/- %6.2e", value, value.uncertainty())
                               : format("% 6f +/- %6f", value, value.uncertainty()));
    }
    ++idx;
  }
  os << "\n"
     << format("%17s", delatexify(dr.xAxis().label()).c_str()) << ":" << std::string(width_, '.')
     << ":\n";  // 2nd abscissa axis
}

void TextDrawer::drawValues(
    std::ostream& os, const Drawable& dr, const Drawable::dual_axis_t& axes, const Mode& mode, bool effects) const {
  const std::string sep(17, ' ');
  if (!dr.yAxis().label().empty()) {
    const auto y_label = delatexify(dr.yAxis().label());
    os << sep;
    if (2 + width_ > y_label.size())
      os << std::string(2 + width_ - y_label.size(), ' ');
    os << y_label << "\n";
  }
  // find the maximum element of the graph
  double min_val = -Limits::INVALID, max_val = Limits::INVALID;
  double min_log_value = -3.;
  for (const auto& [x_coordinate, y_values] : axes) {
    min_val = std::min(min_val,
                       static_cast<double>(std::min_element(y_values.begin(), y_values.end(), map_elements())->second));
    max_val = std::max(max_val,
                       static_cast<double>(std::max_element(y_values.begin(), y_values.end(), map_elements())->second));
    if (mode & Mode::logz)
      for (const auto& [y_coordinate, y_value] : y_values)
        if (y_value > 0.)
          min_log_value = std::min(min_log_value, std::log(y_value / max_val));
  }
  const auto& y_axis = axes.begin()->second;
  os << sep << format("%-5.2f", y_axis.begin()->first.value);
  if (axes.size() > 11)
    os << std::string(axes.size() - 11, ' ');
  os << format("%5.2e", y_axis.rbegin()->first.value) << "\n"
     << format("%17s", delatexify(dr.xAxis().label()).c_str())
     << std::string(1 + y_axis.size() + 1, '.');  // abscissa axis
  size_t idx = 0;
  const auto empty_plot = dr.empty();
  for (const auto& [x_coordinate, y_values] : axes) {
    os << "\n" << (x_coordinate.label.empty() ? format("%16g ", x_coordinate.value) : x_coordinate.label) << ":";
    if (empty_plot) {  // empty plot; display the "EMPTY" banner
      os << (idx == axes.size() / 2 ? emptyBanner(width_) : std::string(width_, ' '));
    } else {
      for (const auto& [y_coordinate, y_value] : y_values) {
        double val_norm = 0.;
        if (mode & Mode::logz)
          val_norm = positive(static_cast<double>(y_value))
                         ? std::max(0., (std::log(y_value / max_val) - min_log_value) / fabs(min_log_value))
                         : 0.;
        else
          val_norm = y_value / max_val;
        if (std::isnan(val_norm)) {
          os << (effects ? colourise("!", kColours.at(0)) : "!");
          continue;
        }
        const short sign = (val_norm == 0. ? 0 : val_norm / fabs(val_norm));
        val_norm *= sign;
        if (sign == -1)
          os << (effects ? colourise(std::string(1, NEG_CHAR), kColours.at(0)) : std::string(1, NEG_CHAR));
        else {
          size_t ch_id = std::ceil(val_norm * (VALUES_CHAR.size() - 1));
          size_t col_id = 1 + val_norm * (kColours.size() - 2);
          os << (effects ? colourise(std::string(1, VALUES_CHAR.at(ch_id)),
                                     kColours.at(col_id),
                                     val_norm > 0.75 ? Modifier::bold : Modifier::reset)
                         : std::string(1, VALUES_CHAR.at(ch_id)));
        }
      }
    }
    os << ":";
    ++idx;
  }
  std::vector<std::string> y_label;
  std::transform(y_axis.begin(), y_axis.end(), std::back_inserter(y_label), [](auto& bin) {
    return bin.first.label.empty() ? format("%+g", bin.first.value) : bin.first.label;
  });
  struct string_length {
    bool operator()(const std::string& a, const std::string& b) const { return a.size() < b.size(); }
  };
  for (size_t i = 0; i < std::max_element(y_label.begin(), y_label.end(), string_length())->size(); ++i) {
    os << "\n" << sep << ":";
    for (const auto& lab : y_label)
      os << (lab.size() > i ? lab.at(i) : ' ');
    os << ":";
  }
  os << "\n" << sep << ":";
  if (!y_axis.empty())
    os << std::string(y_axis.size(), '.');
  os << ": "  // 2nd abscissa axis
     << delatexify(dr.yAxis().label()) << "\n\t"
     << "(scale: \"" << VALUES_CHAR << "\", ";
  for (size_t i = 0; i < kColours.size(); ++i)
    os << (effects ? colourise("*", kColours.at(i)) : ""s) << (i == 0 ? "|" : ""s);
  os << ")\n";
}
REGISTER_DRAWER("text", TextDrawer);
