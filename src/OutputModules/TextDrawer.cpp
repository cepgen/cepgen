/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

namespace cepgen::utils {
  class TextDrawer : public Drawer {
  public:
    explicit TextDrawer(const ParametersList& params)
        : Drawer(params),
          CHAR('*'),
          ERR_CHAR('-'),
          NEG_CHAR('-'),
          MARKERS_CHAR("o.#@"),
          VALUES_CHAR(" .:oO0@%#"),
          width_(steerAs<int, size_t>("width")),
          colourise_(steer<bool>("colourise")) {}

    static ParametersDescription description() {
      auto desc = Drawer::description();
      desc.setDescription("Text-based drawing module");
      desc.add<int>("width", 50);
      desc.add<bool>("colourise", true).setDescription("colourise the output (for TTY-compatible displays)");
      return desc;
    }

    inline const TextDrawer& draw(const Graph1D& graph, const Mode& mode) const override {
      CG_LOG.log([this, &graph, &mode](auto& log) {
        if (!graph.name().empty())
          log << "plot of \"" << graph.name() << "\"\n";
        drawValues(log.stream(), graph, graph.points(), mode, colourise_);
      });
      return *this;
    }
    inline const TextDrawer& draw(const Graph2D& graph, const Mode& mode) const override {
      CG_LOG.log([this, &graph, &mode](auto& log) {
        if (!graph.name().empty())
          log << "plot of \"" << graph.name() << "\"\n";
        drawValues(log.stream(), graph, graph.points(), mode, colourise_);
      });
      return *this;
    }
    const TextDrawer& draw(const Hist1D&, const Mode&) const override;
    const TextDrawer& draw(const Hist2D&, const Mode&) const override;

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

    inline static std::string delatexify(const std::string& tok) { return utils::replaceAll(tok, {{"$", ""}}); }

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

  const TextDrawer& TextDrawer::draw(const Hist1D& hist, const Mode& mode) const {
    CG_LOG.log([this, &hist, &mode](auto& log) {
      if (!hist.name().empty())
        log << "plot of \"" << hist.name() << "\"\n";
      drawValues(log.stream(), hist, hist.axis(), mode, colourise_);
      const double bin_width = hist.range().range() / hist.nbins();
      log << "\t"
          << "bin width=" << utils::s("unit", bin_width, true) << ", "
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

  const TextDrawer& TextDrawer::draw(const Hist2D& hist, const Mode& mode) const {
    CG_LOG.log([this, &hist, &mode](auto& log) {
      if (!hist.name().empty())
        log << "plot of \"" << hist.name() << "\"\n";
      Drawable::dual_axis_t axes;
      for (size_t bin_x = 0; bin_x < hist.nbinsX(); ++bin_x) {
        const auto& range_x = hist.binRangeX(bin_x);
        auto& axis_x = axes[Drawable::coord_t{
            range_x.x(0.5), 0.5 * range_x.range(), utils::format("[%7.2f,%7.2f)", range_x.min(), range_x.max())}];
        for (size_t bin_y = 0; bin_y < hist.nbinsY(); ++bin_y) {
          const auto& range_y = hist.binRangeY(bin_y);
          axis_x[Drawable::coord_t{range_y.x(0.5), 0.5 * range_y.range(), utils::format("%+g", range_y.min())}] =
              hist.value(bin_x, bin_y);
        }
      }
      drawValues(log.stream(), hist, axes, mode, colourise_);
      const auto &x_range = hist.rangeX(), &y_range = hist.rangeY();
      const double bin_width_x = x_range.range() / hist.nbinsX(), bin_width_y = y_range.range() / hist.nbinsY();
      log << "\t"
          << " x-axis: "
          << "bin width=" << utils::s("unit", bin_width_x, true) << ", "
          << "mean=" << hist.meanX() << ","
          << "st.dev.=" << hist.rmsX() << "\n\t"
          << " y-axis: "
          << "bin width=" << utils::s("unit", bin_width_y, true) << ", "
          << "mean=" << hist.meanY() << ","
          << "st.dev.=" << hist.rmsY() << ",\n\t"
          << " integral=" << hist.integral();
      const auto& cnt = hist.outOfRange();
      if (cnt.total() > 0) {
        log << ", outside range (in/overflow):\n" << cnt;
      }
    });
    return *this;
  }

  const TextDrawer& TextDrawer::draw(const DrawableColl& objs,
                                     const std::string& name,
                                     const std::string&,
                                     const Mode& mode) const {
    CG_WARNING("TextDrawer:draw") << "Multi-plots is now only partially supported (no axes rescaling).";
    auto inside_plot = [](const std::string& str) -> std::string {
      std::istringstream ss(str);
      std::ostringstream out;
      for (std::string line; std::getline(ss, line);) {
        const auto tok = utils::split(line, ':');
        if (tok.size() == 3)
          out << tok.at(1) << "\n";
      }
      return out.str();
    };
    auto replace_plot = [](const std::string& orig, const std::string& new_plot) -> std::string {
      std::istringstream ss(orig), ssnew(new_plot);
      std::ostringstream out;
      for (std::string line; std::getline(ss, line);) {
        auto tok = utils::split(line, ':');
        if (tok.size() == 3) {
          std::getline(ssnew, tok[1]);
          tok[2].clear();
          out << utils::merge(tok, ":") << "\n";
        } else
          out << line << "\n";
      }
      return out.str();
    };
    std::stringstream buf, os_base;
    size_t num_plots = 0;
    auto add_plot = [this, &buf, &num_plots](const std::string& plt) {
      ++num_plots;
      if (plt.empty())
        return;
      std::istringstream ss(plt);
      std::ostringstream out;
      for (std::string line; std::getline(ss, line);) {
        std::string base(line.size(), ' ');
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
    for (const auto* obj : objs)
      if (obj->isHist1D()) {
        if (const auto* hist = dynamic_cast<const Hist1D*>(obj); hist) {
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
      } else if (obj->isGraph1D()) {
        if (const auto* gr = dynamic_cast<const Graph1D*>(obj); gr) {
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
      } else {
        CG_WARNING("TextDrawer:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
        continue;
      }
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
    const double max_val = std::max_element(axis.begin(), axis.end(), map_elements())->second *
                           ((mode & Mode::logy) ? 5. : 1.2),
                 min_val = std::min_element(axis.begin(), axis.end(), map_elements())->second;
    const double min_val_log = std::log(std::max(min_val, 1.e-10));
    const double max_val_log = std::log(std::min(max_val, 1.e+10));
    if (!dr.yAxis().label().empty()) {
      const auto y_label = delatexify(dr.yAxis().label());
      os << sep << std::string(std::max(0., 2. + width_ - y_label.size()), ' ') << y_label << "\n";
    }
    os << sep << utils::format("%-5.2f ", (mode & Mode::logy) ? std::exp(min_val_log) : min_val)
       << std::setw(width_ - 11) << std::left << ((mode & Mode::logy) ? "logarithmic scale" : "linear scale")
       << utils::format("%5.2e", (mode & Mode::logy) ? std::exp(max_val_log) : max_val) << "\n"
       << sep << std::string(width_ + 2, '.');  // abscissa axis
    size_t idx = 0;
    for (const auto& coord_set : axis) {
      const auto left_label =
          coord_set.first.label.empty() ? utils::format("%17g", coord_set.first.value) : coord_set.first.label;
      if (min_val == max_val) {
        os << "\n" << left_label << ":";
        if (idx == axis.size() / 2)
          os << std::string((width_ - kEmptyLabel.size()) / 2, ' ') << kEmptyLabel
             << std::string((width_ - kEmptyLabel.size()) / 2, ' ');
        else
          os << std::string(width_, ' ');
        os << ":";
      } else {
        const auto& val = coord_set.second;
        size_t i_value = 0ull, i_uncertainty = 0ull;
        {
          double val_dbl = width_, unc_dbl = width_;
          if (mode & Mode::logy) {
            val_dbl *= (val > 0. && max_val > 0.)
                           ? std::max((std::log(val) - min_val_log) / (max_val_log - min_val_log), 0.)
                           : 0.;
            unc_dbl *= (val > 0. && max_val > 0.)
                           ? std::max((std::log(val.uncertainty()) - min_val_log) / (max_val_log - min_val_log), 0.)
                           : 0.;
          } else if (max_val > 0.) {
            val_dbl *= (val - min_val) / (max_val - min_val);
            unc_dbl *= val.uncertainty() / (max_val - min_val);
          }
          i_value = std::ceil(val_dbl);
          i_uncertainty = std::ceil(unc_dbl);
        }
        os << "\n"
           << left_label << ":" << (i_value > i_uncertainty ? std::string(i_value - i_uncertainty, ' ') : "")
           << (i_uncertainty > 0 ? std::string(i_uncertainty, ERR_CHAR) : "")
           << (effects ? utils::boldify(std::string(1, CHAR)) : std::string(1, CHAR))
           << (i_uncertainty > 0 ? std::string(std::min(width_ - i_value - 1, i_uncertainty), ERR_CHAR) : "")
           << (i_value + i_uncertainty < width_ + 1 ? std::string(width_ - i_value - i_uncertainty - 1, ' ') : "")
           << ": " << utils::format("%6.2e +/- %6.2e", val, val.uncertainty());
      }
      ++idx;
    }
    os << "\n"
       << utils::format("%17s", delatexify(dr.xAxis().label()).c_str()) << ":" << std::string(width_, '.')
       << ":\n";  // 2nd abscissa axis
  }

  void TextDrawer::drawValues(
      std::ostream& os, const Drawable& dr, const Drawable::dual_axis_t& axes, const Mode& mode, bool effects) const {
    const std::string sep(17, ' ');
    if (!dr.yAxis().label().empty()) {
      const auto y_label = delatexify(dr.yAxis().label());
      os << sep << std::string(std::max(0., 2. + width_ - y_label.size()), ' ') << y_label << "\n";
    }
    // find the maximum element of the graph
    double min_val = -Limits::INVALID, max_val = Limits::INVALID;
    double min_log_value = -3.;
    for (const auto& x_value : axes) {
      min_val = std::min(
          min_val,
          static_cast<double>(std::min_element(x_value.second.begin(), x_value.second.end(), map_elements())->second));
      max_val = std::max(
          max_val,
          static_cast<double>(std::max_element(x_value.second.begin(), x_value.second.end(), map_elements())->second));
      if (mode & Mode::logz)
        for (const auto& y_value : x_value.second)
          if (y_value.second > 0.)
            min_log_value = std::min(min_log_value, std::log(y_value.second / max_val));
    }
    const auto& y_axis = axes.begin()->second;
    os << sep << utils::format("%-5.2f", y_axis.begin()->first.value) << std::string(axes.size() - 11, ' ')
       << utils::format("%5.2e", y_axis.rbegin()->first.value) << "\n"
       << utils::format("%17s", delatexify(dr.xAxis().label()).c_str())
       << std::string(1 + y_axis.size() + 1, '.');  // abscissa axis
    size_t idx = 0;
    for (const auto& x_value : axes) {
      os << "\n"
         << (x_value.first.label.empty() ? utils::format("%16g ", x_value.first.value) : x_value.first.label) << ":";
      if (min_val == max_val) {
        if (idx == axes.size() / 2)
          os << std::string((width_ - kEmptyLabel.size()) / 2, ' ') << kEmptyLabel
             << std::string((width_ - kEmptyLabel.size()) / 2, ' ');
        else
          os << std::string(width_, ' ');
      } else {
        for (const auto& y_value : x_value.second) {
          const double val = y_value.second;
          double val_norm = 0.;
          if (mode & Mode::logz)
            val_norm =
                positive(val) ? std::max(0., (std::log(val / max_val) - min_log_value) / fabs(min_log_value)) : 0.;
          else
            val_norm = val / max_val;
          if (std::isnan(val_norm)) {
            os << (effects ? utils::colourise("!", kColours.at(0)) : "!");
            continue;
          }
          const short sign = (val_norm == 0. ? 0 : val_norm / fabs(val_norm));
          val_norm *= sign;
          if (sign == -1)
            os << (effects ? utils::colourise(std::string(1, NEG_CHAR), kColours.at(0)) : std::string(1, NEG_CHAR));
          else {
            size_t ch_id = ceil(val_norm * (VALUES_CHAR.size() - 1));
            size_t col_id = 1 + val_norm * (kColours.size() - 2);
            os << (effects ? utils::colourise(std::string(1, VALUES_CHAR.at(ch_id)),
                                              kColours.at(col_id),
                                              (val_norm > 0.75 ? utils::Modifier::bold : utils::Modifier::reset))
                           : std::string(1, VALUES_CHAR.at(ch_id)));
          }
        }
      }
      os << ":";
      ++idx;
    }
    std::vector<std::string> y_label;
    std::transform(y_axis.begin(), y_axis.end(), std::back_inserter(y_label), [](auto& bin) {
      return bin.first.label.empty() ? utils::format("%+g", bin.first.value) : bin.first.label;
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
    os << "\n"
       << sep << ":" << std::string(y_axis.size(), '.') << ": "  // 2nd abscissa axis
       << delatexify(dr.yAxis().label()) << "\n\t"
       << "(scale: \"" << VALUES_CHAR << "\", ";
    for (size_t i = 0; i < kColours.size(); ++i)
      os << (effects ? utils::colourise("*", kColours.at(i)) : "") << (i == 0 ? "|" : "");
    os << ")\n";
  }

}  // namespace cepgen::utils
using cepgen::utils::TextDrawer;
REGISTER_DRAWER("text", TextDrawer);
