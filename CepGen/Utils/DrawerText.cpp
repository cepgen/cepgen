/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include <cstring>
#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawable.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class DrawerText : public Drawer {
    public:
      explicit DrawerText(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("Text-based drawing module");
        desc.add<int>("width", 50);
        return desc;
      }

      const DrawerText& draw(const Graph1D&, const Mode&) const override;
      const DrawerText& draw(const Graph2D&, const Mode&) const override;
      const DrawerText& draw(const Hist1D&, const Mode&) const override;
      const DrawerText& draw(const Hist2D&, const Mode&) const override;

      const DrawerText& draw(const DrawableColl&,
                             const std::string& name = "",
                             const Mode& mode = Mode::none) const override;

    private:
      friend class Drawable;

      void drawValues(std::ostream&, const Drawable&, const Drawable::axis_t&, const Mode&, bool effects = true) const;
      void drawValues(
          std::ostream&, const Drawable&, const Drawable::dualaxis_t&, const Mode&, bool effects = true) const;

      static constexpr char CHAR = '*', ERR_CHAR = '-';
      static constexpr const char* CHAR_ALT = "o.#@";
      // greyscale ascii art from http://paulbourke.net/dataformats/asciiart/
      //static constexpr const char* CHARS = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
      //static constexpr const char* CHARS = " .:-=+*#%@";
      static constexpr const char* CHARS = " .:oO0@%#";
      static const int kColours[];
      static constexpr const char NEG_CHAR = '-';

      /// Sorting helper for the axis metadata container
      struct map_elements {
        /// Sort two values in an axis
        bool operator()(const std::pair<Drawable::coord_t, Drawable::value_t>& lhs,
                        const std::pair<Drawable::coord_t, Drawable::value_t>& rhs) {
          return lhs.second.value < rhs.second.value;
        }
      };

      const size_t width_{0};
    };

    DrawerText::DrawerText(const ParametersList& params) : Drawer(params), width_(steerAs<int, size_t>("width")) {}

    const DrawerText& DrawerText::draw(const Graph1D& graph, const Mode& mode) const {
      CG_LOG.log([&](auto& log) {
        if (!graph.name().empty())
          log << "plot of \"" << graph.name() << "\"\n";
        drawValues(log.stream(), graph, graph.points(), mode);
      });
      return *this;
    }

    const DrawerText& DrawerText::draw(const Graph2D& graph, const Mode& mode) const {
      CG_LOG.log([&](auto& log) {
        if (!graph.name().empty())
          log << "plot of \"" << graph.name() << "\"\n";
        drawValues(log.stream(), graph, graph.points(), mode);
      });
      return *this;
    }

    const DrawerText& DrawerText::draw(const Hist1D& hist, const Mode& mode) const {
      CG_LOG.log([&](auto& log) {
        if (!hist.name().empty())
          log << "plot of \"" << hist.name() << "\"\n";
        drawValues(log.stream(), hist, hist.axis(), mode);
        const double bin_width = hist.range().range() / hist.nbins();
        log << "\tbin width=" << utils::s("unit", bin_width, true) << ", "
            << "mean=" << hist.mean() << ", "
            << "st.dev.=" << hist.rms() << "\n\t"
            << "integr.=" << hist.integral();
        if (hist.underflow() > 0ull)
          log << ", underflow: " << hist.underflow();
        if (hist.overflow() > 0ull)
          log << ", overflow: " << hist.overflow();
      });
      return *this;
    }

    const DrawerText& DrawerText::draw(const Hist2D& hist, const Mode& mode) const {
      CG_LOG.log([&](auto& log) {
        if (!hist.name().empty())
          log << "plot of \"" << hist.name() << "\"\n";
        Drawable::dualaxis_t axes;
        for (size_t binx = 0; binx < hist.nbinsX(); ++binx) {
          const auto& range_x = hist.binRangeX(binx);
          auto& axis_x =
              axes[Drawable::coord_t{range_x.x(0.5), utils::format("[%7.2f,%7.2f)", range_x.min(), range_x.max())}];
          for (size_t biny = 0; biny < hist.nbinsY(); ++biny) {
            const auto& range_y = hist.binRangeY(biny);
            axis_x[Drawable::coord_t{range_y.x(0.5), utils::format("%+g", range_y.min())}] =
                Drawable::value_t{hist.value(binx, biny), hist.valueUnc(binx, biny)};
          }
        }
        drawValues(log.stream(), hist, axes, mode);
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
        const auto& cnt = hist.content();
        if (cnt.total() > 0) {
          log << ", outside range (in/overflow):\n"
              << utils::format(
                     "%10zu | %10zu | %10zu\n"
                     "%10zu | %10s | %10zu\n"
                     "%10zu | %10zu | %10zu",
                     cnt.LT_LT,
                     cnt.LT_IN,
                     cnt.LT_GT,
                     cnt.IN_LT,
                     "-",
                     cnt.IN_GT,
                     cnt.GT_LT,
                     cnt.GT_IN,
                     cnt.GT_GT);
        }
      });
      return *this;
    }

    const DrawerText& DrawerText::draw(const DrawableColl& objs, const std::string& name, const Mode& mode) const {
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
      size_t num_plts = 0;
      auto add_plot = [&buf, &num_plts](const std::string& plt) {
        ++num_plts;
        if (plt.empty())
          return;
        std::istringstream ss(plt);
        std::ostringstream out;
        for (std::string line; std::getline(ss, line);) {
          std::string base(line.size(), ' ');
          if (!buf.str().empty() && !std::getline(buf, base)) {
            CG_WARNING("DrawerText:draw") << "Invalid plot to be produced... Aborting the multiplot.";
            return;
          }
          for (size_t j = 0; j < line.size(); ++j) {
            if (line[j] == CHAR)
              base[j] = (num_plts > 1 ? CHAR_ALT[num_plts - 2] : CHAR);
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
          const auto* hist = dynamic_cast<const Hist1D*>(obj);
          if (os_base.str().empty()) {
            drawValues(os_base, *hist, hist->axis(), mode, false);
            add_plot(inside_plot(os_base.str()));
          } else {
            std::ostringstream os;
            drawValues(os, *hist, hist->axis(), mode, false);
            add_plot(inside_plot(os.str()));
          }
          plt_names.emplace_back(hist->name());
        } else if (obj->isGraph1D()) {
          const auto* gr = dynamic_cast<const Graph1D*>(obj);
          if (os_base.str().empty()) {
            drawValues(os_base, *gr, gr->points(), mode, false);
            add_plot(inside_plot(os_base.str()));
          } else {
            std::ostringstream os;
            drawValues(os, *gr, gr->points(), mode, false);
            add_plot(inside_plot(os.str()));
          }
          plt_names.emplace_back(gr->name());
        } else {
          CG_WARNING("DrawerText:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
          continue;
        }
      CG_LOG.log([&](auto& log) {
        if (!name.empty())
          log << "plot of \"" << name << "\"\n";
        log << replace_plot(os_base.str(), buf.str());
        if (num_plts > 1)
          log << "\tLegend:\n\t  " << CHAR << ": " << plt_names.at(0);
        for (size_t i = 1; i < num_plts; ++i)
          log << "\n\t  " << CHAR_ALT[i - 1] << ": " << plt_names.at(i);
      });
      return *this;
    }

    void DrawerText::drawValues(
        std::ostream& os, const Drawable& dr, const Drawable::axis_t& axis, const Mode& mode, bool effect) const {
      const std::string sep(17, ' ');
      const double max_val = std::max_element(axis.begin(), axis.end(), map_elements())->second.value *
                             (mode & Mode::logy ? 5. : 1.2),
                   min_val = std::min_element(axis.begin(), axis.end(), map_elements())->second.value;
      const double min_val_log = std::log(std::max(min_val, 1.e-10));
      const double max_val_log = std::log(std::min(max_val, 1.e+10));
      if (!dr.yLabel().empty())
        os << sep << std::string(std::max(0., 2. + width_ - dr.yLabel().size()), ' ') << dr.yLabel() << "\n";
      os << sep << utils::format("%-5.2f ", mode & Mode::logy ? std::exp(min_val_log) : min_val)
         << std::setw(width_ - 11) << std::left << (mode & Mode::logy ? "logarithmic scale" : "linear scale")
         << utils::format("%5.2e", mode & Mode::logy ? std::exp(max_val_log) : max_val) << "\n"
         << sep << std::string(width_ + 2, '.');  // abscissa axis
      size_t idx = 0;
      for (const auto& coord_set : axis) {
        const auto left_label =
            coord_set.first.label.empty() ? utils::format("%17g", coord_set.first.value) : coord_set.first.label;
        if (min_val == max_val) {
          os << "\n" << left_label << ":";
          if (idx == axis.size() / 2) {
            const std::string empty = "E M P T Y ";
            os << std::string((width_ - empty.size()) / 2, ' ') << empty
               << std::string((width_ - empty.size()) / 2, ' ');
          } else
            os << std::string(width_, ' ');
          os << ":";
        } else {
          const auto& set = coord_set.second;
          const double val = set.value, unc = set.value_unc;
          size_t ival = 0ull, ierr = 0ull;
          {
            double val_dbl = width_, unc_dbl = width_;
            if (mode & Mode::logy) {
              val_dbl *= (val > 0. && max_val > 0.)
                             ? std::max((std::log(val) - min_val_log) / (max_val_log - min_val_log), 0.)
                             : 0.;
              unc_dbl *= (val > 0. && max_val > 0.)
                             ? std::max((std::log(unc) - min_val_log) / (max_val_log - min_val_log), 0.)
                             : 0.;
            } else if (max_val > 0.) {
              val_dbl *= (val - min_val) / (max_val - min_val);
              unc_dbl *= unc / (max_val - min_val);
            }
            ival = std::ceil(val_dbl);
            ierr = std::ceil(unc_dbl);
          }
          os << "\n"
             << left_label << ":" << (ival > ierr ? std::string(ival - ierr, ' ') : "")
             << (ierr > 0 ? std::string(ierr, ERR_CHAR) : "")
             << (effect ? utils::boldify(std::string(1, CHAR)) : std::string(1, CHAR))
             << (ierr > 0 ? std::string(std::min(width_ - ival - 1, ierr), ERR_CHAR) : "")
             << (ival + ierr < width_ + 1 ? std::string(width_ - ival - ierr - 1, ' ') : "") << ": "
             << utils::format("%6.2e +/- %6.2e", val, unc);
        }
        ++idx;
      }
      os << "\n"
         << utils::format("%17s", dr.xLabel().c_str()) << ":" << std::string(width_, '.')
         << ":\n";  // 2nd abscissa axis
    }

    const int DrawerText::kColours[] = {(int)Colour::red,
                                        (int)Colour::cyan,
                                        (int)Colour::blue,
                                        (int)Colour::magenta,
                                        (int)Colour::green,
                                        (int)Colour::yellow,
                                        (int)Colour::reset};

    void DrawerText::drawValues(
        std::ostream& os, const Drawable& dr, const Drawable::dualaxis_t& axes, const Mode& mode, bool effects) const {
      const std::string sep(17, ' ');
      if (!dr.yLabel().empty())
        os << sep << std::string(std::max(0., 2. + width_ - dr.yLabel().size()), ' ') << dr.yLabel() << "\n";
      // find the maximum element of the graph
      double min_val = -Limits::INVALID, max_val = Limits::INVALID;
      double min_logval = -3.;
      for (const auto& xval : axes) {
        min_val =
            std::min(min_val, std::min_element(xval.second.begin(), xval.second.end(), map_elements())->second.value);
        max_val =
            std::max(max_val, std::max_element(xval.second.begin(), xval.second.end(), map_elements())->second.value);
        if (mode & Mode::logz)
          for (const auto& yval : xval.second)
            if (yval.second.value > 0.)
              min_logval = std::min(min_logval, std::log(yval.second.value / max_val));
      }
      const auto& y_axis = axes.begin()->second;
      os << sep << utils::format("%-5.2f", y_axis.begin()->first.value) << std::string(axes.size() - 11, ' ')
         << utils::format("%5.2e", y_axis.rbegin()->first.value) << "\n"
         << utils::format("%17s", dr.xLabel().c_str()) << std::string(1 + y_axis.size() + 1, '.');  // abscissa axis
      size_t idx = 0;
      for (const auto& xval : axes) {
        os << "\n" << (xval.first.label.empty() ? utils::format("%16g ", xval.first.value) : xval.first.label) << ":";
        if (min_val == max_val) {
          if (idx == axes.size() / 2)
            os << std::string((width_ - 10) / 2, ' ') << "E M P T Y " << std::string((width_ - 10) / 2, ' ');
          else
            os << std::string(width_, ' ');
        } else {
          for (const auto& yval : xval.second) {
            const double val = yval.second.value;
            double val_norm = 0.;
            if (mode & Mode::logz)
              val_norm = val <= 0. ? 0. : std::max(0., (std::log(val / max_val) - min_logval) / fabs(min_logval));
            else
              val_norm = val / max_val;
            if (std::isnan(val_norm)) {
              os << (effects ? utils::colourise("!", (utils::Colour)kColours[0]) : "!");
              continue;
            }
            const short sign = (val_norm == 0. ? 0 : val_norm / fabs(val_norm));
            val_norm *= sign;
            if (sign == -1)
              os << (effects ? utils::colourise(std::string(1, NEG_CHAR), (utils::Colour)kColours[0])
                             : std::string(1, NEG_CHAR));
            else {
              size_t ch_id = ceil(val_norm * (strlen(CHARS) - 1));
              size_t col_id = 1 + (val_norm * (sizeof(kColours) / (sizeof(int)) - 2));
              os << (effects ? utils::colourise(std::string(1, CHARS[ch_id]),
                                                (utils::Colour)kColours[col_id],
                                                (val_norm > 0.75 ? utils::Modifier::bold : utils::Modifier::reset))
                             : std::string(1, CHARS[ch_id]));
            }
          }
        }
        os << ":";
        ++idx;
      }
      std::vector<std::string> ylabels;
      std::transform(y_axis.begin(), y_axis.end(), std::back_inserter(ylabels), [](auto& bin) {
        return bin.first.label.empty() ? utils::format("%+g", bin.first.value) : bin.first.label;
      });
      struct stringlen {
        bool operator()(const std::string& a, const std::string& b) { return a.size() < b.size(); }
      };
      for (size_t i = 0; i < std::max_element(ylabels.begin(), ylabels.end(), stringlen())->size(); ++i) {
        os << "\n" << sep << ":";
        for (const auto& lab : ylabels)
          os << (lab.size() > i ? lab.at(i) : ' ');
        os << ":";
      }
      os << "\n"
         << sep << ":" << std::string(y_axis.size(), '.') << ": "  // 2nd abscissa axis
         << dr.yLabel() << "\n\t"
         << "(scale: \"" << std::string(CHARS) << "\", ";
      for (size_t i = 0; i < sizeof(kColours) / sizeof(kColours[0]); ++i)
        os << utils::colourise("*", (utils::Colour)kColours[i]) << (i == 0 ? "|" : "");
      os << ")\n";
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("text", DrawerText)
