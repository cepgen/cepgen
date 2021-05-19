#include "CepGen/Utils/Plotter.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Core/Exception.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <iomanip>
#include <cmath>
#include <fstream>

namespace cepgen {
  namespace utils {
    Plotter::Hist1D::Hist1D(size_t num_bins_x, double min_x, double max_x) {
      //info_.log = hvar.get<bool>("log", false);
      auto hist = gsl_histogram_alloc(num_bins_x);
      gsl_histogram_set_ranges_uniform(hist, min_x, max_x);
      hist_ = gsl_histogram_ptr(hist);
      CG_INFO("Plotter:Hist1D") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x) << " between " << min_x
                                << " and " << max_x << ".";
    }

    void Plotter::Hist1D::fill(double x, double weight) { gsl_histogram_accumulate(hist_.get(), x, weight); }

    void Plotter::Hist1D::draw(std::ostream& os, size_t width) const {
      const size_t nbins = gsl_histogram_bins(hist_.get());
      const double max_bin = gsl_histogram_max_val(hist_.get());
      const double min_bin = gsl_histogram_min_val(hist_.get());
      const double min_range_log = std::log(std::max(min_bin, 1.e-10));
      const double max_range_log = std::log(std::min(max_bin, 1.e+10));
      const std::string sep(17, ' ');
      const auto& var = info_.name.at(0);
      os << "plot of \"" << var << "\"\n"
         << sep << std::string(width - 15 - var.size(), ' ') << info_.xlabel << "\n"
         << sep << utils::format("%-5.2f", info_.log ? exp(min_range_log) : min_bin) << std::setw(width - 11)
         << std::left << (info_.log ? "logarithmic scale" : "linear scale")
         << utils::format("%5.2e", info_.log ? exp(max_range_log) : max_bin) << "\n"
         << sep << std::string(width + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins; ++i) {
        double min, max;
        gsl_histogram_get_range(hist_.get(), i, &min, &max);
        const double value = gsl_histogram_get(hist_.get(), i), unc = sqrt(value);
        size_t val = 0ull;
        {
          double val_dbl = width;
          if (info_.log)
            val_dbl *= (value > 0. && max_bin > 0.)
                           ? std::max((log(value) - min_range_log) / (max_range_log - min_range_log), 0.)
                           : 0.;
          else if (max_bin > 0.)
            val_dbl *= (value > 0. && max_bin > 0.) ? value / max_bin : 0.;
          val = std::ceil(val_dbl);
        }
        os << "\n"
           << utils::format("[%7.2f,%7.2f):", min, max) << std::string(val, CHAR) << std::string(width - val, ' ')
           << ": " << utils::format("%6.2e", value) << " +/- " << utils::format("%6.2e", unc);
      }
      const double bin_width = (gsl_histogram_max(hist_.get()) - gsl_histogram_min(hist_.get())) / nbins;
      os << "\n"
         << utils::format("%17s", var.c_str()) << ":" << std::string(width, '.') << ":\n"  // 2nd abscissa axis
         << "\t("
         << "bin width=" << bin_width << utils::s(" unit", (int)bin_width, false) << ", "
         << "mean=" << gsl_histogram_mean(hist_.get()) << ", "
         << "st.dev.=" << gsl_histogram_sigma(hist_.get()) << ")";
    }

    Plotter::Hist2D::Hist2D(
        size_t num_bins_x, double min_x, double max_x, size_t num_bins_y, double min_y, double max_y) {
      auto hist = gsl_histogram2d_alloc(num_bins_x, num_bins_y);
      gsl_histogram2d_set_ranges_uniform(hist, min_x, max_x, min_y, max_y);
      hist_ = gsl_histogram2d_ptr(hist);
      CG_INFO("TextHandler") << "Booking a 2D correlation plot with " << utils::s("bin", num_bins_x + num_bins_y)
                             << " between (" << min_x << ", " << max_x << ") and (" << min_y << ", " << max_y << ").";
    }

    void Plotter::Hist2D::fill(double x, double y, double weight) {
      gsl_histogram2d_accumulate(hist_.get(), x, y, weight);
    }

    void Plotter::Hist2D::draw(std::ostream& os, size_t width) const {
      const size_t nbins_x = gsl_histogram2d_nx(hist_.get());
      const size_t nbins_y = gsl_histogram2d_ny(hist_.get());
      const double max_bin = gsl_histogram2d_max_val(hist_.get());
      const std::string sep(17, ' ');
      const auto& vars = info_.name;
      const auto var = utils::merge(vars, "/");
      os << "plot of \"" << info_.xlabel << "\"\n"
         << sep << std::string((size_t)std::max(0., nbins_y - 15. - var.size()), ' ');
      os << info_.ylabel << "\n"
         << sep << utils::format("%-5.2f", gsl_histogram2d_ymin(hist_.get())) << std::string(nbins_y - 11, ' ')
         << utils::format("%5.2e", gsl_histogram2d_ymax(hist_.get())) << "\n"
         << utils::format("%17s", vars.at(0).c_str()) << std::string(nbins_y + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins_x; ++i) {
        double min_x, max_x;
        gsl_histogram2d_get_xrange(hist_.get(), i, &min_x, &max_x);
        os << "\n" << utils::format("[%7.2f,%7.2f):", min_x, max_x);
        for (size_t j = 0; j < nbins_y; ++j) {
          const double value = gsl_histogram2d_get(hist_.get(), i, j);
          const double value_norm = info_.log ? (value == 0. ? 0. : log(value) / log(max_bin)) : value / max_bin;
          os << CHARS[(size_t)ceil(value_norm * (strlen(CHARS) - 1))];
        }
        os << ":";
      }
      std::vector<std::string> ylabels;
      for (size_t j = 0; j < nbins_y; ++j) {
        double min_y, max_y;
        gsl_histogram2d_get_yrange(hist_.get(), j, &min_y, &max_y);
        ylabels.emplace_back(utils::format("%+g", min_y));
      }
      struct stringlen {
        bool operator()(const std::string& a, const std::string& b) { return a.size() < b.size(); }
      };
      for (size_t i = 0; i < std::max_element(ylabels.begin(), ylabels.end(), stringlen())->size(); ++i) {
        os << "\n" << sep << ":";
        for (const auto& lab : ylabels)
          os << (lab.size() > i ? lab.at(i) : ' ');
        os << ":";
      }
      const double bin_width_x = (gsl_histogram2d_xmax(hist_.get()) - gsl_histogram2d_xmin(hist_.get())) / nbins_x;
      const double bin_width_y = (gsl_histogram2d_ymax(hist_.get()) - gsl_histogram2d_ymin(hist_.get())) / nbins_y;
      os << "\n"
         << sep << ":" << std::string(nbins_y, '.') << ": "  // 2nd abscissa axis
         << vars.at(1) << "\n\t"
         << "(scale: \"" << std::string(CHARS) << "\"\n\t"
         << " x-axis: "
         << "bin width=" << bin_width_x << utils::s(" unit", (int)bin_width_x, false) << ", "
         << "mean=" << gsl_histogram2d_xmean(hist_.get()) << ","
         << "st.dev.=" << gsl_histogram2d_xsigma(hist_.get()) << "\n\t"
         << " y-axis: "
         << "bin width=" << bin_width_y << utils::s(" unit", (int)bin_width_y, false) << ", "
         << "mean=" << gsl_histogram2d_ymean(hist_.get()) << ","
         << "st.dev.=" << gsl_histogram2d_ysigma(hist_.get()) << ")";
    }
  }  // namespace utils
}  // namespace cepgen
