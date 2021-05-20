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
    Plotter::Hist::Hist() : log_(false) {}

    Plotter::Hist::Hist(const Hist& oth)
        : name_(oth.name_), xlabel_(oth.xlabel_), ylabel_(oth.ylabel_), log_(oth.log_) {}

    Plotter::Hist::~Hist() {}

    Plotter::Hist1D::Hist1D(size_t num_bins_x, double min_x, double max_x) {
      //info_.log = hvar.get<bool>("log", false);
      auto hist = gsl_histogram_alloc(num_bins_x);
      gsl_histogram_set_ranges_uniform(hist, min_x, max_x);
      hist_ = gsl_histogram_ptr(hist);
      CG_INFO("Plotter:Hist1D") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x) << " between " << min_x
                                << " and " << max_x << ".";
    }

    Plotter::Hist1D::Hist1D(const Hist1D& oth) : Hist(oth), hist_(gsl_histogram_clone(oth.hist_.get())) {}

    void Plotter::Hist1D::fill(double x, double weight) { gsl_histogram_accumulate(hist_.get(), x, weight); }

    void Plotter::Hist1D::scale(double scaling) { gsl_histogram_scale(hist_.get(), scaling); }

    void Plotter::Hist1D::draw(std::ostream& os, size_t width) const {
      const size_t nbins = gsl_histogram_bins(hist_.get());
      const double max_bin = gsl_histogram_max_val(hist_.get());
      const double min_bin = gsl_histogram_min_val(hist_.get());
      const double min_range_log = std::log(std::max(min_bin, 1.e-10));
      const double max_range_log = std::log(std::min(max_bin, 1.e+10));
      const std::string sep(17, ' ');
      if (!name_.empty())
        os << "plot of \"" << name_ << "\"\n";
      os << sep << std::string(std::max(0., 2. + width - ylabel_.size()), ' ') << ylabel_ << "\n"
         << sep << utils::format("%-5.2f", log_ ? std::exp(min_range_log) : min_bin) << std::setw(width - 11)
         << std::left << (log_ ? "logarithmic scale" : "linear scale")
         << utils::format("%5.2e", log_ ? std::exp(max_range_log) : max_bin) << "\n"
         << sep << std::string(width + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins; ++i) {
        double min, max;
        gsl_histogram_get_range(hist_.get(), i, &min, &max);
        const double value = gsl_histogram_get(hist_.get(), i), unc = sqrt(value);
        size_t val = 0ull;
        {
          double val_dbl = width;
          if (log_)
            val_dbl *= (value > 0. && max_bin > 0.)
                           ? std::max((std::log(value) - min_range_log) / (max_range_log - min_range_log), 0.)
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
         << utils::format("%17s", name_.c_str()) << ":" << std::string(width, '.') << ":\n"  // 2nd abscissa axis
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

    Plotter::Hist2D::Hist2D(const Hist2D& oth) : Hist(oth), hist_(gsl_histogram2d_clone(oth.hist_.get())) {}

    void Plotter::Hist2D::fill(double x, double y, double weight) {
      gsl_histogram2d_accumulate(hist_.get(), x, y, weight);
    }

    void Plotter::Hist2D::scale(double scaling) { gsl_histogram2d_scale(hist_.get(), scaling); }

    void Plotter::Hist2D::draw(std::ostream& os, size_t) const {
      const size_t nbins_x = gsl_histogram2d_nx(hist_.get());
      const size_t nbins_y = gsl_histogram2d_ny(hist_.get());
      const double max_bin = gsl_histogram2d_max_val(hist_.get());
      const std::string sep(17, ' ');
      if (!name_.empty())
        os << "plot of \"" << name_ << "\"\n";
      os << sep << std::string((size_t)std::max(0., 2. + nbins_y - ylabel_.size()), ' ') << ylabel_ << "\n"
         << sep << utils::format("%-5.2f", gsl_histogram2d_ymin(hist_.get())) << std::string(nbins_y - 11, ' ')
         << utils::format("%5.2e", gsl_histogram2d_ymax(hist_.get())) << "\n"
         << utils::format("%17s", xlabel_.c_str()) << std::string(nbins_y + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins_x; ++i) {
        double min_x, max_x;
        gsl_histogram2d_get_xrange(hist_.get(), i, &min_x, &max_x);
        os << "\n" << utils::format("[%7.2f,%7.2f):", min_x, max_x);
        for (size_t j = 0; j < nbins_y; ++j) {
          const double value = gsl_histogram2d_get(hist_.get(), i, j);
          const double value_norm = log_ ? (value == 0. ? 0. : std::log(value) / std::log(max_bin)) : value / max_bin;
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
         << ylabel_ << "\n\t"
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
