#include "CepGen/Utils/Plotter.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Core/Exception.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_errno.h>

#include <iomanip>
#include <cmath>

namespace cepgen {
  namespace utils {
    Plotter::Hist::Hist() : log_(false) {}

    Plotter::Hist::Hist(const Hist& oth)
        : name_(oth.name_), xlabel_(oth.xlabel_), ylabel_(oth.ylabel_), log_(oth.log_) {}

    Plotter::Hist::~Hist() {}

    Plotter::Hist1D::Hist1D(size_t num_bins_x, const Limits& xrange) : underflow_(0ull), overflow_(0ull) {
      //info_.log = hvar.get<bool>("log", false);
      auto hist = gsl_histogram_alloc(num_bins_x);
      auto ret = gsl_histogram_set_ranges_uniform(hist, xrange.min(), xrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      CG_INFO("Plotter:Hist1D") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x) << " in range "
                                << xrange << ".";
    }

    Plotter::Hist1D::Hist1D(const Hist1D& oth)
        : Hist(oth),
          hist_(gsl_histogram_clone(oth.hist_.get())),
          underflow_(oth.underflow_),
          overflow_(oth.overflow_) {}

    void Plotter::Hist1D::fill(double x, double weight) {
      auto ret = gsl_histogram_accumulate(hist_.get(), x, weight);
      if (ret == GSL_SUCCESS)
        return;
      if (ret != GSL_EDOM)
        throw CG_FATAL("Hist1D:fill") << gsl_strerror(ret);
      const auto& rng = xrange();
      if (x < rng.min())
        underflow_++;
      else
        overflow_++;
    }

    void Plotter::Hist1D::add(Hist1D oth, double scaling) {
      oth.scale(scaling);
      auto ret = gsl_histogram_add(hist_.get(), oth.hist_.get());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:add") << gsl_strerror(ret);
    }

    void Plotter::Hist1D::scale(double scaling) {
      auto ret = gsl_histogram_scale(hist_.get(), scaling);
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:scale") << gsl_strerror(ret);
    }

    Limits Plotter::Hist1D::xrange() const {
      return Limits{gsl_histogram_min(hist_.get()), gsl_histogram_max(hist_.get())};
    }

    double Plotter::Hist1D::mean() const { return gsl_histogram_mean(hist_.get()); }
    double Plotter::Hist1D::rms() const { return gsl_histogram_sigma(hist_.get()); }

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
      const double bin_width = xrange().range() / nbins;
      os << "\n"
         << utils::format("%17s", name_.c_str()) << ":" << std::string(width, '.') << ":\n"  // 2nd abscissa axis
         << "\t("
         << "bin width=" << utils::s("unit", bin_width, true) << ", "
         << "mean=" << mean() << ", "
         << "st.dev.=" << rms();
      if (underflow_ > 0ull)
        os << ", underflow: " << underflow_;
      if (overflow_ > 0ull)
        os << ", overflow: " << overflow_;
      os << ")";
    }

    Plotter::Hist2D::Hist2D(size_t num_bins_x, const Limits& xrange, size_t num_bins_y, const Limits& yrange) {
      auto hist = gsl_histogram2d_alloc(num_bins_x, num_bins_y);
      auto ret = gsl_histogram2d_set_ranges_uniform(hist, xrange.min(), xrange.max(), yrange.min(), yrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D") << gsl_strerror(ret);
      hist_ = gsl_histogram2d_ptr(hist);
      CG_INFO("TextHandler") << "Booking a 2D correlation plot with " << utils::s("bin", num_bins_x + num_bins_y, true)
                             << " in ranges " << xrange << " and " << yrange << ".";
    }

    Plotter::Hist2D::Hist2D(const Hist2D& oth)
        : Hist(oth), hist_(gsl_histogram2d_clone(oth.hist_.get())), values_(oth.values_) {}

    void Plotter::Hist2D::fill(double x, double y, double weight) {
      auto ret = gsl_histogram2d_accumulate(hist_.get(), x, y, weight);
      if (ret == GSL_SUCCESS)
        return;
      if (ret != GSL_EDOM)
        throw CG_FATAL("Hist2D:fill") << gsl_strerror(ret);
      const auto &xrng = xrange(), &yrng = yrange();
      if (xrng.contains(x)) {
        if (y < yrng.min())
          values_.IN_LT++;
        else
          values_.IN_GT++;
      } else if (x < xrng.min()) {
        if (yrng.contains(y))
          values_.LT_IN++;
        else if (y < yrng.min())
          values_.LT_LT++;
        else
          values_.LT_GT++;
      } else {
        if (yrng.contains(y))
          values_.GT_IN++;
        else if (y < yrng.min())
          values_.GT_LT++;
        else
          values_.GT_GT++;
      }
    }

    void Plotter::Hist2D::add(Hist2D oth, double scaling) {
      oth.scale(scaling);
      auto ret = gsl_histogram2d_add(hist_.get(), oth.hist_.get());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D:add") << gsl_strerror(ret);
    }

    void Plotter::Hist2D::scale(double scaling) {
      auto ret = gsl_histogram2d_scale(hist_.get(), scaling);
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D:scale") << gsl_strerror(ret);
    }

    Limits Plotter::Hist2D::xrange() const {
      return Limits{gsl_histogram2d_xmin(hist_.get()), gsl_histogram2d_xmax(hist_.get())};
    }

    Limits Plotter::Hist2D::yrange() const {
      return Limits{gsl_histogram2d_ymin(hist_.get()), gsl_histogram2d_ymax(hist_.get())};
    }

    double Plotter::Hist2D::meanX() const { return gsl_histogram2d_xmean(hist_.get()); }
    double Plotter::Hist2D::rmsX() const { return gsl_histogram2d_xsigma(hist_.get()); }
    double Plotter::Hist2D::meanY() const { return gsl_histogram2d_ymean(hist_.get()); }
    double Plotter::Hist2D::rmsY() const { return gsl_histogram2d_ysigma(hist_.get()); }

    void Plotter::Hist2D::draw(std::ostream& os, size_t) const {
      const size_t nbins_x = gsl_histogram2d_nx(hist_.get());
      const size_t nbins_y = gsl_histogram2d_ny(hist_.get());
      const double max_bin = gsl_histogram2d_max_val(hist_.get());
      const std::string sep(17, ' ');
      if (!name_.empty())
        os << "plot of \"" << name_ << "\"\n";
      const auto x_range = xrange(), y_range = yrange();
      os << sep << std::string((size_t)std::max(0., 2. + nbins_y - ylabel_.size()), ' ') << ylabel_ << "\n"
         << sep << utils::format("%-5.2f", y_range.min()) << std::string(nbins_y - 11, ' ')
         << utils::format("%5.2e", y_range.max()) << "\n"
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
      const double bin_width_x = x_range.range() / nbins_x, bin_width_y = y_range.range() / nbins_y;
      os << "\n"
         << sep << ":" << std::string(nbins_y, '.') << ": "  // 2nd abscissa axis
         << ylabel_ << "\n\t"
         << "(scale: \"" << std::string(CHARS) << "\"\n\t"
         << " x-axis: "
         << "bin width=" << utils::s("unit", bin_width_x, true) << ", "
         << "mean=" << meanX() << ","
         << "st.dev.=" << rmsX() << "\n\t"
         << " y-axis: "
         << "bin width=" << utils::s("unit", bin_width_y, true) << ", "
         << "mean=" << meanY() << ","
         << "st.dev.=" << rmsY() << ")";
    }
  }  // namespace utils
}  // namespace cepgen
