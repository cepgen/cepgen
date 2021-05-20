#include "CepGen/Utils/Plotter.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Core/Exception.h"

#include <gsl/gsl_errno.h>

#include <iomanip>
#include <cmath>

namespace cepgen {
  namespace utils {
    Hist::Hist() : log_(false) {}

    Hist::Hist(const Hist& oth) : name_(oth.name_), xlabel_(oth.xlabel_), ylabel_(oth.ylabel_), log_(oth.log_) {}

    Hist::~Hist() {}

    Hist1D::Hist1D(size_t num_bins_x, const Limits& xrange) : underflow_(0ull), overflow_(0ull) {
      //info_.log = hvar.get<bool>("log", false);
      auto hist = gsl_histogram_alloc(num_bins_x);
      auto ret = gsl_histogram_set_ranges_uniform(hist, xrange.min(), xrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
      CG_INFO("Plotter:Hist1D") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x, true) << " in range "
                                << xrange << ".";
    }

    Hist1D::Hist1D(const Hist1D& oth)
        : Hist(oth),
          hist_(gsl_histogram_clone(oth.hist_.get())),
          hist_w2_(gsl_histogram_clone(oth.hist_w2_.get())),
          underflow_(oth.underflow_),
          overflow_(oth.overflow_) {}

    void Hist1D::clear() {
      gsl_histogram_reset(hist_.get());
      gsl_histogram_reset(hist_w2_.get());
    }

    void Hist1D::fill(double x, double weight) {
      auto ret = gsl_histogram_accumulate(hist_.get(), x, weight);
      if (ret == GSL_SUCCESS) {
        gsl_histogram_accumulate(hist_w2_.get(), x, weight * weight);
        return;
      }
      if (ret != GSL_EDOM)
        throw CG_FATAL("Hist1D:fill") << gsl_strerror(ret);
      if (x < range().min())
        underflow_ += weight;
      else
        overflow_ += weight;
    }

    void Hist1D::add(Hist1D oth, double scaling) {
      if (oth.integral() == 0.) {
        CG_WARNING("Hist1D:add") << "Other histogram is empty.";
        return;
      }
      const double scl = std::pow(oth.integral(), -2);
      oth.scale(scaling);
      gsl_histogram_scale(oth.hist_w2_.get(), scl);
      auto ret = gsl_histogram_add(hist_.get(), oth.hist_.get());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:add") << gsl_strerror(ret);
      gsl_histogram_add(hist_w2_.get(), oth.hist_w2_.get());
    }

    void Hist1D::scale(double scaling) {
      auto ret = gsl_histogram_scale(hist_.get(), scaling);
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:scale") << gsl_strerror(ret);
      gsl_histogram_scale(hist_w2_.get(), scaling * scaling);
    }

    size_t Hist1D::nbins() const { return gsl_histogram_bins(hist_.get()); }
    Limits Hist1D::range() const { return Limits{gsl_histogram_min(hist_.get()), gsl_histogram_max(hist_.get())}; }
    Limits Hist1D::binRange(size_t bin) const {
      Limits range;
      auto ret = gsl_histogram_get_range(hist_.get(), bin, &range.min(), &range.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
      return range;
    }

    double Hist1D::value(size_t bin) const { return gsl_histogram_get(hist_.get(), bin); }
    double Hist1D::valueUnc(size_t bin) const { return std::sqrt(gsl_histogram_get(hist_w2_.get(), bin)); }

    double Hist1D::mean() const { return gsl_histogram_mean(hist_.get()); }
    double Hist1D::rms() const { return gsl_histogram_sigma(hist_.get()); }
    double Hist1D::minimum() const { return gsl_histogram_min_val(hist_.get()); }
    double Hist1D::maximum() const { return gsl_histogram_max_val(hist_.get()); }
    double Hist1D::integral() const { return gsl_histogram_sum(hist_.get()); }

    void Hist1D::draw(std::ostream& os, size_t width) const {
      const double max_val = maximum() * (log_ ? 5. : 1.2), min_val = minimum();
      const double min_range_log = std::log(std::max(min_val, 1.e-10));
      const double max_range_log = std::log(std::min(max_val, 1.e+10));
      const std::string sep(17, ' ');
      if (!name_.empty())
        os << "plot of \"" << name_ << "\"\n";
      os << sep << std::string(std::max(0., 2. + width - ylabel_.size()), ' ') << ylabel_ << "\n"
         << sep << utils::format("%-5.2f", log_ ? std::exp(min_range_log) : min_val) << std::setw(width - 11)
         << std::left << (log_ ? "logarithmic scale" : "linear scale")
         << utils::format("%5.2e", log_ ? std::exp(max_range_log) : max_val) << "\n"
         << sep << std::string(width + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins(); ++i) {
        const auto range_i = binRange(i);
        const double val = value(i), unc = valueUnc(i);
        size_t ival = 0ull, ierr = 0ull;
        {
          double val_dbl = width, unc_dbl = width;
          if (log_) {
            val_dbl *= (val > 0. && max_val > 0.)
                           ? std::max((std::log(val) - min_range_log) / (max_range_log - min_range_log), 0.)
                           : 0.;
            unc_dbl *= (val > 0. && max_val > 0.)
                           ? std::max((std::log(unc) - min_range_log) / (max_range_log - min_range_log), 0.)
                           : 0.;
          } else if (max_val > 0.) {
            val_dbl *= (val > 0. && max_val > 0.) ? val / max_val : 0.;
            unc_dbl *= (unc > 0. && max_val > 0.) ? unc / max_val : 0.;
          }
          ival = std::ceil(val_dbl);
          ierr = std::ceil(unc_dbl);
        }
        os << "\n"
           << utils::format("[%7.2f,%7.2f):", range_i.min(), range_i.max())
           << (ival > ierr ? std::string(ival - ierr, ' ') : "") << (ierr > 0 ? std::string(ierr, ERR_CHAR) : "")
           << CHAR << (ierr > 0 ? std::string(std::min(width - ival - 1, ierr), ERR_CHAR) : "")
           << (ival + ierr < width ? std::string(width - ival - ierr - 1, ' ') : "") << ": "
           << utils::format("%6.2e +/- %6.2e", val, unc);
      }
      const double bin_width = range().range() / nbins();
      os << "\n"
         << utils::format("%17s", name_.c_str()) << ":" << std::string(width, '.') << ":\n"  // 2nd abscissa axis
         << "\t("
         << "bin width=" << utils::s("unit", bin_width, true) << ", "
         << "mean=" << mean() << ", "
         << "st.dev.=" << rms() << ", "
         << "integr.=" << integral();
      if (underflow_ > 0ull)
        os << ", underflow: " << underflow_;
      if (overflow_ > 0ull)
        os << ", overflow: " << overflow_;
      os << ")";
    }

    Hist2D::Hist2D(size_t num_bins_x, const Limits& xrange, size_t num_bins_y, const Limits& yrange) {
      auto hist = gsl_histogram2d_alloc(num_bins_x, num_bins_y);
      auto ret = gsl_histogram2d_set_ranges_uniform(hist, xrange.min(), xrange.max(), yrange.min(), yrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D") << gsl_strerror(ret);
      hist_ = gsl_histogram2d_ptr(hist);
      hist_w2_ = gsl_histogram2d_ptr(gsl_histogram2d_clone(hist_.get()));
      CG_INFO("TextHandler") << "Booking a 2D correlation plot with " << utils::s("bin", num_bins_x + num_bins_y, true)
                             << " in ranges " << xrange << " and " << yrange << ".";
    }

    Hist2D::Hist2D(const Hist2D& oth)
        : Hist(oth),
          hist_(gsl_histogram2d_clone(oth.hist_.get())),
          hist_w2_(gsl_histogram2d_clone(oth.hist_w2_.get())),
          values_(oth.values_) {}

    void Hist2D::clear() {
      gsl_histogram2d_reset(hist_.get());
      gsl_histogram2d_reset(hist_w2_.get());
    }

    void Hist2D::fill(double x, double y, double weight) {
      auto ret = gsl_histogram2d_accumulate(hist_.get(), x, y, weight);
      if (ret == GSL_SUCCESS) {
        gsl_histogram2d_accumulate(hist_w2_.get(), x, y, weight * weight);
        return;
      }
      if (ret != GSL_EDOM)
        throw CG_FATAL("Hist2D:fill") << gsl_strerror(ret);
      const auto &xrng = rangeX(), &yrng = rangeY();
      if (xrng.contains(x)) {
        if (y < yrng.min())
          values_.IN_LT += weight;
        else
          values_.IN_GT += weight;
      } else if (x < xrng.min()) {
        if (yrng.contains(y))
          values_.LT_IN += weight;
        else if (y < yrng.min())
          values_.LT_LT += weight;
        else
          values_.LT_GT += weight;
      } else {
        if (yrng.contains(y))
          values_.GT_IN += weight;
        else if (y < yrng.min())
          values_.GT_LT += weight;
        else
          values_.GT_GT += weight;
      }
    }

    void Hist2D::add(Hist2D oth, double scaling) {
      if (oth.integral() == 0.) {
        CG_WARNING("Hist1D:add") << "Other histogram is empty.";
        return;
      }
      const double scl = std::pow(oth.integral(), -2);
      oth.scale(scaling);
      gsl_histogram2d_scale(oth.hist_w2_.get(), scl);
      auto ret = gsl_histogram2d_add(hist_.get(), oth.hist_.get());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D:add") << gsl_strerror(ret);
      gsl_histogram2d_add(hist_w2_.get(), oth.hist_w2_.get());
    }

    void Hist2D::scale(double scaling) {
      auto ret = gsl_histogram2d_scale(hist_.get(), scaling);
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D:scale") << gsl_strerror(ret);
      gsl_histogram2d_scale(hist_w2_.get(), scaling * scaling);
    }

    size_t Hist2D::nbinsX() const { return gsl_histogram2d_nx(hist_.get()); }
    Limits Hist2D::rangeX() const {
      return Limits{gsl_histogram2d_xmin(hist_.get()), gsl_histogram2d_xmax(hist_.get())};
    }
    Limits Hist2D::binRangeX(size_t bin) const {
      Limits range;
      auto ret = gsl_histogram2d_get_xrange(hist_.get(), bin, &range.min(), &range.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
      return range;
    }

    size_t Hist2D::nbinsY() const { return gsl_histogram2d_ny(hist_.get()); }
    Limits Hist2D::rangeY() const {
      return Limits{gsl_histogram2d_ymin(hist_.get()), gsl_histogram2d_ymax(hist_.get())};
    }
    Limits Hist2D::binRangeY(size_t bin) const {
      Limits range;
      auto ret = gsl_histogram2d_get_yrange(hist_.get(), bin, &range.min(), &range.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
      return range;
    }

    double Hist2D::value(size_t bin_x, size_t bin_y) const { return gsl_histogram2d_get(hist_.get(), bin_x, bin_y); }
    double Hist2D::valueUnc(size_t bin_x, size_t bin_y) const {
      return std::sqrt(gsl_histogram2d_get(hist_w2_.get(), bin_x, bin_y));
    }

    double Hist2D::meanX() const { return gsl_histogram2d_xmean(hist_.get()); }
    double Hist2D::rmsX() const { return gsl_histogram2d_xsigma(hist_.get()); }
    double Hist2D::meanY() const { return gsl_histogram2d_ymean(hist_.get()); }
    double Hist2D::rmsY() const { return gsl_histogram2d_ysigma(hist_.get()); }
    double Hist2D::minimum() const { return gsl_histogram2d_min_val(hist_.get()); }
    double Hist2D::maximum() const { return gsl_histogram2d_max_val(hist_.get()); }
    double Hist2D::integral() const { return gsl_histogram2d_sum(hist_.get()); }

    void Hist2D::draw(std::ostream& os, size_t) const {
      const size_t nbins_x = gsl_histogram2d_nx(hist_.get());
      const size_t nbins_y = gsl_histogram2d_ny(hist_.get());
      const double max_val = maximum();
      const std::string sep(17, ' ');
      if (!name_.empty())
        os << "plot of \"" << name_ << "\"\n";
      const auto x_range = rangeX(), y_range = rangeY();
      os << sep << std::string((size_t)std::max(0., 2. + nbins_y - ylabel_.size()), ' ') << ylabel_ << "\n"
         << sep << utils::format("%-5.2f", y_range.min()) << std::string(nbins_y - 11, ' ')
         << utils::format("%5.2e", y_range.max()) << "\n"
         << utils::format("%17s", xlabel_.c_str()) << std::string(nbins_y + 2, '.');  // abscissa axis
      for (size_t i = 0; i < nbins_x; ++i) {
        const auto& xrange_i = binRangeX(i);
        os << "\n" << utils::format("[%7.2f,%7.2f):", xrange_i.min(), xrange_i.max());
        for (size_t j = 0; j < nbins_y; ++j) {
          const double val = value(i, j);
          const double val_norm = log_ ? (val == 0. ? 0. : std::log(val) / std::log(max_val)) : val / max_val;
          os << CHARS[(size_t)ceil(val_norm * (strlen(CHARS) - 1))];
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
         << "st.dev.=" << rmsY() << ",\n\t"
         << " integr.=" << integral() << ")";
    }
  }  // namespace utils
}  // namespace cepgen
