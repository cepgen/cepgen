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

#include <gsl/gsl_errno.h>

#include <cmath>
#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Hist1D::Hist1D(size_t num_bins_x, const Limits& xrange, const std::string& name, const std::string& title)
        : Drawable(name, title) {
      auto hist = gsl_histogram_alloc(num_bins_x);
      auto ret = gsl_histogram_set_ranges_uniform(hist, xrange.min(), xrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
      CG_INFO("Hist1D") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x, true) << " in range " << xrange
                        << ".";
    }

    Hist1D::Hist1D(const std::vector<double>& xbins, const std::string& name, const std::string& title)
        : Drawable(name, title) {
      auto hist = gsl_histogram_alloc(xbins.size() - 1);
      auto ret = gsl_histogram_set_ranges(hist, xbins.data(), xbins.size());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
      CG_INFO("Hist1D") << "Booking a 1D histogram with " << utils::s("bin", xbins.size(), true) << " in range "
                        << xbins << ".";
    }

    Hist1D::Hist1D(const Hist1D& oth)
        : Histogram(oth),
          Drawable(oth),
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
        if (auto ret2 = gsl_histogram_accumulate(hist_w2_.get(), x, weight * weight) != GSL_SUCCESS)
          throw CG_FATAL("Hist1D:fill") << "(w2 histogram): " << gsl_strerror(ret2);
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

    Drawable::axis_t Hist1D::axis() const {
      axis_t axis;
      for (size_t bin = 0; bin < nbins(); ++bin) {
        const auto& range_i = binRange(bin);
        axis[coord_t{
            range_i.x(0.5), 0.5 * range_i.range(), utils::format("[%7.2f,%7.2f)", range_i.min(), range_i.max())}] =
            value_t{value(bin), valueUnc(bin)};
      }
      return axis;
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

    Hist2D::Hist2D(size_t num_bins_x,
                   const Limits& xrange,
                   size_t num_bins_y,
                   const Limits& yrange,
                   const std::string& name,
                   const std::string& title)
        : Drawable(name, title) {
      auto hist = gsl_histogram2d_alloc(num_bins_x, num_bins_y);
      auto ret = gsl_histogram2d_set_ranges_uniform(hist, xrange.min(), xrange.max(), yrange.min(), yrange.max());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D") << gsl_strerror(ret);
      hist_ = gsl_histogram2d_ptr(hist);
      hist_w2_ = gsl_histogram2d_ptr(gsl_histogram2d_clone(hist_.get()));
      CG_INFO("Hist2D") << "Booking a 2D correlation plot with " << utils::s("bin", num_bins_x + num_bins_y, true)
                        << " in ranges " << xrange << " and " << yrange << ".";
    }

    Hist2D::Hist2D(const std::vector<double>& xbins,
                   const std::vector<double>& ybins,
                   const std::string& name,
                   const std::string& title)
        : Drawable(name, title) {
      auto hist = gsl_histogram2d_alloc(xbins.size() - 1, ybins.size() - 1);
      auto ret = gsl_histogram2d_set_ranges(hist, xbins.data(), xbins.size(), ybins.data(), ybins.size());
      if (ret != GSL_SUCCESS)
        throw CG_FATAL("Hist2D") << gsl_strerror(ret);
      hist_ = gsl_histogram2d_ptr(hist);
      hist_w2_ = gsl_histogram2d_ptr(gsl_histogram2d_clone(hist_.get()));
      CG_INFO("Hist2D") << "Booking a 2D correlation plot with " << utils::s("bin", xbins.size() + ybins.size(), true)
                        << " in ranges x=(" << xbins << ") and y=" << ybins << ".";
    }

    Hist2D::Hist2D(const Hist2D& oth)
        : Histogram(oth),
          Drawable(oth),
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
  }  // namespace utils
}  // namespace cepgen
