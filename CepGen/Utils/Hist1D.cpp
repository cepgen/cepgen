/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Hist1D::Hist1D(size_t num_bins_x, const Limits& xrange, const std::string& name, const std::string& title)
        : Drawable(name, title) {
      if (num_bins_x == 0)
        throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
      auto hist = gsl_histogram_alloc(num_bins_x);
      if (auto ret = gsl_histogram_set_ranges_uniform(hist, xrange.min(), xrange.max()); ret != GSL_SUCCESS)
        throw CG_ERROR("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
      CG_DEBUG("Hist1D") << "Booking a 1D histogram with " << s("bin", num_bins_x, true) << " in range " << xrange
                         << ".";
    }

    Hist1D::Hist1D(const std::vector<double>& xbins, const std::string& name, const std::string& title)
        : Drawable(name, title) {
      if (xbins.empty())
        throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
      auto hist = gsl_histogram_alloc(xbins.size() - 1);
      if (auto ret = gsl_histogram_set_ranges(hist, xbins.data(), xbins.size()); ret != GSL_SUCCESS)
        throw CG_ERROR("Hist1D") << gsl_strerror(ret);
      hist_ = gsl_histogram_ptr(hist);
      hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
      CG_DEBUG("Hist1D") << "Booking a 1D histogram with " << s("bin", xbins.size(), true) << " in range " << xbins
                         << ".";
    }

    Hist1D::Hist1D(const Hist1D& oth)
        : Histogram(oth),
          Drawable(oth),
          hist_(gsl_histogram_clone(oth.hist_.get())),
          hist_w2_(gsl_histogram_clone(oth.hist_w2_.get())),
          underflow_(oth.underflow_),
          overflow_(oth.overflow_) {}

    void Hist1D::clear() {
      CG_ASSERT(hist_);
      CG_ASSERT(hist_w2_);
      gsl_histogram_reset(hist_.get());
      gsl_histogram_reset(hist_w2_.get());
    }

    void Hist1D::fill(double x, double weight) {
      CG_ASSERT(hist_);
      {  // reduce the scope of 'ret'
        auto ret = gsl_histogram_accumulate(hist_.get(), x, weight);
        if (ret == GSL_SUCCESS) {
          if (auto ret2 = gsl_histogram_accumulate(hist_w2_.get(), x, weight * weight) != GSL_SUCCESS)
            throw CG_ERROR("Hist1D:fill") << "(w2 histogram): " << gsl_strerror(ret2);
          return;
        }
        if (ret != GSL_EDOM)
          throw CG_ERROR("Hist1D:fill") << gsl_strerror(ret);
      }
      if (x < range().min())
        underflow_ += weight;
      else
        overflow_ += weight;
    }

    void Hist1D::add(Hist1D oth, double scaling) {
      CG_ASSERT(hist_);
      CG_ASSERT(hist_w2_);
      CG_ASSERT(oth.hist_);
      CG_ASSERT(oth.hist_w2_);
      if (oth.integral(true) == 0.) {
        CG_WARNING("Hist1D:add") << "Other histogram is empty.";
        return;
      }
      const double scl = std::pow(oth.integral(), -2);
      oth.scale(scaling);
      gsl_histogram_scale(oth.hist_w2_.get(), scl);
      if (auto ret = gsl_histogram_add(hist_.get(), oth.hist_.get()); ret != GSL_SUCCESS)
        throw CG_ERROR("Hist1D:add") << gsl_strerror(ret);
      gsl_histogram_add(hist_w2_.get(), oth.hist_w2_.get());
      underflow_ += oth.underflow_;
      overflow_ += oth.overflow_;
    }

    void Hist1D::scale(double scaling) {
      CG_ASSERT(hist_);
      if (auto ret = gsl_histogram_scale(hist_.get(), scaling); ret != GSL_SUCCESS)
        throw CG_ERROR("Hist1D:scale") << gsl_strerror(ret);
      gsl_histogram_scale(hist_w2_.get(), scaling * scaling);
      underflow_ *= scaling;
      overflow_ *= scaling;
    }

    Drawable::axis_t Hist1D::axis() const {
      axis_t axis;
      for (size_t bin = 0; bin < nbins(); ++bin) {
        const auto& range_i = binRange(bin);
        axis[coord_t{range_i.x(0.5), 0.5 * range_i.range(), format("[%7.2f,%7.2f)", range_i.min(), range_i.max())}] =
            value(bin);
      }
      return axis;
    }

    size_t Hist1D::nbins() const {
      CG_ASSERT(hist_);
      return gsl_histogram_bins(hist_.get());
    }

    Limits Hist1D::range() const {
      CG_ASSERT(hist_);
      return Limits{gsl_histogram_min(hist_.get()), gsl_histogram_max(hist_.get())};
    }

    Limits Hist1D::binRange(size_t bin) const {
      CG_ASSERT(hist_);
      Limits range;
      if (auto ret = gsl_histogram_get_range(hist_.get(), bin, &range.min(), &range.max()); ret != GSL_SUCCESS)
        throw CG_ERROR("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
      return range;
    }

    std::vector<double> Hist1D::bins(BinMode mode) const {
      const auto bins = extractBins(mode, nbins(), std::bind(&Hist1D::binRange, this, std::placeholders::_1));
      return std::vector<double>(bins.begin(), bins.end());
    }

    std::vector<Value> Hist1D::values() const {
      std::vector<Value> values;
      for (size_t i = 0; i < nbins(); ++i)
        values.emplace_back(value(i));
      return values;
    }

    Value Hist1D::value(size_t bin) const {
      CG_ASSERT(hist_);
      CG_ASSERT(hist_w2_);
      return Value{gsl_histogram_get(hist_.get(), bin), std::sqrt(gsl_histogram_get(hist_w2_.get(), bin))};
    }

    double Hist1D::mean() const {
      CG_ASSERT(hist_);
      return gsl_histogram_mean(hist_.get());
    }

    double Hist1D::rms() const {
      CG_ASSERT(hist_);
      return gsl_histogram_sigma(hist_.get());
    }

    double Hist1D::minimum() const {
      CG_ASSERT(hist_);
      return gsl_histogram_min_val(hist_.get());
    }

    double Hist1D::maximum() const {
      CG_ASSERT(hist_);
      return gsl_histogram_max_val(hist_.get());
    }

    double Hist1D::integral(bool include_out_of_range) const {
      CG_ASSERT(hist_);
      auto integr = gsl_histogram_sum(hist_.get());
      if (include_out_of_range)
        integr += underflow_ + overflow_;
      return integr;
    }
  }  // namespace utils
}  // namespace cepgen
