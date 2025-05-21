/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::utils;

Hist2D::Hist2D(const ParametersList& params)
    : Drawable(params.get<std::string>("name"), params.get<std::string>("title")) {
  const auto &x_bins = params.get<std::vector<double> >("xbins"), &y_bins = params.get<std::vector<double> >("ybins");
  const auto &x_range = params.get<Limits>("xrange"), &y_range = params.get<Limits>("yrange");
  const auto &num_bins_x = params.get<int>("nbinsX"), &num_bins_y = params.get<int>("nbinsY");
  if (x_bins.size() > 1 && y_bins.size() > 1)
    buildFromBins(x_bins, y_bins);
  else if (x_range.valid() && y_range.valid() && num_bins_x > 1 && num_bins_y > 1)
    buildFromRange(num_bins_x, x_range, num_bins_y, y_range);
  else
    throw CG_FATAL("Hist2D") << "Failed to build a 2D histogram with user parameters: " << params << ".";
}

Hist2D::Hist2D(size_t num_bins_x,
               const Limits& xrange,
               size_t num_bins_y,
               const Limits& y_range,
               const std::string& name,
               const std::string& title)
    : Drawable(name, title) {
  if (num_bins_x == 0 || num_bins_y == 0)
    throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
  buildFromRange(num_bins_x, xrange, num_bins_y, y_range);
}

Hist2D::Hist2D(const std::vector<double>& x_bins,
               const std::vector<double>& y_bins,
               const std::string& name,
               const std::string& title)
    : Drawable(name, title) {
  if (x_bins.empty() || y_bins.empty())
    throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
  buildFromBins(x_bins, y_bins);
}

Hist2D::Hist2D(const Hist2D& oth)
    : Histogram(oth),
      Drawable(oth),
      hist_(gsl_histogram2d_clone(oth.hist_.get())),
      hist_w2_(gsl_histogram2d_clone(oth.hist_w2_.get())),
      out_of_range_values_(oth.out_of_range_values_) {}

void Hist2D::buildFromBins(const std::vector<double>& x_bins, const std::vector<double>& y_bins) {
  if (x_bins.size() < 1 || y_bins.size() < 1)
    throw CG_ERROR("Hist2D:buildFromBins") << "Building a 2D histogram requires at least 1x1 bin.";
  hist_.reset(gsl_histogram2d_alloc(x_bins.size() - 1, y_bins.size() - 1));
  CG_ASSERT(hist_);
  if (auto ret = gsl_histogram2d_set_ranges(hist_.get(), x_bins.data(), x_bins.size(), y_bins.data(), y_bins.size());
      ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:buildFromBins") << gsl_strerror(ret);
  hist_w2_ = gsl_histogram2d_ptr(gsl_histogram2d_clone(hist_.get()));
  CG_ASSERT(hist_w2_);
  CG_DEBUG("Hist2D:buildFromBins") << "Booking a 2D correlation plot with " << s("bin", x_bins.size(), true)
                                   << " in range x=" << x_bins << " and " << s("bin", y_bins.size(), true)
                                   << " in range y=" << y_bins << ".";
}

void Hist2D::buildFromRange(size_t num_bins_x, const Limits& xrange, size_t num_bins_y, const Limits& y_range) {
  if (xrange.range() <= 0. || y_range.range() <= 0.)
    throw CG_ERROR("Hist2D:buildFromRange") << "Invalid range for binning: " << xrange << "x" << y_range << ".";
  if (num_bins_x < 1 || num_bins_y < 1)
    throw CG_ERROR("Hist2D:buildFromRange") << "Building a 2D histogram requires at least 1x1 bin.";
  hist_.reset(gsl_histogram2d_alloc(num_bins_x, num_bins_y));
  CG_ASSERT(hist_);
  if (auto ret =
          gsl_histogram2d_set_ranges_uniform(hist_.get(), xrange.min(), xrange.max(), y_range.min(), y_range.max());
      ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:buildFromRange") << gsl_strerror(ret);
  hist_w2_ = gsl_histogram2d_ptr(gsl_histogram2d_clone(hist_.get()));
  CG_ASSERT(hist_w2_);
  CG_DEBUG("Hist2D:buildFromRange") << "Booking a 2D correlation plot with " << s("bin", num_bins_x, true)
                                    << " in range " << xrange << " and " << s("bin", num_bins_y, true) << " in range "
                                    << y_range << ".";
}

bool Hist2D::empty() const { return integral(true) == 0.; }

void Hist2D::clear() {
  CG_ASSERT(hist_);
  CG_ASSERT(hist_w2_);
  gsl_histogram2d_reset(hist_.get());
  gsl_histogram2d_reset(hist_w2_.get());
}

void Hist2D::fill(double x, double y, double weight) {
  CG_ASSERT(hist_);
  CG_ASSERT(hist_w2_);
  {  // reduce the scope of 'ret'
    auto ret = gsl_histogram2d_accumulate(hist_.get(), x, y, weight);
    if (ret == GSL_SUCCESS) {
      gsl_histogram2d_accumulate(hist_w2_.get(), x, y, weight * weight);
      return;
    }
    if (ret != GSL_EDOM)
      throw CG_ERROR("Hist2D:fill") << gsl_strerror(ret);
  }
  const auto &x_range = rangeX(), &y_range = rangeY();
  if (x_range.contains(x)) {
    if (y < y_range.min())
      out_of_range_values_[contents_t::IN_LT] += weight;
    else
      out_of_range_values_[contents_t::IN_GT] += weight;
  } else if (x < x_range.min()) {
    if (y_range.contains(y))
      out_of_range_values_[contents_t::LT_IN] += weight;
    else if (y < y_range.min())
      out_of_range_values_[contents_t::LT_LT] += weight;
    else
      out_of_range_values_[contents_t::LT_GT] += weight;
  } else {
    if (y_range.contains(y))
      out_of_range_values_[contents_t::GT_IN] += weight;
    else if (y < y_range.min())
      out_of_range_values_[contents_t::GT_LT] += weight;
    else
      out_of_range_values_[contents_t::GT_GT] += weight;
  }
}

void Hist2D::add(Hist2D oth, double scaling) {
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
  gsl_histogram2d_scale(oth.hist_w2_.get(), scl);
  if (auto ret = gsl_histogram2d_add(hist_.get(), oth.hist_.get()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:add") << gsl_strerror(ret);
  gsl_histogram2d_add(hist_w2_.get(), oth.hist_w2_.get());
  out_of_range_values_ += scaling * oth.out_of_range_values_;
}

void Hist2D::scale(double scaling) {
  CG_ASSERT(hist_);
  if (const auto ret = gsl_histogram2d_scale(hist_.get(), scaling); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:scale") << gsl_strerror(ret);
  gsl_histogram2d_scale(hist_w2_.get(), scaling * scaling);
}

size_t Hist2D::nbinsX() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_nx(hist_.get());
}

Limits Hist2D::rangeX() const {
  CG_ASSERT(hist_);
  return Limits{gsl_histogram2d_xmin(hist_.get()), gsl_histogram2d_xmax(hist_.get())};
}

Limits Hist2D::binRangeX(size_t bin) const {
  CG_ASSERT(hist_);
  Limits range;
  if (const auto ret = gsl_histogram2d_get_xrange(hist_.get(), bin, &range.min(), &range.max()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
  return range;
}

std::vector<double> Hist2D::binsX(BinMode mode) const {
  const auto bins = extractBins(mode, nbinsX(), std::bind(&Hist2D::binRangeX, this, std::placeholders::_1));
  return std::vector(bins.begin(), bins.end());
}

size_t Hist2D::nbinsY() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_ny(hist_.get());
}

Limits Hist2D::rangeY() const {
  CG_ASSERT(hist_);
  return Limits{gsl_histogram2d_ymin(hist_.get()), gsl_histogram2d_ymax(hist_.get())};
}

Limits Hist2D::binRangeY(size_t bin) const {
  CG_ASSERT(hist_);
  Limits range;
  if (const auto ret = gsl_histogram2d_get_yrange(hist_.get(), bin, &range.min(), &range.max()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
  return range;
}

std::vector<double> Hist2D::binsY(BinMode mode) const {
  const auto bins = extractBins(mode, nbinsY(), std::bind(&Hist2D::binRangeY, this, std::placeholders::_1));
  return std::vector(bins.begin(), bins.end());
}

std::pair<size_t, size_t> Hist2D::bin(double x, double y) const {
  std::pair<size_t, size_t> bin_ids;
  if (const auto ret = gsl_histogram2d_find(hist_.get(), x, y, &bin_ids.first, &bin_ids.second); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:bin") << "Failed to retrieve bin index for values (" << x << ", " << y
                                 << "): " << gsl_strerror(ret);
  return bin_ids;
}

Value Hist2D::value(size_t bin_x, size_t bin_y) const {
  CG_ASSERT(hist_);
  return Value{gsl_histogram2d_get(hist_.get(), bin_x, bin_y),
               std::sqrt(gsl_histogram2d_get(hist_w2_.get(), bin_x, bin_y))};
}

void Hist2D::setValue(size_t bin_x, size_t bin_y, Value val) {
  const auto bin_centre_x = binRangeX(bin_x).x(0.5), bin_centre_y = binRangeY(bin_y).x(0.5);
  const auto val_old = value(bin_x, bin_y);
  if (const auto ret = gsl_histogram2d_accumulate(hist_.get(), bin_centre_x, bin_centre_y, val - val_old);
      ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:setValue") << "Failed to accumulate values histogram. GSL error: " << gsl_strerror(ret);
  if (const auto ret = gsl_histogram2d_accumulate(
          hist_w2_.get(),
          bin_centre_x,
          bin_centre_y,
          val.uncertainty() * val.uncertainty() - val_old.uncertainty() * val_old.uncertainty());
      ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:setValue") << "Failed to accumulate square sum histogram. GSL error: " << gsl_strerror(ret);
}

double Hist2D::meanX() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_xmean(hist_.get());
}

double Hist2D::rmsX() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_xsigma(hist_.get());
}

double Hist2D::meanY() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_ymean(hist_.get());
}

double Hist2D::rmsY() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_ysigma(hist_.get());
}

double Hist2D::minimum() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_min_val(hist_.get());
}

double Hist2D::maximum() const {
  CG_ASSERT(hist_);
  return gsl_histogram2d_max_val(hist_.get());
}

double Hist2D::integral(bool include_out_of_range) const {
  CG_ASSERT(hist_);
  auto integral_value = gsl_histogram2d_sum(hist_.get());
  if (include_out_of_range)
    integral_value += out_of_range_values_.total();
  return integral_value;
}

size_t Hist2D::contents_t::total() const { return std::accumulate(begin(), end(), 0); }

Hist2D::contents_t& Hist2D::contents_t::operator+=(const contents_t& oth) {
  for (size_t i = 0; i < num_content; ++i)
    operator[](i) += oth.at(i);
  return *this;
}

std::pair<double, double> Hist2D::sample(RandomGenerator& rng) const {
  if (!pdf_) {
    pdf_.reset(gsl_histogram2d_pdf_alloc(nbinsX(), nbinsY()));
    if (const auto ret = gsl_histogram2d_pdf_init(pdf_.get(), hist_.get()); ret != GSL_SUCCESS)
      throw CG_FATAL("Hist2D:sample") << "Failed to allocate the histogram PDF. GSL yielded: " << gsl_strerror(ret);
  }
  std::pair<double, double> value;
  const auto xi = rng.uniform(), yi = rng.uniform();
  if (const auto ret = gsl_histogram2d_pdf_sample(pdf_.get(), xi, yi, &value.first, &value.second); ret != GSL_SUCCESS)
    throw CG_FATAL("Hist2D:sample") << "Failed to sample point (" << xi << ", " << yi
                                    << ") from the histogram PDF. GSL yielded: " << gsl_strerror(ret);
  return value;
}

namespace cepgen::utils {
  Hist2D::contents_t operator*(double scaling, const Hist2D::contents_t& oth) {
    Hist2D::contents_t tmp = oth;
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), [&scaling](const auto& bin) { return bin * scaling; });
    return tmp;
  }

  std::ostream& operator<<(std::ostream& os, const Hist2D::contents_t& cnt) {
    return os << format(
               "%10zu | %10zu | %10zu\n"
               "%10zu | %10s | %10zu\n"
               "%10zu | %10zu | %10zu",
               cnt.at(Hist2D::contents_t::LT_LT),
               cnt.at(Hist2D::contents_t::LT_IN),
               cnt.at(Hist2D::contents_t::LT_GT),
               cnt.at(Hist2D::contents_t::IN_LT),
               "-",
               cnt.at(Hist2D::contents_t::IN_GT),
               cnt.at(Hist2D::contents_t::GT_LT),
               cnt.at(Hist2D::contents_t::GT_IN),
               cnt.at(Hist2D::contents_t::GT_GT));
  }
}  // namespace cepgen::utils
