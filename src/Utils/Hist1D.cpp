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
using namespace std::string_literals;

Hist1D::Hist1D(const ParametersList& params)
    : Drawable(params.get<std::string>("name"), params.get<std::string>("title")) {
  if (const auto& x_bins = params.get<std::vector<double> >("xbins"s); x_bins.size() > 1)
    buildFromBins(x_bins);
  else if (const auto& xrange = params.get<Limits>("xrange"); xrange.valid())
    buildFromRange(params.get<int>("nbins"s) > 0 ? params.get<int>("nbins"s) : params.get<int>("nbinsX"s), xrange);
  else
    throw CG_FATAL("Hist1D") << "Failed to build a 1D histogram with user parameters: " << params << ".";
}

Hist1D::Hist1D(size_t num_bins_x, const Limits& xrange, const std::string& name, const std::string& title)
    : Drawable(name, title) {
  if (num_bins_x == 0)
    throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
  buildFromRange(num_bins_x, xrange);
}

Hist1D::Hist1D(const std::vector<double>& x_bins, const std::string& name, const std::string& title)
    : Drawable(name, title) {
  if (x_bins.empty())
    throw CG_ERROR("Hist1D") << "Number of bins must be strictly positive!";
  buildFromBins(x_bins);
}

Hist1D::Hist1D(const Hist1D& oth)
    : Histogram(oth),
      Drawable(oth),
      hist_(gsl_histogram_clone(oth.hist_.get())),
      hist_w2_(gsl_histogram_clone(oth.hist_w2_.get())),
      underflow_(oth.underflow_),
      overflow_(oth.overflow_) {}

ParametersDescription Hist1D::description() {
  auto desc = ParametersDescription();
  desc.add("name"s, ""s).setDescription("Histogram filename-safe name");
  desc.add("title"s, ""s).setDescription("Histogram title/description");
  desc.add("log"s, false).setDescription("Plot logarithmic axis?");
  // x-axis attributes
  desc.add("xbins"s, std::vector<double>{}).setDescription("x-axis bins definition");
  desc.add("nbinsX"s, 25).setDescription("Bins multiplicity for x-axis");
  desc.add("xrange"s, Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
  return desc;
}

void Hist1D::buildFromBins(const std::vector<double>& bins) {
  if (bins.size() < 1)
    throw CG_ERROR("Hist1D:buildFromBins") << "Building a 1D histogram requires at least 1 bin.";
  hist_.reset(gsl_histogram_alloc(bins.size() - 1));
  CG_ASSERT(hist_);
  if (auto ret = gsl_histogram_set_ranges(hist_.get(), bins.data(), bins.size()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:buildFromBins") << gsl_strerror(ret);
  hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
  CG_ASSERT(hist_w2_);
  CG_DEBUG("Hist1D:buildFromBins") << "Booking a 1D histogram with " << s("bin", bins.size(), true) << " in range "
                                   << bins << ".";
}

void Hist1D::buildFromRange(size_t num_bins, const Limits& range) {
  if (range.range() <= 0.)
    throw CG_ERROR("Hist1D:buildFromRange") << "Invalid range for binning: " << range << ".";
  if (num_bins < 1)
    throw CG_ERROR("Hist1D:buildFromRange") << "Building a 1D histogram requires at least 1 bin.";
  hist_.reset(gsl_histogram_alloc(num_bins));
  CG_ASSERT(hist_);
  if (auto ret = gsl_histogram_set_ranges_uniform(hist_.get(), range.min(), range.max()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:buildFromRange") << gsl_strerror(ret);
  hist_w2_ = gsl_histogram_ptr(gsl_histogram_clone(hist_.get()));
  CG_ASSERT(hist_w2_);
  CG_DEBUG("Hist1D:buildFromRange") << "Booking a 1D histogram with " << s("bin", num_bins, true) << " in range "
                                    << range << ".";
}

bool Hist1D::empty() const { return integral(true) == 0.; }

void Hist1D::clear() {
  CG_ASSERT(hist_);
  CG_ASSERT(hist_w2_);
  gsl_histogram_reset(hist_.get());
  gsl_histogram_reset(hist_w2_.get());
}

void Hist1D::fill(double x, double weight) {
  CG_ASSERT(hist_);
  CG_ASSERT(hist_w2_);
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
  if (const auto ret = gsl_histogram_get_range(hist_.get(), bin, &range.min(), &range.max()); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:binRange") << "Bin " << bin << ": " << gsl_strerror(ret);
  return range;
}

std::vector<double> Hist1D::bins(BinMode mode) const {
  const auto bins = extractBins(mode, nbins(), std::bind(&Hist1D::binRange, this, std::placeholders::_1));
  return std::vector(bins.begin(), bins.end());
}

size_t Hist1D::bin(double x) const {
  size_t bin_id;
  if (const auto ret = gsl_histogram_find(hist_.get(), x, &bin_id); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist1D:bin") << "Failed to retrieve bin index for value " << x << ": " << gsl_strerror(ret);
  return bin_id;
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

void Hist1D::setValue(size_t bin, Value val) {
  const auto bin_centre = binRange(bin).x(0.5);
  const auto val_old = value(bin);
  if (const auto ret = gsl_histogram_accumulate(hist_.get(), bin_centre, val - val_old); ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:setValue") << "Failed to accumulate values histogram. GSL error: " << gsl_strerror(ret);
  if (const auto ret = gsl_histogram_accumulate(
          hist_w2_.get(),
          bin_centre,
          val.uncertainty() * val.uncertainty() - val_old.uncertainty() * val_old.uncertainty());
      ret != GSL_SUCCESS)
    throw CG_ERROR("Hist2D:setValue") << "Failed to accumulate square sum histogram. GSL error: " << gsl_strerror(ret);
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
  auto integral_value = gsl_histogram_sum(hist_.get());
  if (include_out_of_range)
    integral_value += underflow_ + overflow_;
  return integral_value;
}

double Hist1D::sample(RandomGenerator& rng) const {
  if (!pdf_) {
    pdf_.reset(gsl_histogram_pdf_alloc(nbins()));
    if (const auto ret = gsl_histogram_pdf_init(pdf_.get(), hist_.get()); ret != GSL_SUCCESS)
      throw CG_ERROR("Hist1D:sample") << "Failed to allocate the histogram PDF. GSL yielded: " << gsl_strerror(ret);
  }
  return gsl_histogram_pdf_sample(pdf_.get(), rng.uniform());
}

double Hist1D::chi2test(const Hist1D& oth, size_t& ndf_value) const {
  if (nbins() != oth.nbins())
    return 0.;
  double sum1{0.}, sum2{0.};
  for (size_t i = 0; i < nbins(); ++i) {
    const auto ru1 = value(i).relativeUncertainty(), ru2 = oth.value(i).relativeUncertainty();
    sum1 += ru1 > 0. ? 1. / ru1 / ru1 : 0.;
    sum2 += ru2 > 0. ? 1. / ru2 / ru2 : 0.;
  }
  if (sum1 == 0. || sum2 == 0.) {
    ndf_value = 0;
    return 0.;
  }
  // perform the test
  double chi2val = 0.;
  ndf_value = nbins();
  for (size_t i = 0; i < nbins(); ++i) {
    const auto ru1 = value(i).relativeUncertainty(), ru2 = oth.value(i).relativeUncertainty();
    const auto cnt1 = ru1 > 0. ? 1. / ru1 / ru1 : 0., cnt2 = ru2 > 0. ? 1. / ru2 / ru2 : 0.;
    if (cnt1 == 0. && cnt2 == 0.) {
      --ndf_value;
      continue;
    }
    chi2val += std::pow(sum2 * cnt1 - sum1 * cnt2, 2) / (cnt1 + cnt2);
  }
  chi2val /= sum1 * sum2;
  return chi2val;
}
