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

#ifndef CepGen_Utils_Histogram_h
#define CepGen_Utils_Histogram_h

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <array>
#include <functional>
#include <memory>
#include <set>

#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Utils/Drawable.h"

namespace cepgen::utils {
  class RandomGenerator;
  /// Generic container for binned distributions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class Histogram {
  public:
    Histogram() = default;
    virtual ~Histogram() = default;

    enum BinMode { low = 0, high, both };

    virtual void clear() = 0;                    ///< Reset the histogram
    virtual void scale(double) = 0;              ///< Rescale all histogram bins by a constant factor
    void normalise(double integral_value = 1.);  ///< Normalise the histogram to a given constant
    virtual double integral(bool include_out_of_range = false) const = 0;  ///< Compute the histogram integral

    virtual double minimum() const = 0;  ///< Retrieve the maximum bin value
    virtual double maximum() const = 0;  ///< Retrieve the minimum bin value

  protected:
    /// Extract the list of bin limits
    /// \param[in] mode type of extraction (low/high/low-high)
    /// \param[in] num_bins total number of bins
    /// \param[in] bins_extractor method used to extract range of one single bin
    static std::set<double> extractBins(BinMode mode,
                                        size_t num_bins,
                                        const std::function<Limits(size_t)>& bins_extractor);
  };

  /// 1D histogram container
  class Hist1D : public Histogram, public Drawable {
  public:
    explicit Hist1D(const ParametersList&);  ///< Build a histogram from user-steered parameters
    /// Build a histogram from uniform-width bins
    explicit Hist1D(size_t num_bins_x, const Limits&, const std::string& name = "", const std::string& title = "");
    /// Build a histogram from variable-width bins
    explicit Hist1D(const std::vector<double>& bins, const std::string& name = "", const std::string& title = "");
    Hist1D(const Hist1D&);  ///< Copy constructor

    static ParametersDescription description();

    bool empty() const override;
    void clear() override;
    void fill(double x, double weight = 1.);  ///< Increment the histogram with one value
    void add(Hist1D, double scaling = 1.);    ///< Bin-to-bin addition of another histogram to this one
    void scale(double) override;
    double sample(RandomGenerator&) const;  ///< Sample individual "events" from a distribution

    /// Perform a chi^2 test between two histograms
    /// \param[out] ndf number of degrees of freedom (non-empty bins)
    /// \return chi^2-value of the equivalence test
    double chi2test(const Hist1D&, size_t& ndf) const;

    std::vector<Value> values() const;       ///< Retrieve the value + uncertainty for all bins
    Value value(size_t bin) const;           ///< Retrieve the value + uncertainty for one bin
    void setValue(size_t bin, Value value);  ///< Set the value + uncertainty for one bin

    axis_t axis() const;                      ///< Axis content
    size_t nbins() const;                     ///< Number of histogram bins
    Limits range() const;                     ///< Axis range
    Limits binRange(size_t bin) const;        ///< Range for a single bin
    std::vector<double> bins(BinMode) const;  ///< List of bins limits (nbins+1 values if min-max, nbins otherwise)
    size_t bin(double x) const;               ///< Retrieve the bin index for an x value

    double mean() const;  ///< Compute the mean histogram value over full range
    double rms() const;   ///< Compute the root-mean-square value over full range
    double minimum() const override;
    double maximum() const override;
    double integral(bool = false) const override;
    inline size_t underflow() const { return underflow_; }
    inline size_t overflow() const { return overflow_; }

    bool isHist1D() const final { return true; }

  private:
    void buildFromBins(const std::vector<double>&);
    void buildFromRange(size_t, const Limits&);

    struct gsl_histogram_deleter {
      void operator()(gsl_histogram* h) const { gsl_histogram_free(h); }
    };
    using gsl_histogram_ptr = std::unique_ptr<gsl_histogram, gsl_histogram_deleter>;
    gsl_histogram_ptr hist_, hist_w2_;
    size_t underflow_{0ull}, overflow_{0ull};
    struct gsl_histogram_pdf_deleter {
      void operator()(gsl_histogram_pdf* h) const { gsl_histogram_pdf_free(h); }
    };
    using gsl_histogram_pdf_ptr = std::unique_ptr<gsl_histogram_pdf, gsl_histogram_pdf_deleter>;
    mutable gsl_histogram_pdf_ptr pdf_;
  };

  /// 2D histogram container
  class Hist2D : public Histogram, public Drawable {
  public:
    explicit Hist2D(const ParametersList&);  ///< Build a histogram from user-steered parameters
    /// Build a histogram from uniform-width bins
    explicit Hist2D(size_t num_bins_x,
                    const Limits& x_limits,
                    size_t num_bins_y,
                    const Limits& y_limits,
                    const std::string& name = "",
                    const std::string& title = "");
    /// Build a histogram from variable-width bins
    explicit Hist2D(const std::vector<double>& x_bins,
                    const std::vector<double>& y_bins,
                    const std::string& name = "",
                    const std::string& title = "");
    Hist2D(const Hist2D&);  ///< Copy constructor

    static ParametersDescription description();

    bool empty() const override;
    void clear() override;
    void fill(double x, double y, double weight = 1.);  ///< Fill the histogram with one value
    /// Fill the histogram with one value
    inline void fill(const std::pair<double, double>& xy, double weight = 1.) { fill(xy.first, xy.second, weight); }
    void add(Hist2D, double scaling = 1.);  ///< Bin-by-bin addition of another histogram to this one
    void scale(double) override;
    std::pair<double, double> sample(RandomGenerator&) const;  ///< Sample individual "events" from a distribution

    Value value(size_t bin_x, size_t bin_y) const;           ///< Retrieve the value + uncertainty for one bin
    void setValue(size_t bin_x, size_t bin_y, Value value);  ///< Set the value + uncertainty for one bin

    size_t nbinsX() const;                     ///< Number of x-axis bins
    Limits rangeX() const;                     ///< x-axis range
    Limits binRangeX(size_t bin) const;        ///< Range for a single x-axis bin
    std::vector<double> binsX(BinMode) const;  ///< List of x-bins limits (nbinsX+1 values if min-max, nbins otherwise)
    size_t nbinsY() const;                     ///< Number of y-axis bins
    Limits rangeY() const;                     ///< y-axis range
    Limits binRangeY(size_t bin) const;        ///< Range for a single y-axis bin
    std::vector<double> binsY(BinMode) const;  ///< List of y-bins limits (nbinsY+1 values if min-max, nbins otherwise)
    std::pair<size_t, size_t> bin(double x, double y) const;  ///< Retrieve the bin indices for a (x, y) value

    double meanX() const;  ///< Compute the mean histogram value over full x-axis range
    double rmsX() const;   ///< Compute the root-mean-square value over full x-axis range
    double meanY() const;  ///< Compute the mean histogram value over full y-axis range
    double rmsY() const;   ///< Compute the root-mean-square value over full y-axis range
    double minimum() const override;
    double maximum() const override;
    double integral(bool = false) const override;

    struct contents_t : std::array<size_t, 8> {
      contents_t() { std::fill(begin(), end(), 0ull); }
      size_t total() const;
      friend std::ostream& operator<<(std::ostream&, const contents_t&);
      friend contents_t operator*(double, const contents_t&);
      contents_t& operator+=(const contents_t&);
      enum {
        LT_GT,
        IN_GT,
        GT_GT,
        LT_IN,
        /* INSIDE */ GT_IN,
        LT_LT,
        IN_LT,
        GT_LT,
        num_content
      };
    };
    const contents_t& outOfRange() const { return out_of_range_values_; }

    bool isHist2D() const final { return true; }

  private:
    void buildFromBins(const std::vector<double>&, const std::vector<double>&);
    void buildFromRange(size_t, const Limits&, size_t, const Limits&);

    struct gsl_histogram2d_deleter {
      void operator()(gsl_histogram2d* h) const { gsl_histogram2d_free(h); }
    };
    using gsl_histogram2d_ptr = std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter>;
    gsl_histogram2d_ptr hist_, hist_w2_;
    contents_t out_of_range_values_;
    struct gsl_histogram2d_pdf_deleter {
      void operator()(gsl_histogram2d_pdf* h) const { gsl_histogram2d_pdf_free(h); }
    };
    using gsl_histogram2d_pdf_ptr = std::unique_ptr<gsl_histogram2d_pdf, gsl_histogram2d_pdf_deleter>;
    mutable gsl_histogram2d_pdf_ptr pdf_;
  };
}  // namespace cepgen::utils

#endif
