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

#ifndef CepGen_Utils_Histogram_h
#define CepGen_Utils_Histogram_h

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <array>
#include <functional>
#include <memory>
#include <set>

#include "CepGen/Utils/Drawable.h"

namespace cepgen {
  class ParametersList;
  namespace utils {
    class RandomGenerator;
    /**
     * \brief Generic text-based plotting utility
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class Histogram {
    public:
      Histogram() = default;
      virtual ~Histogram() = default;

      enum BinMode { low = 0, high, both };

      /// Reset the histogram
      virtual void clear() = 0;
      /// Rescale all histogram bins by a constant factor
      virtual void scale(double) = 0;
      /// Compute the histogram integral
      virtual double integral(bool include_out_of_range = false) const = 0;
      /// Retrieve the maximum bin value
      virtual double minimum() const = 0;
      /// Retrieve the minimum bin value
      virtual double maximum() const = 0;
      /// Normalise the histogram to a given constant
      void normalise(double integ = 1.);

    protected:
      /// Extract the list of bin limits
      /// \param[in] mode type of extraction (low/high/low-high)
      /// \param[in] num_bins total number of bins
      /// \param[in] bins_extractor method used to extract range of one single bin
      std::set<double> extractBins(BinMode mode,
                                   size_t num_bins,
                                   const std::function<Limits(size_t)>& bins_extractor) const;
    };

    /// 1D histogram container
    class Hist1D : public Histogram, public Drawable {
    public:
      /// Build a histogram from user-steeed parameters
      explicit Hist1D(const ParametersList&);
      /// Build a histogram from uniform-width bins
      explicit Hist1D(size_t num_bins_x, const Limits&, const std::string& name = "", const std::string& title = "");
      /// Build a histogram from variable-width bins
      explicit Hist1D(const std::vector<double>& bins, const std::string& name = "", const std::string& title = "");
      Hist1D(const Hist1D&);  ///< Copy constructor

      void clear() override;
      /// Increment the histogram with one value
      void fill(double x, double weight = 1.);
      /// Bin-to-bin addition of another histogram to this one
      void add(Hist1D, double scaling = 1.);
      void scale(double) override;
      /// Sample individual "events" from a distribution
      double sample(RandomGenerator&) const;

      /// Perform a chi^2 test between two histograms
      /// \param[out] ndf number of degrees of freedom (non-empty bins)
      /// \return chi^2-value of the equivalence test
      double chi2test(const Hist1D&, size_t& ndf) const;

      /// Retrieve the value + uncertainty for all bins
      std::vector<Value> values() const;
      /// Retrieve the value + uncertainty for one bin
      Value value(size_t bin) const;
      /// Set the value + uncertainty for one bin
      void setValue(size_t bin, Value value);

      /// Axis content
      axis_t axis() const;
      /// Number of histogram bins
      size_t nbins() const;
      /// Axis range
      Limits range() const;
      /// Range for a single bin
      Limits binRange(size_t bin) const;
      /// List of bins limits (nbins + 1 values if min-max, nbins values otherwise)
      std::vector<double> bins(BinMode) const;
      /// Retrieve the bin index for a x value
      size_t bin(double x) const;

      /// Compute the mean histogram value over full range
      double mean() const;
      /// Compute the root-mean-square value over full range
      double rms() const;
      double minimum() const override;
      double maximum() const override;
      double integral(bool = false) const override;
      size_t underflow() const { return underflow_; }
      size_t overflow() const { return overflow_; }

      bool isHist1D() const override final { return true; }

    private:
      void buildFromBins(const std::vector<double>&);
      void buildFromRange(size_t, const Limits&);

      struct gsl_histogram_deleter {
        void operator()(gsl_histogram* h) { gsl_histogram_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram, gsl_histogram_deleter> gsl_histogram_ptr;
      gsl_histogram_ptr hist_, hist_w2_;
      size_t underflow_{0ull}, overflow_{0ull};
      struct gsl_histogram_pdf_deleter {
        void operator()(gsl_histogram_pdf* h) { gsl_histogram_pdf_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram_pdf, gsl_histogram_pdf_deleter> gsl_histogram_pdf_ptr;
      mutable gsl_histogram_pdf_ptr pdf_;
    };

    /// 2D histogram container
    class Hist2D : public Histogram, public Drawable {
    public:
      /// Build a histogram from user-steeed parameters
      explicit Hist2D(const ParametersList&);
      /// Build a histogram from uniform-width bins
      explicit Hist2D(size_t num_bins_x,
                      const Limits& xlim,
                      size_t num_bins_y,
                      const Limits& ylim,
                      const std::string& name = "",
                      const std::string& title = "");
      /// Build a histogram from variable-width bins
      explicit Hist2D(const std::vector<double>& xbins,
                      const std::vector<double>& ybins,
                      const std::string& name = "",
                      const std::string& title = "");
      Hist2D(const Hist2D&);  ///< Copy constructor

      void clear() override;
      /// Fill the histogram with one value
      void fill(double x, double y, double weight = 1.);
      /// Fill the histogram with one value
      inline void fill(const std::pair<double, double>& xy, double weight = 1.) { fill(xy.first, xy.second, weight); }
      /// Bin-by-bin addition of another histogram to this one
      void add(Hist2D, double scaling = 1.);
      void scale(double) override;
      /// Sample individual "events" from a distribution
      std::pair<double, double> sample(RandomGenerator&) const;

      /// Retrieve the value + uncertainty for one bin
      Value value(size_t bin_x, size_t bin_y) const;
      /// Set the value + uncertainty for one bin
      void setValue(size_t bin_x, size_t bin_y, Value value);

      /// Number of x-axis bins
      size_t nbinsX() const;
      /// x-axis range
      Limits rangeX() const;
      /// Range for a single x-axis bin
      Limits binRangeX(size_t bin) const;
      /// List of x-bins limits (nbinsX + 1 values if min-max, nbins values otherwise)
      std::vector<double> binsX(BinMode) const;
      /// Number of y-axis bins
      size_t nbinsY() const;
      /// y-axis range
      Limits rangeY() const;
      /// Range for a single y-axis bin
      Limits binRangeY(size_t bin) const;
      /// List of y-bins limits (nbinsY + 1 values if min-max, nbins values otherwise)
      std::vector<double> binsY(BinMode) const;
      /// Retrieve the bin indices for a (x, y) value
      std::pair<size_t, size_t> bin(double x, double y) const;

      /// Compute the mean histogram value over full x-axis range
      double meanX() const;
      /// Compute the root-mean-square value over full x-axis range
      double rmsX() const;
      /// Compute the mean histogram value over full y-axis range
      double meanY() const;
      /// Compute the root-mean-square value over full y-axis range
      double rmsY() const;
      double minimum() const override;
      double maximum() const override;
      double integral(bool = false) const override;

      struct contents_t : public std::array<size_t, 8> {
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

      bool isHist2D() const override final { return true; }

    private:
      void buildFromBins(const std::vector<double>&, const std::vector<double>&);
      void buildFromRange(size_t, const Limits&, size_t, const Limits&);

      struct gsl_histogram2d_deleter {
        void operator()(gsl_histogram2d* h) { gsl_histogram2d_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter> gsl_histogram2d_ptr;
      gsl_histogram2d_ptr hist_, hist_w2_;
      contents_t out_of_range_values_;
      struct gsl_histogram2d_pdf_deleter {
        void operator()(gsl_histogram2d_pdf* h) { gsl_histogram2d_pdf_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram2d_pdf, gsl_histogram2d_pdf_deleter> gsl_histogram2d_pdf_ptr;
      mutable gsl_histogram2d_pdf_ptr pdf_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
