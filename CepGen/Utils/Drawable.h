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

#ifndef CepGen_Utils_Drawable_h
#define CepGen_Utils_Drawable_h

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace utils {
    /// A generic object which can be drawn in the standard output
    class Drawable {
    public:
      explicit Drawable(const std::string& name = "", const std::string& title = "");
      Drawable(const Drawable&);  ///< Copy constructor

      /// Drawable name
      const std::string& name() const { return name_; }
      /// Set the drawable name
      void setName(const std::string& name) { name_ = name; }

      /// Drawable name
      const std::string& title() const { return name_; }
      /// Set the drawable title
      void setTitle(const std::string& title) { title_ = title; }

      /// x-axis label
      const std::string& xLabel() const { return xlabel_; }
      /// Set the x-axis label
      void setXlabel(const std::string& lab) { xlabel_ = lab; }
      /// y-axis label
      const std::string& yLabel() const { return ylabel_; }
      /// Set the y-axis label
      void setYlabel(const std::string& lab) { ylabel_ = lab; }

      /// Generic bin coordinate and its human-readable label
      struct coord_t {
        /// Sorting helper for axis coordinates
        bool operator<(const coord_t& oth) const { return value < oth.value; }
        double value;            ///< Bin central value
        std::string label = "";  ///< Human-readable description of the bin
      };
      /// Helper view of a pair of bin value and its uncertainty
      struct value_t {
        /// Sorting helper for bin values
        bool operator<(const value_t& oth) const { return value < oth.value; }
        double value;           ///< Single bin content
        double value_unc = 0.;  ///< Uncertainty on bin content
      };
      /// Metadata for an axis (coordinates and bins value)
      typedef std::map<coord_t, value_t> axis_t;
      /// Metadata for a two-dimensional axis definition (coordinates and bins values)
      typedef std::map<coord_t, axis_t> dualaxis_t;

      virtual bool isHist1D() const { return false; }   ///< Is this drawable a one-dimensional histogram?
      virtual bool isHist2D() const { return false; }   ///< Is this drawable a two-dimensional histogram?
      virtual bool isGraph1D() const { return false; }  ///< Is this drawable a one-dimensional graph?
      virtual bool isGraph2D() const { return false; }  ///< Is this drawable a two-dimensional graph?

    protected:
      std::string name_;    ///< Computer-readable name
      std::string title_;   ///< Human-readable title
      std::string xlabel_;  ///< x-axis title
      std::string ylabel_;  ///< y-axis title
    };

    /// 1D histogram container
    class Hist1D : public Histogram, public Drawable {
    public:
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

      /// Retrieve the value for one bin
      double value(size_t bin) const;
      /// Retrieve the absolute uncertainty on one bin value
      double valueUnc(size_t bin) const;

      /// Number of histogram bins
      size_t nbins() const;
      /// Axis range
      Limits range() const;
      /// Range for a single bin
      Limits binRange(size_t bin) const;

      /// Compute the mean histogram value over full range
      double mean() const;
      /// Compute the root-mean-square value over full range
      double rms() const;
      double minimum() const override;
      double maximum() const override;
      double integral() const override;
      size_t underflow() const { return underflow_; }
      size_t overflow() const { return overflow_; }

      bool isHist1D() const override { return true; }

    private:
      struct gsl_histogram_deleter {
        void operator()(gsl_histogram* h) { gsl_histogram_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram, gsl_histogram_deleter> gsl_histogram_ptr;
      gsl_histogram_ptr hist_, hist_w2_;
      size_t underflow_{0ull}, overflow_{0ull};
    };

    /// 2D histogram container
    class Hist2D : public Histogram, public Drawable {
    public:
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
      /// Bin-by-bin addition of another histogram to this one
      void add(Hist2D, double scaling = 1.);
      void scale(double) override;

      /// Retrieve the value for one bin
      double value(size_t bin_x, size_t bin_y) const;
      /// Retrieve the absolute uncertainty on one bin value
      double valueUnc(size_t bin_x, size_t bin_y) const;

      /// Number of x-axis bins
      size_t nbinsX() const;
      /// x-axis range
      Limits rangeX() const;
      /// Range for a single x-axis bin
      Limits binRangeX(size_t bin) const;
      /// Number of y-axis bins
      size_t nbinsY() const;
      /// y-axis range
      Limits rangeY() const;
      /// Range for a single y-axis bin
      Limits binRangeY(size_t bin) const;

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
      double integral() const override;

      struct contents_t {
        inline size_t total() const { return LT_GT + IN_GT + GT_GT + LT_IN + GT_IN + LT_LT + IN_LT + GT_LT; }
        size_t LT_GT{0ull}, IN_GT{0ull}, GT_GT{0ull};
        size_t LT_IN{0ull}, /* INSIDE */ GT_IN{0ull};
        size_t LT_LT{0ull}, IN_LT{0ull}, GT_LT{0ull};
      };
      const contents_t& content() const { return values_; }

      bool isHist2D() const override { return true; }

    private:
      struct gsl_histogram2d_deleter {
        void operator()(gsl_histogram2d* h) { gsl_histogram2d_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter> gsl_histogram2d_ptr;
      gsl_histogram2d_ptr hist_, hist_w2_;
      contents_t values_;
    };

    /// A one-dimensional graph object
    class Graph1D : public Drawable {
    public:
      explicit Graph1D(const std::string& name = "", const std::string& title = "");

      /// Add one value to the graph
      void addPoint(double x, double y);
      /// Retrieve all values in the graph
      const axis_t& points() const { return values_; }

      bool isGraph1D() const override { return true; }

    private:
      axis_t values_;
    };

    /// A two-dimensional graph object
    class Graph2D : public Drawable {
    public:
      explicit Graph2D(const std::string& name = "", const std::string& title = "");

      /// Add one value to the graph
      void addPoint(double x, double y, double z);
      /// Retrieve all values in the graph
      const dualaxis_t& points() const { return values_; }
      /// List all values registered in the graph
      void dumpPoints(std::ostream&) const;

      bool isGraph2D() const override { return true; }

    private:
      dualaxis_t values_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
