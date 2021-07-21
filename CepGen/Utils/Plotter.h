#ifndef CepGen_Utils_Plotter_h
#define CepGen_Utils_Plotter_h

#include "CepGen/Utils/Limits.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <vector>
#include <map>
#include <string>
#include <memory>

namespace cepgen {
  namespace utils {
    /**
     * \brief Generic text-based plotting utility
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class Hist {
    public:
      Hist() = default;
      explicit Hist(const Hist&);  ///< Copy constructor
      virtual ~Hist();

      /// Reset the histogram
      virtual void clear() = 0;
      /// Rescale all histogram bins by a constant factor
      virtual void scale(double) = 0;
      /// Compute the histogram integral
      virtual double integral() const = 0;
      /// Retrieve the maximum bin value
      virtual double minimum() const = 0;
      /// Retrieve the minimum bin value
      virtual double maximum() const = 0;

      /// Set the histogram name
      void setName(const std::string& name) { name_ = name; }

    protected:
      std::string name_;  ///< Histogram human-readable name
    };

    /// A generic object which can be drawn in the standard output
    class Drawable {
    public:
      Drawable() : width_(50ul), log_(false) {}
      Drawable(const Drawable&);  ///< Copy constructor

      /// Main drawing method and its standard output
      /// \param[out] os Output stream to use for drawing
      virtual void draw(std::ostream& os) const = 0;

      /// Set the output width
      void setWidth(size_t width) { width_ = width; }
      /// Set the x-axis label
      void setXlabel(const std::string& lab) { xlabel_ = lab; }
      /// Set the y-axis label
      void setYlabel(const std::string& lab) { ylabel_ = lab; }
      /// Switch logarithmic view
      void setLog(bool log = true) { log_ = log; }

    protected:
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
      /// Sorting helper for the axis metadata container
      struct map_elements {
        /// Sort two values in an axis
        bool operator()(const std::pair<coord_t, value_t>& lhs, const std::pair<coord_t, value_t>& rhs) {
          return lhs.second.value < rhs.second.value;
        }
      };
      size_t width_;        ///< Plot width, in TTY characters
      std::string xlabel_;  ///< x-axis title
      std::string ylabel_;  ///< y-axis title
      bool log_;            ///< Switch on/off the logarithmic z-axis
    };

    /// Any drawable with one axis
    class Drawable1D : public Drawable {
    protected:
      /// Draw all values for one axis
      void drawValues(std::ostream&, const axis_t&) const;

    private:
      static constexpr char CHAR = '*', ERR_CHAR = '-';
    };

    /// Any drawable with two axes
    class Drawable2D : public Drawable {
    protected:
      /// Metadata for a two-dimensional axis definition (coordinates and bins values)
      typedef std::map<coord_t, axis_t> dualaxis_t;

      /// Draw all values for two axes
      void drawValues(std::ostream&, const dualaxis_t&) const;

    private:
      // greyscale ascii art from http://paulbourke.net/dataformats/asciiart/
      //static constexpr const char* CHARS = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
      //static constexpr const char* CHARS = " .:-=+*#%@";
      static constexpr const char* CHARS = " .:oO0@%#";
      static const int kColours[];
      static constexpr const char NEG_CHAR = '-';
    };

    /// 1D histogram container
    class Hist1D : public Hist, public Drawable1D {
    public:
      /// Build a histogram from uniform-width bins
      Hist1D(size_t num_bins_x, const Limits&);
      /// Build a histogram from variable-width bins
      explicit Hist1D(const std::vector<double>& bins);
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
      void draw(std::ostream&) const override;

    private:
      struct gsl_histogram_deleter {
        void operator()(gsl_histogram* h) { gsl_histogram_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram, gsl_histogram_deleter> gsl_histogram_ptr;
      gsl_histogram_ptr hist_, hist_w2_;
      size_t underflow_, overflow_;
    };

    /// 2D histogram container
    class Hist2D : public Hist, public Drawable2D {
    public:
      /// Build a histogram from uniform-width bins
      Hist2D(size_t num_bins_x, const Limits& xlim, size_t num_bins_y, const Limits& ylim);
      /// Build a histogram from variable-width bins
      Hist2D(const std::vector<double>& xbins, const std::vector<double>& ybins);
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
      void draw(std::ostream&) const override;

    private:
      struct gsl_histogram2d_deleter {
        void operator()(gsl_histogram2d* h) { gsl_histogram2d_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter> gsl_histogram2d_ptr;
      gsl_histogram2d_ptr hist_, hist_w2_;
      struct contents_t {
        inline size_t total() const { return LT_GT + IN_GT + GT_GT + LT_IN + GT_IN + LT_LT + IN_LT + GT_LT; }
        std::string summary() const;
        size_t LT_GT = 0ull, IN_GT = 0ull, GT_GT = 0ull;
        size_t LT_IN = 0ull, /* INSIDE  */ GT_IN = 0ull;
        size_t LT_LT = 0ull, IN_LT = 0ull, GT_LT = 0ull;
      } values_;
    };

    /// A one-dimensional graph object
    class Graph1D : public Drawable1D {
    public:
      Graph1D() = default;

      /// Add one value to the graph
      void addPoint(double x, double y);

      void draw(std::ostream&) const override;

    private:
      axis_t values_;
    };

    /// A two-dimensional graph object
    class Graph2D : public Drawable2D {
    public:
      Graph2D() = default;

      /// Add one value to the graph
      void addPoint(double x, double y, double z);
      /// List all values registered in the graph
      void dumpPoints(std::ostream&) const;

      void draw(std::ostream&) const override;

    private:
      dualaxis_t values_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
