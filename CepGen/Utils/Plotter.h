#ifndef CepGen_Utils_Plotter_h
#define CepGen_Utils_Plotter_h

#include "CepGen/Utils/Limits.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <vector>
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
      explicit Hist(const Hist&);
      virtual ~Hist();

      virtual void clear() = 0;
      virtual void scale(double) = 0;
      virtual double integral() const = 0;
      virtual double minimum() const = 0;
      virtual double maximum() const = 0;

      void setName(const std::string& name) { name_ = name; }

    protected:
      std::string name_;
    };

    class Drawable {
    public:
      Drawable() : log_(false) {}
      Drawable(const Drawable&);
      virtual void draw(std::ostream&, size_t width) const = 0;

      void setXlabel(const std::string& lab) { xlabel_ = lab; }
      void setYlabel(const std::string& lab) { ylabel_ = lab; }
      void setLog(bool log = true) { log_ = log; }

    protected:
      std::string xlabel_, ylabel_;
      bool log_;
    };

    class Drawable1D : public Drawable {
    public:
      explicit Drawable1D(size_t width) : width_(width) {}

    protected:
      struct value_t {
        double coord, value, value_unc = 0.;
        std::string label = "";
      };
      typedef std::vector<value_t> values_t;

      void drawValues(std::ostream&, const values_t&) const;

      const size_t width_;

    private:
      static constexpr char CHAR = '*', ERR_CHAR = '-';
    };

    class Drawable2D : public Drawable {};

    /// 1D histogram container
    class Hist1D : public Hist, public Drawable1D {
    public:
      Hist1D(size_t num_bins_x, const Limits&);
      Hist1D(const std::vector<double>&);
      Hist1D(const Hist1D&);

      void clear() override;
      void fill(double x, double weight = 1.);
      void add(Hist1D, double scaling = 1.);
      void scale(double) override;

      double value(size_t bin) const;
      double valueUnc(size_t bin) const;

      size_t nbins() const;
      Limits range() const;
      Limits binRange(size_t bin) const;

      double mean() const;
      double rms() const;
      double minimum() const override;
      double maximum() const override;
      double integral() const override;
      void draw(std::ostream&, size_t width = 50) const override;

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
      Hist2D(size_t num_bins_x, const Limits& xlim, size_t num_bins_y, const Limits& ylim);
      Hist2D(const std::vector<double>&, const std::vector<double>&);
      Hist2D(const Hist2D&);

      void clear() override;
      void fill(double x, double y, double weight = 1.);
      void add(Hist2D, double scaling = 1.);
      void scale(double) override;

      double value(size_t bin_x, size_t bin_y) const;
      double valueUnc(size_t bin_x, size_t bin_y) const;

      size_t nbinsX() const;
      Limits rangeX() const;
      Limits binRangeX(size_t bin) const;
      size_t nbinsY() const;
      Limits rangeY() const;
      Limits binRangeY(size_t bin) const;

      double meanX() const;
      double rmsX() const;
      double meanY() const;
      double rmsY() const;
      double minimum() const override;
      double maximum() const override;
      double integral() const override;
      void draw(std::ostream&, size_t width = 50) const override;

    private:
      // greyscale ascii art from http://paulbourke.net/dataformats/asciiart/
      //static constexpr const char* CHARS = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
      //static constexpr const char* CHARS = " .:-=+*#%@";
      static constexpr const char* CHARS = " .:oO0@%#";

      struct gsl_histogram2d_deleter {
        void operator()(gsl_histogram2d* h) { gsl_histogram2d_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter> gsl_histogram2d_ptr;
      gsl_histogram2d_ptr hist_, hist_w2_;
      struct values_t {
        size_t LT_GT = 0ull, IN_GT = 0ull, GT_GT = 0ull;
        size_t LT_IN = 0ull, /* INSIDE  */ GT_IN = 0ull;
        size_t LT_LT = 0ull, IN_LT = 0ull, GT_LT = 0ull;
      } values_;
    };

    class Graph1D : public Drawable1D {
    public:
      Graph1D();

      void addPoint(double x, double y);

      void draw(std::ostream&, size_t width = 50) const override;

    private:
      values_t values_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
