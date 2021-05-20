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
      Hist();
      explicit Hist(const Hist&);
      virtual ~Hist();

      virtual void scale(double) = 0;
      virtual double integral() const = 0;
      virtual void draw(std::ostream&, size_t width) const = 0;

      void setName(const std::string& name) { name_ = name; }
      void setXlabel(const std::string& lab) { xlabel_ = lab; }
      void setYlabel(const std::string& lab) { ylabel_ = lab; }
      void setLog(bool log = true) { log_ = log; }

    protected:
      std::string name_;
      std::string xlabel_, ylabel_;
      bool log_;
    };

    /// 1D histogram container
    class Hist1D : public Hist {
    public:
      Hist1D(size_t num_bins_x, const Limits&);
      Hist1D(const Hist1D&);

      void fill(double x, double weight = 1.);
      void add(Hist1D, double scaling = 1.);
      void scale(double) override;

      Limits xrange() const;

      double mean() const;
      double rms() const;
      double integral() const override;
      void draw(std::ostream&, size_t width = 50) const override;

    private:
      static constexpr char CHAR = '*';

      struct gsl_histogram_deleter {
        void operator()(gsl_histogram* h) { gsl_histogram_free(h); }
      };
      typedef std::unique_ptr<gsl_histogram, gsl_histogram_deleter> gsl_histogram_ptr;
      gsl_histogram_ptr hist_;
      size_t underflow_, overflow_;
    };
    /// 2D histogram container
    class Hist2D : public Hist {
    public:
      Hist2D(size_t num_bins_x, const Limits& xlim, size_t num_bins_y, const Limits& ylim);
      Hist2D(const Hist2D&);

      void fill(double x, double y, double weight = 1.);
      void add(Hist2D, double scaling = 1.);
      void scale(double) override;

      Limits xrange() const;
      Limits yrange() const;

      double meanX() const;
      double rmsX() const;
      double meanY() const;
      double rmsY() const;
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
      gsl_histogram2d_ptr hist_;
      struct values_t {
        size_t LT_GT = 0ull, IN_GT = 0ull, GT_GT = 0ull;
        size_t LT_IN = 0ull, /* INSIDE  */ GT_IN = 0ull;
        size_t LT_LT = 0ull, IN_LT = 0ull, GT_LT = 0ull;
      } values_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
