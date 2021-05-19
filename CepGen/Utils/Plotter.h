#ifndef CepGen_Utils_Plotter_h
#define CepGen_Utils_Plotter_h

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
    class Plotter {
    public:
      explicit Plotter();
      ~Plotter();

    private:
      struct info_t {
        std::vector<std::string> name;
        std::string xlabel, ylabel;
        bool log;
      };
      /// 1D histogram container
      class Hist1D {
      public:
        Hist1D(size_t num_bins_x, double min_x, double max_x);

        void setXlabel(const std::string& lab) { info_.xlabel = lab; }
        void fill(double x, double weight = 1.);
        void draw(std::ostream&, size_t width = 50) const;

      private:
        static constexpr char CHAR = '*';

        info_t info_;
        struct gsl_histogram_deleter {
          void operator()(gsl_histogram* h) { gsl_histogram_free(h); }
        };
        typedef std::unique_ptr<gsl_histogram, gsl_histogram_deleter> gsl_histogram_ptr;
        gsl_histogram_ptr hist_;
      };
      /// 2D histogram container
      class Hist2D {
      public:
        Hist2D(size_t num_bins_x, double min_x, double max_x, size_t num_bins_y, double min_y, double max_y);

        void setXlabel(const std::string& lab) { info_.xlabel = lab; }
        void setYlabel(const std::string& lab) { info_.ylabel = lab; }
        void fill(double x, double y, double weight = 1.);
        void draw(std::ostream&, size_t width = 50) const;

      private:
        // greyscale ascii art from http://paulbourke.net/dataformats/asciiart/
        //static constexpr const char* CHARS = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
        //static constexpr const char* CHARS = " .:-=+*#%@";
        static constexpr const char* CHARS = " .:oO0@%#";

        info_t info_;
        struct gsl_histogram2d_deleter {
          void operator()(gsl_histogram2d* h) { gsl_histogram2d_free(h); }
        };
        typedef std::unique_ptr<gsl_histogram2d, gsl_histogram2d_deleter> gsl_histogram2d_ptr;
        gsl_histogram2d_ptr hist_;
      };

      static Hist1D newHist1D(size_t num_bins, double min_x, double max_x);
    };
  }  // namespace utils
}  // namespace cepgen

#endif
