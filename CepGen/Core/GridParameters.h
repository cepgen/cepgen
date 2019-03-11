#ifndef CepGen_Core_GridParameters_h
#define CepGen_Core_GridParameters_h

#include <vector>
#include <map>
#include <gsl/gsl_rng.h>

namespace cepgen
{
  /// A parameters placeholder for the grid integration helper
  class GridParameters
  {
    public:
      typedef std::vector<unsigned short> coord_t;

      GridParameters( unsigned short ndim );

      size_t size() const { return max_; }
      const coord_t& n( size_t coord ) const;
      double globalMax() const { return f_max_global_; }
      double maxValue( size_t coord ) const;
      void setValue( size_t coord, double val );

      void shoot( const gsl_rng* rng, size_t coord, std::vector<double>& out ) const;

      /// Maximal number of dimensions handled by this integrator instance
      static constexpr unsigned short MAX_DIM = 15;
      /// Integration grid size parameter
      static constexpr unsigned short M_BIN = 3;
      static constexpr double INV_M_BIN = 1./M_BIN;

      /// Has the generation been prepared?
      bool gen_prepared;
      double correc;
      double correc2;
      double f_max2;
      double f_max_diff;
      double f_max_old;
      std::vector<unsigned int> num;
      double r_boxes;

    private:
      size_t max_;
      /// List of grid coordinates
      std::vector<coord_t> n_map_;
      /// Maximal value of the function at one given point
      std::vector<double> f_max_;
      /// Maximal value of the function in the considered integration range
      double f_max_global_;
  };
}

#endif
