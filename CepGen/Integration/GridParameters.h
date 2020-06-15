#ifndef CepGen_Core_GridParameters_h
#define CepGen_Core_GridParameters_h

#include <vector>

namespace cepgen
{
  class Integrator;
  /// A parameters placeholder for the grid integration helper
  class GridParameters
  {
    public:
      GridParameters( size_t ndim );

      /// Coordinates definition
      typedef std::vector<unsigned short> coord_t;

      /// Dump the grid coordinates
      void dump() const;

      /// Grid multiplicity
      size_t size() const { return max_; }
      const coord_t& n( size_t coord ) const;
      /// Global function maximum
      double globalMax() const { return f_max_global_; }
      /// Maximal function value for a given grid coordinate
      double maxValue( size_t coord ) const;
      /// Set the function value for a given grid coordinate
      void setValue( size_t coord, double val );
      /// Shoot a phase space point for a grid coordinate
      void shoot( const Integrator* integ, size_t coord, std::vector<double>& out ) const;
      /// Specify a new trial has been attempted for bin
      void increment( size_t coord );
      /// Number of points already shot for a given grid coordinate
      size_t numPoints( size_t coord ) const;

      /// Maximal number of dimensions handled by this integrator instance
      static constexpr unsigned short MAX_DIM = 15;
      /// Integration grid size parameter
      static constexpr unsigned short M_BIN = 3;
      /// Weight of each grid coordinate
      static constexpr double INV_M_BIN = 1./M_BIN;

      /// Has the grid been already prepared?
      bool gen_prepared;
      double correc;
      double correc2;
      double f_max2;
      double f_max_diff;
      double f_max_old;

    private:
      size_t max_;
      /// List of grid coordinates
      std::vector<coord_t> n_map_;
      /// Number of functions values evaluated for this point
      std::vector<size_t> num_points_;
      /// Maximal value of the function at one given point
      std::vector<double> f_max_;
      /// Maximal value of the function in the considered integration range
      double f_max_global_;
  };
}

#endif
