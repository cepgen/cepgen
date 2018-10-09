#ifndef CepGen_Core_GridParameters_h
#define CepGen_Core_GridParameters_h

#include <vector>
#include <map>

namespace cepgen
{
  /// A parameters placeholder for the grid integration helper
  class GridParameters
  {
    public:
      /// Maximal number of dimensions handled by this integrator instance
      static const unsigned short max_dimensions_;
      /// Integration grid size parameter
      static const unsigned short mbin_;
      static const double inv_mbin_;

      GridParameters();

      std::map<unsigned int,std::vector<unsigned short> > n_map;

      unsigned int max;
      /// Has the generation been prepared?
      bool gen_prepared;
      double correc;
      double correc2;
      /// Maximal value of the function at one given point
      std::vector<double> f_max;
      /// Maximal value of the function in the considered integration range
      double f_max_global;
      double f_max2;
      double f_max_diff;
      double f_max_old;
      std::vector<unsigned int> num;
  };
}

#endif

