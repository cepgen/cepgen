#ifndef CepGen_Core_GridParameters_h
#define CepGen_Core_GridParameters_h

#include <vector>
#include <map>

namespace CepGen
{
  class GridParameters
  {
    public:
      /// Maximal number of dimensions handled by this integrator instance
      static constexpr unsigned short max_dimensions_ = 15;
      /// Integration grid size parameter
      static constexpr unsigned short mbin_ = 3;
      static constexpr double inv_mbin_ = 1./mbin_;

      GridParameters();

      std::map<unsigned int,std::vector<unsigned short> > n_map;

      unsigned int max;
      /// Has the generation been prepared?
      bool gen_prepared;
      /// Maximal value of the function at one given point
      std::vector<double> f_max;
      /// Maximal value of the function in the considered integration range
      double f_max_global;
      double f_max_diff;
      //std::set<unsigned int> probed_bins;
  };
}

#endif

