#ifndef CepGen_Physics_GluonGrid_h
#define CepGen_Physics_GluonGrid_h

#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Core/ParametersList.h"

#define DEFAULT_KMR_GRID_PATH "gluon_mmht2014nlo_Watt.dat"

/// Kimber-Martin-Ryskin unintegrated gluon densities
namespace kmr
{
  /// A KMR unintegrated gluon densities grid interpolator
  class GluonGrid : private cepgen::GridHandler<3,1>
  {
    public:
      GluonGrid( const GluonGrid& ) = delete;
      void operator=( const GridHandler& ) = delete;

      /// Retrieve the grid interpolator (singleton)
      static GluonGrid& get( const char* path = DEFAULT_KMR_GRID_PATH );

      /// Compute the gluon flux
      double operator()( double x, double kt2, double mu2 ) const;

    private:
      explicit GluonGrid( const cepgen::ParametersList& = cepgen::ParametersList() );
      /// Location of the grid to be interpolated
      const std::string grid_path_;
  };
}

#endif
