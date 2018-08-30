#ifndef CepGen_Physics_GluonGrid_h
#define CepGen_Physics_GluonGrid_h

#include "CepGen/IO/GridHandler.h"

#define DEFAULT_KMR_GRID_PATH "gluon_mmht2014nlo_Watt.dat"

/// Kimber-Martin-Ryskin unintegrated gluon densities
namespace kmr
{
  /// A KMR unintegrated gluon densities grid interpolator
  class GluonGrid : private CepGen::GridHandler<3,1>
  {
    public:
      struct Parameterisation {
        Parameterisation() : grid_path( DEFAULT_KMR_GRID_PATH ) {}
        std::string grid_path;
      };

    public:
      /// Retrieve the grid interpolator (singleton)
      static GluonGrid& get( const char* path = DEFAULT_KMR_GRID_PATH );

      /// Compute the gluon flux
      double operator()( double x, double kt2, double mu2 ) const;
      Parameterisation params;

    public:
      GluonGrid( const GluonGrid& ) = delete;
      void operator=( const GridHandler& ) = delete;

    private:
      explicit GluonGrid( const Parameterisation& = Parameterisation() );
  };
}

#undef DEFAULT_KMR_GRID_PATH

#endif

