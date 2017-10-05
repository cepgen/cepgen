#ifndef CepGen_StructureFunctions_ChristyBosted_h
#define CepGen_StructureFunctions_ChristyBosted_h

#include "StructureFunctions.h"
#include "CepGen/Physics/Constants.h"
#include <array>

namespace CepGen
{
  namespace SF
  {
    
    struct ChristyBostedParameterisation
    {
      static ChristyBostedParameterisation standard();

      struct ResonanceParameters
      {
        struct BranchingRatios
        {
          BranchingRatios() : singlepi( 0. ), doublepi( 0. ), eta( 0. ) {}
          BranchingRatios( double singlepi, double doublepi, double eta ) : singlepi( singlepi ), doublepi( doublepi ), eta( eta ) {}
          bool valid() const { return ( singlepi+doublepi+eta == 1. ); }
          /// single pion branching ratio
          double singlepi;
          /// double pion branching ratio
          double doublepi;
          /// eta meson branching ratio
          double eta;
        };
        ResonanceParameters() : angular_momentum( 0. ), x0( 0. ), mass( 0. ), width( 0. ), A0_T( 0. ), A0_L( 0. ) {}
        double kr() const { return 0.5 * ( mass*mass-Constants::mp*Constants::mp )/Constants::mp; }
        double kcmr() const { return ( mass == 0. ) ? 0. : 0.5 * ( mass*mass-Constants::mp*Constants::mp ) / mass; }
        double epicmr() const { return ( mass == 0. ) ? 0. : 0.5 * ( mass*mass+Constants::mpi*Constants::mpi-Constants::mp*Constants::mp ) / mass; }
        double ppicmr() const { return sqrt( std::max( 0., epicmr()*epicmr()-Constants::mpi*Constants::mpi ) ); }
        double epi2cmr() const { return ( mass == 0. ) ? 0. : 0.5 * ( mass*mass+4.*Constants::mpi*Constants::mpi-Constants::mp*Constants::mp ) / mass; }
        double ppi2cmr() const { return sqrt( std::max( 0., epi2cmr()*epi2cmr()-4.*Constants::mpi*Constants::mpi ) ); }
        double eetacmr() const {
          const double meta = 0.547862;
          return ( mass == 0. ) ? 0. : 0.5 * ( mass*mass+meta*meta-Constants::mp*Constants::mp ) / mass;
        }
        double petacmr() const {
          const double meta = 0.547862;
          return sqrt( std::max( 0., eetacmr()*eetacmr()-meta*meta ) );
        }
        BranchingRatios br;
        /// meson angular momentum
        double angular_momentum;
        /// damping parameter
        double x0;
        /// mass, in GeV/c2
        double mass;
        /// full width, in GeV
        double width;
        double A0_T;
        double A0_L;
        std::array<double,5> fit_parameters;
      };
      struct ContinuumParameters
      {
        struct DirectionParameters
        {
          DirectionParameters() : sig0( 0. ) {}
          DirectionParameters( double sig0, const std::vector<double>& params ) : sig0( sig0 ), fit_parameters( params ) {}
          double sig0;
          std::vector<double> fit_parameters;
        };
        std::array<DirectionParameters,2> transverse;
        std::array<DirectionParameters,1> longitudinal;
      };
      double m0;
      std::vector<ResonanceParameters> resonances;
      ContinuumParameters continuum;
    };
    StructureFunctions ChristyBosted( double q2, double xbj );
  }
}

#endif
