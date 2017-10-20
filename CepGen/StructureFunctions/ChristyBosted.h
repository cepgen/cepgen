#ifndef CepGen_StructureFunctions_ChristyBosted_h
#define CepGen_StructureFunctions_ChristyBosted_h

#include "StructureFunctions.h"
#include "CepGen/Physics/Constants.h"
#include <array>
#include <vector>

namespace CepGen
{
  namespace SF
  {
    class ChristyBosted : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard();

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
            double kr() const;
            double ecmr( double m2 ) const;
            double kcmr() const { return ecmr( 0. ); }
            double pcmr( double m2 ) const { return sqrt( std::max( 0., ecmr( m2 )*ecmr( m2 )-m2 ) ); }
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

        ChristyBosted( const ChristyBosted::Parameterisation& params = ChristyBosted::Parameterisation::standard() ) : params_( params ) {}

        ChristyBosted operator()( double q2, double xbj ) const;

      private:
        double resmod507( char sf, double w2, double q2 ) const;
        Parameterisation params_;
    };
  }
}

#endif
