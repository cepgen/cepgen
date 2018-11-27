#ifndef CepGen_StructureFunctions_ChristyBosted_h
#define CepGen_StructureFunctions_ChristyBosted_h

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Physics/Constants.h"

#include <array>
#include <vector>

namespace cepgen
{
  namespace strfun
  {
    /// \f$F_{2,L}\f$ parameterisation by Christy and Bosted \cite Bosted:2007xd
    class ChristyBosted : public Parameterisation
    {
      public:
        struct Parameters
        {
          static Parameters standard();

          struct Resonance
          {
            /// Branching ratios container for resonance decay into single, double pion or eta states
            struct BR
            {
              BR() : singlepi( 0. ), doublepi( 0. ), eta( 0. ) {}
              BR( double singlepi, double doublepi, double eta ) : singlepi( singlepi ), doublepi( doublepi ), eta( eta ) {}
              bool valid() const { return ( singlepi+doublepi+eta == 1. ); }
              /// single pion branching ratio
              double singlepi;
              /// double pion branching ratio
              double doublepi;
              /// eta meson branching ratio
              double eta;
            };
            Resonance() : angular_momentum( 0. ), x0( 0. ), mass( 0. ), width( 0. ), A0_T( 0. ), A0_L( 0. ) {}
            double kr() const;
            double ecmr( double m2 ) const;
            double kcmr() const { return ecmr( 0. ); }
            double pcmr( double m2 ) const { return sqrt( std::max( 0., ecmr( m2 )*ecmr( m2 )-m2 ) ); }
            BR br;
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
          struct Continuum
          {
            struct Direction
            {
              Direction() : sig0( 0. ) {}
              Direction( double sig0, const std::vector<double>& params ) : sig0( sig0 ), fit_parameters( params ) {}
              double sig0;
              std::vector<double> fit_parameters;
            };
            std::array<Direction,2> transverse;
            std::array<Direction,1> longitudinal;
          };
          double m0;
          std::vector<Resonance> resonances;
          Continuum continuum;
        };

        explicit ChristyBosted( const ParametersList& params = ParametersList() );
        ChristyBosted& operator()( double xbj, double q2 ) override;

        //--- already computed internally during F2 computation
        ChristyBosted& computeFL( double xbj, double q2 ) override { return *this; }
        ChristyBosted& computeFL( double xbj, double q2, double r ) override { return *this; }

      private:
        double resmod507( char sf, double w2, double q2 ) const;
        Parameters params_;
    };
  }
}

#endif
