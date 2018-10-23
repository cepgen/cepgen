#ifndef CepGen_StructureFunctions_ALLM_h
#define CepGen_StructureFunctions_ALLM_h

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include <vector>

namespace cepgen
{
  namespace strfun
  {
    /// \f$F_{2,L}\f$ parameterisation by Abramowicz, Levin, Levy, and Maor \cite Abramowicz:1991xz\cite Abramowicz:1997ms
    class ALLM : public Parameterisation
    {
      public:
        class Parameters
        {
          private:
            struct Trajectory
            {
              Trajectory( const ParametersList& params = ParametersList() );
              std::vector<double> a, b, c;
            };

          public:
            Parameters( const ParametersList& params = ParametersList() );
            /// Pre-HERA data fit (694 data points)
            static Parameters allm91();
            /// Fixed target and HERA photoproduction total cross sections (1356 points)
            static Parameters allm97();
            static Parameters hht_allm();
            static Parameters hht_allm_ft();
            static Parameters gd07p();
            static Parameters gd11p();

            Trajectory pomeron, reggeon;
            /// Effective photon squared mass
            double m02;
            /// Effective pomeron squared mass
            double mp2;
            /// Effective reggeon squared mass
            double mr2;
            double q02;
            /// Squared QCD scale
            double lambda2;
            Type type;
        };

        explicit ALLM( const ParametersList& params = ParametersList() );
        ALLM& operator()( double xbj, double q2 ) override;

      private:
        Parameters params_;
    };
  }
}

#endif
