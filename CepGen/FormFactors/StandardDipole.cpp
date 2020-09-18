#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include <cmath>

namespace cepgen
{
  namespace formfac
  {
    class StandardDipole : public Parameterisation
    {
      public:
        using Parameterisation::Parameterisation;
        static std::string description() { return "Standard dipole"; }

      private:
        void compute( double q2 ) override {
          GE = pow( 1.+q2/0.71, -2. );
          GM = MU*GE;
        }
    };
  }
}

REGISTER_FF_MODEL( "StandardDipole", StandardDipole )
