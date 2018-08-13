#ifndef CepGen_StructureFunctions_Schaefer_h
#define CepGen_StructureFunctions_Schaefer_h

#include "StructureFunctions.h"
#include <memory>

namespace CepGen
{
  namespace SF
  {
    class Schaefer : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard();
          static Parameterisation cteq();
          double q2_cut, w2_lo, w2_hi;
          std::shared_ptr<StructureFunctions> resonances_model, perturbative_model, continuum_model;
          bool higher_twist;
        };
        Schaefer( const Parameterisation& param = Parameterisation::standard() );
        Schaefer& operator()( double q2, double xbj ) override;

        Parameterisation params;

      private:
        double rho( double w2 ) const;
        void initialise();
        bool initialised_;
        double inv_omega_range_;
    };
  }
}

#endif

