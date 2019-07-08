#ifndef CepGen_StructureFunctions_SzczurekUleshchenko_h
#define CepGen_StructureFunctions_SzczurekUleshchenko_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

namespace cepgen
{
  namespace strfun
  {
    /// Szcurek and Uleshchenko modelling of \f$F_{1,2}\f$ \cite Szczurek:1999wp
    class SzczurekUleshchenko : public Parameterisation
    {
      public:
        SzczurekUleshchenko( const ParametersList& params = ParametersList() );
        SzczurekUleshchenko& operator()( double xbj, double q2 ) override;

      private:
        /// \f$Q^2\f$ scale shift
        const float q2_shift_;

      public:
        double F1;
    };
  }
}

#endif
