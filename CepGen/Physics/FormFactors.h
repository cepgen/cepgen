#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include <math.h>
#include <array>

#include "CepGen/Core/utils.h"
#include "CepGen/Physics/Constants.h"

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  /// Form factors collection (electric and magnetic parts)
  class FormFactors
  {
    public:
      /// Initialise a collection of electric/magnetic form factors
      FormFactors( double fe = 0., double fm = 0. );
      /// Trivial, spin-0 form factors
      static FormFactors pointlikeScalar();
      /// Trivial, spin-1/2 form factors
      static FormFactors pointlikeFermion();
      /// Composite pion form factors
      static FormFactors compositeScalar( double q2 );
      /// Elastic proton form factors
      static FormFactors protonElastic( double q2 );
      /// Generate the form factors according to the proton structure functions set
      static FormFactors protonInelastic( double q2, double mi2, double mf2, strfun::Parameterisation& );

      /// Electric form factor
      double FE;
      /// Magnetic form factor
      double FM;
      /// Dumping operator for standard output streams
      friend std::ostream& operator<<( std::ostream&, const FormFactors& );
  };
}

#endif
