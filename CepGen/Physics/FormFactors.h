#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include <math.h>

#include "CepGen/Core/utils.h"

#include "Constants.h"
#include "StructureFunctions.h"

namespace CepGen
{
  /// Form factors collection (electric and magnetic parts)
  struct FormFactors {
    /// Initialise a collection of electric/magnetic form factors
    FormFactors( double fe=0.0, double fm=0.0 ) : FE( fe ), FM( fm ) {}
    /// Electric form factor
    double FE;
    /// Magnetic form factor
    double FM;
    /// Dumping operator for standard output streams
    friend std::ostream& operator<<( std::ostream&, const FormFactors& );
  };

  /// Trivial, spin-0 form factors (e.g. pion)
  FormFactors TrivialFormFactors();
  /// Elastic form factors
  FormFactors ElasticFormFactors( double q2, double mi2 );
  /// Suri-Yennie inelastic form factors
  FormFactors SuriYennieFormFactors( double q2, double mi2, double mf2 );
  /// Brasse et al. inelastic form factors
  /// \cite Brasse1976413
  FormFactors FioreBrasseFormFactors( double q2, double mi2, double mf2 );
  /// Szczurek-Uleschenko inelastic form factors
  FormFactors SzczurekUleshchenkoFormFactors( double q2, double mi2, double mf2 );
}

#endif
