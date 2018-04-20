#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include <math.h>
#include <array>

#include "CepGen/Core/utils.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Physics/Constants.h"

namespace CepGen
{
  /// Form factors collection (electric and magnetic parts)
  class FormFactors
  {
    public:
      /// Initialise a collection of electric/magnetic form factors
      FormFactors( double fe = 0., double fm = 0. ) : FE( fe ), FM( fm ) {}
      // compute x from w2/m2
      double x( double q2, double w2, double m2 = 0. ) const;
      /// Trivial, spin-0 form factors (e.g. pion)
      static FormFactors Trivial();
      /// Elastic proton form factors
      static FormFactors ProtonElastic( double q2 );
      /// Suri-Yennie inelastic form factors
      static FormFactors SuriYennie( double q2, double mi2, double mf2 );
      /// Brasse et al. inelastic form factors
      /// \cite Brasse1976413
      static FormFactors FioreBrasse( double q2, double mi2, double mf2 );
      /// Szczurek-Uleschenko inelastic form factors
      static FormFactors SzczurekUleshchenko( double q2, double mi2, double mf2 );
      /// Generate the form factors according to the proton structure functions set
      static FormFactors ProtonInelastic( const StructureFunctions::Type& sf, double q2, double mi2, double mf2 );

      /// Electric form factor
      double FE;
      /// Magnetic form factor
      double FM;
      /// Dumping operator for standard output streams
      friend std::ostream& operator<<( std::ostream&, const FormFactors& );

    private:
      static const double mp_, mp2_;
  };
}

#endif
