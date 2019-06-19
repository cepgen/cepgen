#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include "CepGen/Core/utils.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Physics/Constants.h"

#include <memory>

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  namespace ff
  {
    /// Proton form factors to be used in the outgoing state description
    enum struct Type {
      StandardDipole = 0,
      ArringtonEtAl  = 1, ///< \cite Arrington:2007ux
      BrashEtAl      = 2, ///< \cite Brash:2001qq
    };
    class Parameterisation
    {
      public:
        explicit Parameterisation( double q2 );

        /// Build a SF parameterisation for a given type
        static std::shared_ptr<Parameterisation> build( const Type& type, const ParametersList& params = ParametersList() );
        /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
        virtual Parameterisation& operator()( double /*q2*/ ) { return *this; }

        double GE;
        double GM;
    };

    class StandardDipole : public Parameterisation
    {
      public:
        explicit StandardDipole();
        StandardDipole& operator()( double q2 );

      private:
        static constexpr double MU = 2.79;
    };
  }

  /// Form factors collection (electric and magnetic parts)
  class FormFactors
  {
    public:
      /// Initialise a collection of electric/magnetic form factors
      FormFactors( double fe = 0., double fm = 0. ) : FE( fe ), FM( fm ) {}
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

    private:
      static const double mp_, mp2_;
  };
}

#endif
