#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include "CepGen/Core/utils.h"
#include "CepGen/Physics/Constants.h"

#include <memory>

namespace cepgen
{
  class ParametersList;
  namespace strfun { class Parameterisation; }
  namespace ff
  {
    enum struct Type {
      Invalid          = 0,
      ProtonElastic    = 1, ///< Elastic proton form factors
      PointLikeScalar  = 2, ///< Trivial, spin-0 form factors
      PointLikeFermion = 3, ///< Trivial, spin-1/2 form factors
      CompositeScalar  = 4, ///< Composite pion form factors
      ProtonInelastic  = 5, ///< Inelastic proton form factors (according to the proton structure functions set)
    };
    /// Proton form factors to be used in the outgoing state description
    enum struct Model {
      Invalid        = 0,
      StandardDipole = 1,
      ArringtonEtAl  = 2, ///< \cite Arrington:2007ux
      BrashEtAl      = 3, ///< \cite Brash:2001qq
    };
    /// Form factors parameterisation (electric and magnetic parts)
    class Parameterisation
    {
      public:
        explicit Parameterisation();
        explicit Parameterisation( const ParametersList& );
        /// Build a SF parameterisation for a given type
        static std::shared_ptr<Parameterisation> build( const ParametersList& params );
        /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
        Parameterisation& operator()( double /*q2*/, double mi2 = 0., double mf2 = 0. );
        /// Dumping operator for standard output streams
        friend std::ostream& operator<<( std::ostream&, const Parameterisation& );

      protected:
        static const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        static const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
        virtual void compute( double q2 ) {}
        const Model model_;

      private:
        const Type type_;
        std::shared_ptr<strfun::Parameterisation> str_fun_;

      public:
        double FE; ///< Electric form factor
        double FM; ///< Magnetic form factor

        double GE;
        double GM;
    };

    class StandardDipole : public Parameterisation
    {
      public:
        using Parameterisation::Parameterisation;
        void compute( double q2 ) override;

      private:
        static constexpr double MU = 2.79;
    };

    class ArringtonEtAl : public Parameterisation
    {
      public:
        ArringtonEtAl( const ParametersList& );
        void compute( double q2 ) override;

      private:
        static constexpr double MU = 2.79;
        const int mode_;
        std::vector<double> a_e_, b_e_;
        std::vector<double> a_m_, b_m_;
    };

    class BrashEtAl : public Parameterisation
    {
      public:
        using Parameterisation::Parameterisation;
        void compute( double q2 ) override;

      private:
        static constexpr double MU = 2.79;
    };
  }
}

#endif
