#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include "CepGen/Modules/NamedModule.h"

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
    std::ostream& operator<<( std::ostream&, const Type& );
    /// Proton form factors to be used in the outgoing state description
    enum struct Model {
      Invalid        = 0,
      StandardDipole = 1,
      ArringtonEtAl  = 2, ///< \cite Arrington:2007ux
      BrashEtAl      = 3, ///< \cite Brash:2001qq
      MergellEtAl    = 4, ///< \cite Mergell:1995bf
    };
    std::ostream& operator<<( std::ostream&, const Model& );
    /// Form factors parameterisation (electric and magnetic parts)
    class Parameterisation : public NamedModule<int>
    {
      public:
        explicit Parameterisation();
        Parameterisation( const ParametersList& );
        Parameterisation( const Parameterisation& );

        static std::string description() { return "Unnamed form factors parameterisation"; }
        /// Dumping operator for standard output streams
        friend std::ostream& operator<<( std::ostream&, const Parameterisation& );

        /// Specify the structure functions modelling where applicable
        void setStructureFunctions( strfun::Parameterisation* );
        strfun::Parameterisation* structureFunctions() const { return str_fun_.get(); }

        const Type& type() const { return type_; }
        void setType( const Type& type ) { type_ = type; }

        double tau( double q2 ) const;

        /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
        Parameterisation& operator()( double /*q2*/, double mi2 = 0., double mf2 = 0. );

      protected:
        static constexpr double MU = 2.79;

        virtual void compute( double ) {}

        Type type_;

        const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

      private:
        std::shared_ptr<strfun::Parameterisation> str_fun_;
        double last_q2_;

      public:
        double FE; ///< Electric form factor
        double FM; ///< Magnetic form factor

        double GE;
        double GM;
    };
  }
}

#endif
