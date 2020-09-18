#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/Modes.h"

namespace cepgen
{
  class ParametersList;
  namespace strfun { class Parameterisation; }
  namespace formfac
  {
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
        friend std::ostream& operator<<( std::ostream&, const Parameterisation* );
        /// Dumping operator for standard output streams
        friend std::ostream& operator<<( std::ostream&, const Parameterisation& );

        /// Specify the structure functions modelling where applicable
        void setStructureFunctions( strfun::Parameterisation* );
        strfun::Parameterisation* structureFunctions() const { return str_fun_; }

        double tau( double q2 ) const;

        /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
        Parameterisation& operator()( const mode::Beam& /*type*/, double /*q2*/, double mi2 = 0., double mf2 = 0. );

      protected:
        static constexpr double MU = 2.79;

        virtual void compute( double ) {}

        const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

      private:
        strfun::Parameterisation* str_fun_;
        double last_q2_;

      public:
        double FE; ///< Electric form factor
        double FM; ///< Magnetic form factor

        double GE;
        double GM;
    };
    /// Human-readable dump of the form factor parameterisation
    std::ostream& operator<<( std::ostream&, const formfac::Parameterisation* );
  }
}

#endif
