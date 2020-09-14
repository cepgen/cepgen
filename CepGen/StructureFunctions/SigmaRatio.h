#ifndef CepGen_StructureFunctions_SigmaRatio_h
#define CepGen_StructureFunctions_SigmaRatio_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen
{
  /// A collector namespace for modellings of the \f$R=\sigma_L/\sigma_T\f$ ratio
  namespace sigrat
  {
    /// \f$R=\sigma_L/\sigma_T\f$ ratio modelling type
    enum struct Type { Invalid = 0, E143 = 1, R1990 = 2, CLAS = 3, SibirtsevBlunden = 4 };
    /// A generic modelling of the \f$R=\sigma_L/\sigma_T\f$ ratio
    class Parameterisation : public NamedModule<int>
    {
      public:
        /// \f$R=\sigma_L/\sigma_T\f$ ratio computation algorithm constructor
        explicit Parameterisation( const ParametersList& params = ParametersList() );
        /// Extract the longitudinal/transverse cross section ratio and associated error for a given \f$(x_{\rm Bj},Q^2)\f$ couple.
        virtual double operator()( double xbj, double q2, double& err ) const = 0;

      protected:
        /// \f$x_{\rm Bj}\f$ dependence for QCD-matching of R at high-\f$Q^2\f$
        double theta( double xbj, double q2 ) const;
        const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
    };
  }
  /// Human-readable description of this R-ratio parameterisation type
  std::ostream& operator<<( std::ostream&, const sigrat::Type& );
}

#endif
