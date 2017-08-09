#ifndef CepGen_Physics_FormFactors_h
#define CepGen_Physics_FormFactors_h

#include <math.h>

#include "CepGen/Core/utils.h"
#include "Constants.h"
#include "Particle.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& ); // xpart,q2part,uv,dv,us,ds,ss,wg
}

namespace CepGen
{
  /**
   * Compute the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
   *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
   * \param[in] q2 Squared 4-momentum transfer
   * \param[in] mx2 Squared mass of the proton remnant
   * \param[out] sigma_t ...
   * \param[out] w1 First proton structure function: \f$\mathcal W_1\f$
   * \param[out] w2 Second proton structure function: \f$\mathcal W_2\f$
   * \cite Brasse1976413
   */
  bool PSF( double q2, double mx2, double& sigma_t, double& w1, double& w2 );

  /// Form factors collection (electric and magnetic parts)
  struct FormFactors {
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
