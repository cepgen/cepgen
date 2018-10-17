#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/GluonGrid.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Core/Exception.h"

namespace cepgen
{
  const double KTFluxParameters::kMinKTFlux = 1.e-20;
  const double KTFluxParameters::kMP = particleproperties::mass( PDG::proton );
  const double KTFluxParameters::kMP2 = KTFluxParameters::kMP*KTFluxParameters::kMP;

  double
  ktFlux( const KTFlux& type, double x, double kt2, strfun::Parameterisation& sf, double mx )
  {
    double flux = 0.;
    const double mp2 = KTFluxParameters::kMP2;
    switch ( type ) {
      case KTFlux::P_Photon_Elastic: {
        const double x2 = x*x;
        const double q2min = x2*mp2/( 1.-x ), q2 = q2min + kt2/( 1.-x );
        //--- proton electromagnetic form factors
        const auto& ff = FormFactors::protonElastic( q2 );
        flux = constants::alphaEM*M_1_PI/( 1.-x )/q2*( ( 1.-x )*( 1.-q2min/q2 )*ff.FE + 0.25*x2*ff.FM );
      } break;
      case KTFlux::P_Photon_Inelastic_Budnev: {
        const double mx2 = mx*mx, x2 = x*x;
        const double q2min = ( x2*mp2+x*( mx2-mp2 ) )/( 1.-x ), q2 = q2min + kt2/( 1.-x );
        const double xbj = q2 / ( q2+mx2-mp2 );
        //--- proton structure functions
        auto& str_fun = sf( xbj, q2 );
        str_fun.computeFL( xbj, q2 );
        const double f_D = str_fun.F2/( q2+mx2-mp2 )*( 1.-x )*( 1.-q2min/q2 );
        const double f_C = str_fun.F1( xbj, q2 ) * 2./q2;
        flux = constants::alphaEM*M_1_PI*( 1.-x )/q2*( f_D+0.5*x2*f_C );
      } break;
      case KTFlux::P_Gluon_KMR: {
        flux = kmr::GluonGrid::get()( log10( x ), log10( kt2 ), 2.*log10( mx ) );
      } break;
      default:
        throw CG_FATAL( "GenericKTProcess:flux" ) << "Invalid flux type: " << type;
    }
    if ( flux < KTFluxParameters::kMinKTFlux )
      return 0.;
    return flux;
  }

  double
  ktFlux( const KTFlux& type, double x, double kt2, const HeavyIon& hi )
  {
    double flux = 0.;
    switch ( type ) {
      case KTFlux::HI_Photon_Elastic: {
        const double r_a = 1.1*std::pow( hi.A, 1./3 ), a0 = 0.7, m_a = hi.A*KTFluxParameters::kMP;
        const double q2_ela = ( kt2+x*x*m_a*m_a )/( 1.-x ), cons = sqrt( q2_ela )/0.1973;
        const double tau = cons*r_a, tau1 = cons*a0;
        // "Realistic nuclear form-factor" as used in STARLIGHT
        const double ff1 = 3.*( sin( tau )-tau*cos( tau ) )/pow( tau+1.e-10, 3 );
        const double ff2 = 1./( 1.+tau1*tau1 );
        const double ela1 = pow( kt2/( kt2+x*x*m_a*m_a ), 2 );
        const double ela2 = pow( ff1*ff2, 2 )/*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
        const unsigned int z = (unsigned short)hi.Z;
        flux = constants::alphaEM*M_1_PI*z*z*ela1*ela2/q2_ela;
      } break;
      default:
        throw CG_FATAL("GenericKTProcess:flux") << "Invalid flux type: " << type;
    }
    if ( flux < KTFluxParameters::kMinKTFlux )
      return 0.;
    return flux;
  }

  std::ostream&
  operator<<( std::ostream& os, const KTFlux& type )
  {
    switch ( type ) {
      case KTFlux::P_Photon_Elastic:
        return os << "elastic photon from proton";
      case KTFlux::P_Photon_Inelastic:
        return os << "inelastic photon from proton";
      case KTFlux::P_Photon_Inelastic_Budnev:
        return os << "inelastic photon from proton (Budnev)";
      case KTFlux::P_Gluon_KMR:
        return os << "elastic gluon from proton (KMR)";
      case KTFlux::HI_Photon_Elastic:
        return os << "elastic photon from HI";
      case KTFlux::invalid: default:
        return os << "unrecognized flux (" << (int)type << ")";
    }
  }
}
