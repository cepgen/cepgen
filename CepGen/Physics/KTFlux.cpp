#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Core/Exception.h"

namespace CepGen
{
  const double KTFluxParameters::kMinKTFlux = 1.e-20;
  const double KTFluxParameters::kMP = ParticleProperties::mass( PDG::proton );
  const double KTFluxParameters::kMP2 = KTFluxParameters::kMP*KTFluxParameters::kMP;

  double
  ktFlux( const KTFlux& type, double x, double kt2, StructureFunctions& sf, double mx )
  {
    double flux = 0.;
    const double mp2 = KTFluxParameters::kMP2;
    switch ( type ) {
      case KTFlux::P_Photon_Elastic: {
        const double x2 = x*x;
        const double q2min = x2*mp2/( 1.-x ), q2 = q2min + kt2/( 1.-x );

        // electromagnetic form factors
        const auto& ff = FormFactors::protonElastic( q2 );

        const double ela1 = ( 1.-x )*( 1.-q2min/q2 );
        //const double ela3 = 1.-( q2-kt2 )/q2;

        flux = Constants::alphaEM*M_1_PI*( 1.-x )/q2*( ela1*ff.FE + 0.5*x2*ff.FM );
      } break;
      case KTFlux::P_Photon_Inelastic_Budnev: {
        const double mx2 = mx*mx, x2 = x*x;
        const double q2min = ( x*( mx2-mp2 ) + x2*mp2 )/( 1.-x ), q2 = q2min + kt2/( 1.-x );
        const double xbj = q2 / ( q2+mx2-mp2 );

        // structure functions
        auto& str_fun = sf( xbj, q2 );
        str_fun.computeFL( xbj, q2 );

        const double f_D = str_fun.F2/( q2+mx2-mp2 )* ( 1.-x )*( 1.-q2min/q2 );
        const double f_C = str_fun.F1( xbj, q2 ) * 2./q2;

        flux = Constants::alphaEM*M_1_PI*( 1.-x )/q2*( f_D+0.5*x2*f_C );
      } break;
      default:
        throw CG_FATAL( "GenericKTProcess:flux" ) << "Invalid flux type: " << type;
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
      case KTFlux::invalid: default:
        return os << "unrecognized flux (" << (int)type << ")";
    }
  }
}
