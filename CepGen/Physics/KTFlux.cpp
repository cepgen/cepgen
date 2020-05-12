#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/GluonGrid.h"

#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Timer.h"

namespace
{
  extern "C"
  {
    void f_inter_kmr_fg_( double& logx, double& logkt2, double& logmu2, int& mode, double& fg );
  }
}

namespace cepgen
{
  const double KTFluxParameters::kMinKTFlux = 1.e-20;

  double
  ktFlux( const KTFlux& type, double x, double kt2, strfun::Parameterisation& sf, double mi2, double mf2 )
  {
    switch ( type ) {
      case KTFlux::P_Photon_Elastic:
      case KTFlux::P_Photon_Elastic_Budnev: {
        const double x2 = x*x;
        const double q2min = x2*mi2/( 1.-x ), q2 = q2min + kt2/( 1.-x );
        //--- proton electromagnetic form factors
        const auto& ff = FormFactors::protonElastic( q2 );
        if ( type == KTFlux::P_Photon_Elastic )
          return constants::ALPHA_EM*M_1_PI*ff.FE*std::pow( kt2/( kt2+x2*mi2 ), 2 )/q2;
        else {
          const double f_D = ff.FE*( 1.-x )*( 1.-q2min/q2 );
          const double f_C = ff.FM;
          return constants::ALPHA_EM*M_1_PI*( 1.-x )/q2*( f_D+0.5*x2*f_C );
        }
      } break;
      case KTFlux::P_Photon_Inelastic:
      case KTFlux::P_Photon_Inelastic_Budnev: {
        const double x2 = x*x;
        const double q2min = ( x*( mf2-mi2 ) + x2*mi2 )/( 1.-x );
        const double q2 = q2min + kt2/( 1.-x );
        const double denom = 1./( q2+mf2-mi2 );
        const double xbj = denom*q2;
        //--- proton structure functions
        auto& str_fun = sf( xbj, q2 );
        if ( type == KTFlux::P_Photon_Inelastic ) {
          const double f_aux = str_fun.F2*denom*( 1.-( q2-kt2 )/q2 )*std::pow( kt2/( kt2+x*( mf2-mi2 )+x2*mi2 ), 2 );
          return constants::ALPHA_EM*M_1_PI*( 1.-x )*f_aux/kt2;
        }
        else {
          str_fun.computeFL( xbj, q2 );
          const double f_D = str_fun.F2*denom*( 1.-x )*( 1.-q2min/q2 );
          const double f_C = str_fun.F1( xbj, q2 ) * 2./q2;
          return constants::ALPHA_EM*M_1_PI*( 1.-x )/q2*( f_D+0.5*x2*f_C );
        }
      } break;
      case KTFlux::P_Gluon_KMR_legacy: {
        static bool built = false;
        double lx = log10( x ), lkt2 = log10( kt2 ), lmx2 = log10( mf2 ), fg;
        if ( !built ) {
          utils::Timer tmr;
          CG_INFO( "KTFlux:KMR_legacy" )
            << "Building the legacy KMR interpolation grid.";
          int zero = 0;
          f_inter_kmr_fg_( lx, lkt2, lmx2, zero, fg );
          CG_INFO( "KTFlux:KMR_legacy" )
            << "Legacy KMR interpolation grid built in " << tmr.elapsed() << " s.";
          built = true;
        }
        int one = 1;
        f_inter_kmr_fg_( lx, lkt2, lmx2, one, fg );
        return fg;
      }
      case KTFlux::P_Gluon_KMR: {
        return kmr::GluonGrid::get()( x, kt2, mf2 );
      } break;
      default:
        throw CG_FATAL( "KTFlux" ) << "Invalid flux type: " << type;
    }
  }

  double
  ktFlux( const KTFlux& type, double x, double kt2, const HeavyIon& hi )
  {
    const double& mp = PDG::get().mass( PDG::proton );
    double flux = 0.;
    switch ( type ) {
      case KTFlux::HI_Photon_Elastic: {
        const double r_a = 1.1*cbrt( hi.A ), a0 = 0.7, m_a = hi.A*mp;
        const double q2_ela = ( kt2+x*x*m_a*m_a )/( 1.-x ), cons = sqrt( q2_ela )/0.1973;
        const double tau = cons*r_a, tau1 = cons*a0;
        // "Realistic nuclear form-factor" as used in STARLIGHT
        const double ff1 = 3.*( sin( tau )-tau*cos( tau ) )/pow( tau+1.e-10, 3 );
        const double ff2 = 1./( 1.+tau1*tau1 );
        const double ela1 = pow( kt2/( kt2+x*x*m_a*m_a ), 2 );
        const double ela2 = pow( ff1*ff2, 2 )/*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
        const unsigned int z = (unsigned short)hi.Z;
        flux = constants::ALPHA_EM*M_1_PI*z*z*ela1*ela2/q2_ela;
      } break;
      default:
        throw CG_FATAL("KTFlux") << "Invalid flux type: " << type;
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
      case KTFlux::P_Photon_Elastic_Budnev:
        return os << "elastic photon from proton (Budnev)";
      case KTFlux::P_Photon_Inelastic:
        return os << "inelastic photon from proton";
      case KTFlux::P_Photon_Inelastic_Budnev:
        return os << "inelastic photon from proton (Budnev)";
      case KTFlux::P_Gluon_KMR:
        return os << "elastic gluon from proton (KMR)";
      case KTFlux::HI_Photon_Elastic:
        return os << "elastic photon from HI";
      case KTFlux::invalid: default:
        return os << "unrecognised flux (" << (int)type << ")";
    }
  }
}
