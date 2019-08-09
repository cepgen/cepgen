#include "CepGen/Core/Process2to4.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/ProcessesHandler.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#include <assert.h>

namespace cepgen
{
  namespace proc
  {
    /// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
    /// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
    class PPtoWW : public Process2to4
    {
      public:
        PPtoWW( const ParametersList& params = ParametersList() );
        ProcessPtr clone( const ParametersList& params ) const override {
          return ProcessPtr( new PPtoWW( *this ) );
        }
        enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };

      private:
        static const double mw_, mw2_;

        void prepareKinematics() override;
        double computeCentralMatrixElement() const override;

        double amplitudeWW( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 ) const;
        double onShellME( double shat, double that, double uhat ) const;
        double offShellME( double shat, double that, double uhat, double phi_sum, double phi_diff ) const;

        const int method_;

        std::vector<short> pol_w1_, pol_w2_;
    };

    const double PPtoWW::mw_ = PDG::get().mass( PDG::W );
    const double PPtoWW::mw2_ = PPtoWW::mw_*PPtoWW::mw_;

    PPtoWW::PPtoWW( const ParametersList& params ) :
      Process2to4( params, "pptoww", "ɣɣ → W⁺W¯", { PDG::photon, PDG::photon }, PDG::W ),
      method_( params.get<int>( "method", 1 ) )
    {
      switch ( (Polarisation)params.get<int>( "polarisationStates", 0 ) ) {
        case Polarisation::LL: pol_w1_ = pol_w2_ = { 0 }; break;
        case Polarisation::LT: pol_w1_ = { 0 }; pol_w2_ = { -1, 1 }; break;
        case Polarisation::TL: pol_w1_ = { -1, 1 }; pol_w2_ = { 0 }; break;
        case Polarisation::TT: pol_w1_ = pol_w2_ = { -1, 1 }; break;
        default:
        case Polarisation::full: pol_w1_ = pol_w2_ = { -1, 0, 1 }; break;
      }
      CG_DEBUG( "PPtoWW:mode" )
        << "matrix element computation method: " << method_ << ".";
    }

    void
    PPtoWW::prepareKinematics()
    {
      Cuts single_w_cuts;
      if ( kin_.cuts.central_particles.count( PDG::W ) > 0 )
        single_w_cuts = kin_.cuts.central_particles.at( PDG::W );
      setCuts( single_w_cuts );
    }

    double
    PPtoWW::computeCentralMatrixElement() const
    {
      //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2_, 2))/s_*s_));

      //--- first compute a few Mendelstam variables
      const double shat = ( q1_+q2_ ).mass2();
      const double that1 = ( q1_-p_c1_ ).mass2(), that2 = ( q2_-p_c2_ ).mass2(), that = 0.5*( that1+that2 );
      const double uhat1 = ( q1_-p_c2_ ).mass2(), uhat2 = ( q2_-p_c1_ ).mass2(), uhat = 0.5*( uhat1+uhat2 );

      //--- matrix element computation

      CG_DEBUG_LOOP( "PPtoWW" )
        << "matrix element mode: " << method_ << ".";

      double amat2 = 0.;
      if ( method_ == 0 ) // on-shell matrix element
        // (Denner+Dittmaier+Schuster, + work in collaboration with C. Royon)
        amat2 = onShellME( shat, that, uhat );
      else if ( method_ == 1 ) // off-shell Nachtmann formulae
        amat2 = offShellME( shat, that, uhat, phi_qt1_+phi_qt2_, phi_qt1_-phi_qt2_ );

      const double g_em = 4.*M_PI*constants::ALPHA_EM;
      return std::max( g_em*g_em*amat2, 0. );
    }

    double
    PPtoWW::onShellME( double shat, double that, double uhat ) const
    {
      const double mw4 = mw2_*mw2_;

      const double term1 = 2.*shat * ( 2.*shat+3.*mw2_ ) / ( 3.*( mw2_-that )*( mw2_-uhat ) );
      const double term2 = 2.*shat*shat * ( shat*shat + 3.*mw4 ) / ( 3.*pow( mw2_-that, 2 )*pow( mw2_-uhat, 2 ) );

      return 6.*( 1.-term1+term2 );
    }

    double
    PPtoWW::offShellME( double shat, double that, double uhat, double phi_sum, double phi_diff ) const
    {
      double amat2_0 = 0., amat2_1 = 0., amat2_interf = 0.;
      for ( const auto lam3 : pol_w1_ )
        for ( const auto lam4 : pol_w2_ ) {
          double ampli_pp = amplitudeWW( shat, that, uhat, +1, +1, lam3, lam4 );
          double ampli_mm = amplitudeWW( shat, that, uhat, -1, -1, lam3, lam4 );
          double ampli_pm = amplitudeWW( shat, that, uhat, +1, -1, lam3, lam4 );
          double ampli_mp = amplitudeWW( shat, that, uhat, -1, +1, lam3, lam4 );

          amat2_0 += ampli_pp*ampli_pp + ampli_mm*ampli_mm + 2.*cos( 2.*phi_diff )*ampli_pp*ampli_mm;
          amat2_1 += ampli_pm*ampli_pm + ampli_mp*ampli_mp + 2.*cos( 2.*phi_sum  )*ampli_pm*ampli_mp;
          amat2_interf -= 2.*( cos( phi_sum+phi_diff )*( ampli_pp*ampli_pm+ampli_mm*ampli_mp )
                              +cos( phi_sum-phi_diff )*( ampli_pp*ampli_mp+ampli_mm*ampli_pm ) );
        }
      return amat2_0+amat2_1+amat2_interf;
    }

    double
    PPtoWW::amplitudeWW( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 ) const
    {
      //--- first compute some kinematic variables
      const double cos_theta = ( that-uhat ) / shat / sqrt( 1.+1.e-10-4.*mw2_/shat ),
                   cos_theta2 = cos_theta*cos_theta;
      const double sin_theta2 = 1.-cos_theta2,
                   sin_theta = sqrt( sin_theta2 );
      const double beta = sqrt( 1.-4.*mw2_/shat ), beta2 = beta*beta;
      const double inv_gamma = sqrt( 1.-beta2 ), gamma = 1./inv_gamma,
                   gamma2 = gamma*gamma, inv_gamma2 = inv_gamma*inv_gamma;
      const double invA = 1./( 1.-beta2*cos_theta2 );

      //--- per-helicity amplitude

      if ( lam3 == 0 && lam4 == 0 ) // longitudinal-longitudinal
        return invA*inv_gamma2*( ( gamma2+1. )*( 1.-lam1*lam2 )*sin_theta2 - ( 1.+lam1*lam2 ) );

      if ( lam4 == 0 )              // transverse-longitudinal
        return invA*( -M_SQRT2*inv_gamma*( lam1-lam2 )*( 1.+lam1*lam3*cos_theta )*sin_theta );

      if ( lam3 == 0 )              // longitudinal-transverse
        return invA*( -M_SQRT2*inv_gamma*( lam2-lam1 )*( 1.+lam2*lam4*cos_theta )*sin_theta );

      if ( lam3 != 0 && lam4 != 0 ) // transverse-transverse
        return -0.5*invA*( 2.*beta*( lam1+lam2 )*( lam3+lam4 )
                          -inv_gamma2*( 1.+lam3*lam4 )*( 2.*lam1*lam2+( 1.-lam1*lam2 ) * cos_theta2 )
                          +( 1.+lam1*lam2*lam3*lam4 )*( 3.+lam1*lam2 )
                          +2.*( lam1-lam2 )*( lam3-lam4 )*cos_theta
                          +( 1.-lam1*lam2 )*( 1.-lam3*lam4 )*cos_theta2 );
      return 0.;
    }
  }
}
// register process
REGISTER_PROCESS( "pptoww", PPtoWW )

