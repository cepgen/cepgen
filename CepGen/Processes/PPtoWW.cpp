#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include <assert.h>

#include "CepGen/Processes/ProcessesHandler.h"

namespace CepGen
{
  namespace process
  {
    const double PPtoWW::mw_ = part::mass( PDG::W );
    const double PPtoWW::mw2_ = PPtoWW::mw_*PPtoWW::mw_;

    PPtoWW::PPtoWW( const ParametersList& params ) :
      GenericKTProcess( params, "pptoww", "ɣɣ → W⁺W¯", { { PDG::photon, PDG::photon } }, { PDG::W, PDG::W } ),
      method_( params.get<int>( "method", 1 ) ),
      pol_state_( (Polarisation)params.get<int>( "polarisationStates", 0 ) ),
      y1_( 0. ), y2_( 0. ), pt_diff_( 0. ), phi_pt_diff_( 0. )
    {}

    void
    PPtoWW::preparePhaseSpace()
    {
      registerVariable( y1_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "First outgoing W rapidity" );
      registerVariable( y2_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "Second outgoing W rapidity" );
      registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central.pt_diff, { 0., 500. }, "Ws transverse momentum difference" );
      registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central.phi_pt_diff, { 0., 2.*M_PI }, "Ws azimuthal angle difference" );

      switch ( pol_state_ ) {
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

    double
    PPtoWW::computeKTFactorisedMatrixElement()
    {
      //=================================================================
      //     matrix element computation
      //=================================================================
      //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2_, 2))/s_*s_));

      // Inner photons
      const double q1tx = qt1_*cos( phi_qt1_ ), q1ty = qt1_*sin( phi_qt1_ ),
                   q2tx = qt2_*cos( phi_qt2_ ), q2ty = qt2_*sin( phi_qt2_ );
      CG_DEBUG_LOOP( "PPtoWW:qt" )
        << "q1t(x/y) = " << q1tx << " / " << q1ty << "\n\t"
        << "q2t(x/y) = " << q2tx << " / " << q2ty << ".";

      // Two-photon system
      const double ptsumx = q1tx+q2tx,
                   ptsumy = q1ty+q2ty,
                   ptsum = sqrt( ptsumx*ptsumx+ptsumy*ptsumy );

      const double ptdiffx = pt_diff_*cos( phi_pt_diff_ ),
                   ptdiffy = pt_diff_*sin( phi_pt_diff_ );

      // Outgoing leptons
      const double pt1x = ( ptsumx+ptdiffx )*0.5, pt1y = ( ptsumy+ptdiffy )*0.5, pt1 = std::hypot( pt1x, pt1y ),
                   pt2x = ( ptsumx-ptdiffx )*0.5, pt2y = ( ptsumy-ptdiffy )*0.5, pt2 = std::hypot( pt2x, pt2y );

      if ( cuts_.cuts.central_particles.count( PDG::W ) > 0
        && cuts_.cuts.central_particles.at( PDG::W ).pt_single.valid() ) {
        const Limits pt_limits = cuts_.cuts.central_particles.at( PDG::W ).pt_single;
        if ( !pt_limits.passes( pt1 ) || !pt_limits.passes( pt2 ) )
          return 0.;
      }

      // transverse mass for the two leptons
      const double amt1 = sqrt( pt1*pt1+mw2_ ),
                   amt2 = sqrt( pt2*pt2+mw2_ );

      //=================================================================
      //     a window in two-boson invariant mass
      //=================================================================

      const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh( y1_-y2_ ) - ptsum*ptsum );
      if ( !cuts_.cuts.central.mass_sum.passes( invm ) )
        return 0.;

      //=================================================================
      //     a window in transverse momentum difference
      //=================================================================

      if ( !cuts_.cuts.central.pt_diff.passes( fabs( pt1-pt2 ) ) )
        return 0.;

      //=================================================================
      //     a window in rapidity distance
      //=================================================================

      if ( !cuts_.cuts.central.rapidity_diff.passes( fabs( y1_-y2_ ) ) )
        return 0.;

      //=================================================================
      //     auxiliary quantities
      //=================================================================

      const double alpha1 = amt1/sqs_*exp( y1_ ), beta1  = amt1/sqs_*exp( -y1_ ),
                   alpha2 = amt2/sqs_*exp( y2_ ), beta2  = amt2/sqs_*exp( -y2_ );

      CG_DEBUG_LOOP( "PPtoWW:sudakov" )
        << "Sudakov parameters:\n\t"
        << "  alpha1/2 = " << alpha1 << " / " << alpha2 << "\n\t"
        << "   beta1/2 = " << beta1 << " / " << beta2 << ".";

      const double q1t2 = q1tx*q1tx+q1ty*q1ty, q2t2 = q2tx*q2tx+q2ty*q2ty;

      const double x1 = alpha1+alpha2, x2 = beta1+beta2;

      const double z1p = alpha1/x1, z1m = alpha2/x1,
                   z2p = beta1 /x2, z2m = beta2 /x2;
      CG_DEBUG_LOOP( "PPtoWW:zeta" )
        << "z(1/2)p = " << z1p << " / " << z2p << "\n\t"
        << "z(1/2)m = " << z1m << " / " << z2m << ".";

      if ( x1 > 1. || x2 > 1. )
        return 0.; // sanity check

      // FIXME FIXME FIXME
      const double ak10 = event_->getOneByRole( Particle::IncomingBeam1 ).energy(),
                   ak1z = event_->getOneByRole( Particle::IncomingBeam1 ).momentum().pz(),
                   ak20 = event_->getOneByRole( Particle::IncomingBeam2 ).energy(),
                   ak2z = event_->getOneByRole( Particle::IncomingBeam2 ).momentum().pz();
      CG_DEBUG_LOOP( "PPtoWW:incoming" )
        << "incoming particles: p1: " << ak1z << " / " << ak10 << "\n\t"
        << "                    p2: " << ak2z << " / " << ak20 << ".";

      //=================================================================
      //     additional conditions for energy-momentum conservation
      //=================================================================

      const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
      CG_DEBUG_LOOP( "PPtoWW:central" )
        << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
        << "diboson invariant mass = " << invm << " GeV";

      if ( ( cuts_.mode == KinematicsMode::ElasticInelastic
          || cuts_.mode == KinematicsMode::InelasticInelastic )
        && ( sqrt( s1_eff ) <= ( MY_+invm ) ) )
        return 0.;
      if ( ( cuts_.mode == KinematicsMode::InelasticElastic
          || cuts_.mode == KinematicsMode::InelasticInelastic )
        && ( sqrt( s2_eff ) <= ( MX_+invm ) ) )
        return 0.;

      //const double qcaptx = pcaptx, qcapty = pcapty;

      //=================================================================
      //     four-momenta of the outgoing protons (or remnants)
      //=================================================================

      const double px_plus  = ( 1.-x1 )*fabs( ak1z )*M_SQRT2,
                   px_minus = ( MX_*MX_ + q1t2 )*0.5/px_plus;

      const double py_minus = ( 1.-x2 )*fabs( ak2z )*M_SQRT2, // warning! sign of pz??
                   py_plus  = ( MY_*MY_ + q2t2 )*0.5/py_minus;

      CG_DEBUG_LOOP( "PPtoWW:pxy" )
        << "px± = " << px_plus << " / " << px_minus << "\n\t"
        << "py± = " << py_plus << " / " << py_minus << ".";

      PX_ = Particle::Momentum( -q1tx, -q1ty, ( px_plus-px_minus )*M_SQRT1_2, ( px_plus+px_minus )*M_SQRT1_2 );
      PY_ = Particle::Momentum( -q2tx, -q2ty, ( py_plus-py_minus )*M_SQRT1_2, ( py_plus+py_minus )*M_SQRT1_2 );

      CG_DEBUG_LOOP( "PPtoWW:remnants" )
        << "First remnant:  " << PX_ << ", mass = " << PX_.mass() << "\n\t"
        << "Second remnant: " << PY_ << ", mass = " << PY_.mass() << ".";

      /*assert( fabs( PX_.mass()-MX_ ) < 1.e-6 );
      assert( fabs( PY_.mass()-MY_ ) < 1.e-6 );*/

      //=================================================================
      //     four-momenta squared of the virtual photons
      //=================================================================

      const double ww = 0.5 * ( 1.+sqrt( 1.-4.*mp2_/s_ ) );

      // FIXME FIXME FIXME /////////////////////
      const Particle::Momentum q1(
        q1tx, q1ty,
        +0.5 * x1*ww*sqs_*( 1.-q1t2/x1/x1/ww/ww/s_ ),
        +0.5 * x1*ww*sqs_*( 1.+q1t2/x1/x1/ww/ww/s_ ) );
      const Particle::Momentum q2(
        q2tx, q2ty,
        -0.5 * x2*ww*sqs_*( 1.-q2t2/x2/x2/ww/ww/s_ ),
        +0.5 * x2*ww*sqs_*( 1.+q2t2/x2/x2/ww/ww/s_ ) );
      //////////////////////////////////////////

      CG_DEBUG_LOOP( "PPtoWW:partons" )
        << "First photon*:  " << q1 << ", mass2 = " << q1.mass2() << "\n\t"
        << "Second photon*: " << q2 << ", mass2 = " << q2.mass2() << ".";
      //const double q12 = q1.mass2(), q22 = q2.mass2();

      //=================================================================
      //     four-momenta of the outgoing W^+ and W^-
      //=================================================================

      p_w1_ = Particle::Momentum( pt1x, pt1y, alpha1*ak1z + beta1*ak2z, alpha1*ak10 + beta1*ak20 );
      p_w2_ = Particle::Momentum( pt2x, pt2y, alpha2*ak1z + beta2*ak2z, alpha2*ak10 + beta2*ak20 );

      CG_DEBUG_LOOP( "PPtoWW:central" )
        << "First W:  " << p_w1_ << ", mass = " << p_w1_.mass() << "\n\t"
        << "Second W: " << p_w2_ << ", mass = " << p_w2_.mass() << ".";

      //assert( fabs( p_w1_.mass()-event_->getByRole( Particle::CentralSystem )[0].mass() ) < 1.e-6 );
      //assert( fabs( p_w2_.mass()-event_->getByRole( Particle::CentralSystem )[1].mass() ) < 1.e-6 );

      //=================================================================
      //     Mendelstam variables
      //=================================================================

      //const double shat = s_*x1*x2; // ishat = 1 (approximation)
      const double shat = ( q1+q2 ).mass2(); // ishat = 2 (exact formula)

      const double that1 = ( q1-p_w1_ ).mass2(), that2 = ( q2-p_w2_ ).mass2();
      const double uhat1 = ( q1-p_w2_ ).mass2(), uhat2 = ( q2-p_w1_ ).mass2();
      CG_DEBUG_LOOP( "PPtoWW" )
        << "that(1/2) = " << that1 << " / " << that2 << "\n\t"
        << "uhat(1/2) = " << uhat1 << " / " << uhat2 << ".";

      //const double mll = sqrt( shat );

      const double that = 0.5*( that1+that2 ), uhat = 0.5*( uhat1+uhat2 );

      //=================================================================
      //     matrix elements
      //=================================================================
      double amat2 = 0.;

      //=================================================================
      //     How matrix element is calculated
      //=================================================================

      CG_DEBUG_LOOP( "PPtoWW" )
        << "matrix element mode: " << method_ << ".";

      //=================================================================
      //     matrix element for gamma gamma --> W^+ W^-
      //     (Denner+Dittmaier+Schuster)
      //     (work in collaboration with C. Royon)
      //=================================================================
      if ( method_ == 0 )
        amat2 = onShellME( shat, that, uhat );

      //=================================================================
      //     off-shell Nachtmann formulae
      //=================================================================
      else if ( method_ == 1 )
        amat2 = offShellME( shat, that, uhat, phi_qt1_+phi_qt2_, phi_qt1_-phi_qt2_ );

      if ( amat2 <= 0. )
        return 0.;

      //============================================
      //     unintegrated photon distributions
      //============================================

      const std::pair<double,double> fluxes
        = GenericKTProcess::incomingFluxes( x1, q1t2, x2, q2t2 );

      CG_DEBUG_LOOP( "PPtoWW:fluxes" )
        << "Incoming photon fluxes for (x/kt2) = "
        << "(" << x1 << "/" << q1t2 << "), "
        << "(" << x2 << "/" << q2t2 << "):\n\t"
        << fluxes.first << ", " << fluxes.second << ".";

      //=================================================================
      //     factor 2.*pi from integration over phi_sum
      //     factor 1/4 from jacobian of transformations
      //     factors 1/pi and 1/pi due to integration over
      //       d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
      //=================================================================

      const double aintegral = amat2 / ( 16.*M_PI*M_PI*( x1*x2*s_ )*( x1*x2*s_ ) )
                             * fluxes.first*M_1_PI * fluxes.second*M_1_PI * 0.25
                             * constants::GeV2toBarn;
      /*const double aintegral = amat2 / ( 16.*M_PI*M_PI*x1*x1*x2*x2*s_*s_ )
                             * fluxes.first*M_1_PI * fluxes.second*M_1_PI
                             * constants::GeV2toBarn * 0.25;*/

      //=================================================================
      return aintegral*qt1_*qt2_*pt_diff_;
      //=================================================================
    }

    void
    PPtoWW::fillCentralParticlesKinematics()
    {
      // randomise the charge of the outgoing leptons
      short sign = ( drand() > 0.5 ) ? +1 : -1;

      //=================================================================
      //     first outgoing lepton
      //=================================================================
      Particle& ow1 = event_->getByRole( Particle::CentralSystem )[0];
      ow1.setPdgId( ow1.pdgId(), sign );
      ow1.setStatus( Particle::Status::Undecayed );
      ow1.setMomentum( p_w1_ );

      //=================================================================
      //     second outgoing lepton
      //=================================================================
      Particle& ow2 = event_->getByRole( Particle::CentralSystem )[1];
      ow2.setPdgId( ow2.pdgId(), -sign );
      ow2.setStatus( Particle::Status::Undecayed );
      ow2.setMomentum( p_w2_ );
    }

    double
    PPtoWW::onShellME( double shat, double that, double uhat )
    {
      const double mw4 = mw2_*mw2_;

      const double term1 = 2.*shat * ( 2.*shat+3.*mw2_ ) / ( 3.*( mw2_-that )*( mw2_-uhat ) );
      const double term2 = 2.*shat*shat * ( shat*shat + 3.*mw4 ) / ( 3.*pow( mw2_-that, 2 )*pow( mw2_-uhat, 2 ) );

      const double auxil_gamgam = 1.-term1+term2;
      const double beta = sqrt( 1.-4.*mw2_/shat );

      return 3.*constants::alphaEM*constants::alphaEM*beta / ( 2.*shat ) * auxil_gamgam / ( beta/( 64.*M_PI*M_PI*shat ) );
    }

    double
    PPtoWW::offShellME( double shat, double that, double uhat, double phi_sum, double phi_diff )
    {
      const double e2 = 4.*M_PI*constants::alphaEM;

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
      return e2*e2*( amat2_0+amat2_1+amat2_interf );
    }

    double
    PPtoWW::amplitudeWW( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 )
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
      // longitudinal-longitudinal
      if ( lam3 == 0 && lam4 == 0 )
        return invA*inv_gamma2*( ( gamma2+1. )*( 1.-lam1*lam2 )*sin_theta2 - ( 1.+lam1*lam2 ) );
      // transverse-longitudinal
      if ( lam4 == 0 )
        return invA*( -M_SQRT2*inv_gamma*( lam1-lam2 )*( 1.+lam1*lam3*cos_theta )*sin_theta );
      // longitudinal-transverse
      if ( lam3 == 0 )
        return invA*( -M_SQRT2*inv_gamma*( lam2-lam1 )*( 1.+lam2*lam4*cos_theta )*sin_theta );
      // transverse-transverse
      if ( lam3 != 0 && lam4 != 0 )
        return -0.5*invA*( 2.*beta*( lam1+lam2 )*( lam3+lam4 )
                          -inv_gamma2*( 1.+lam3*lam4 )*( 2.*lam1*lam2+( 1.-lam1*lam2 ) * cos_theta2 )
                          +( 1.+lam1*lam2*lam3*lam4 )*( 3.+lam1*lam2 )
                          +2.*( lam1-lam2 )*( lam3-lam4 )*cos_theta
                          +( 1.-lam1*lam2 )*( 1.-lam3*lam4 )*cos_theta2 );
      return 0.;
    }
    // register process and define aliases
    REGISTER_PROCESS( pptoww, PPtoWW )
  }
}
