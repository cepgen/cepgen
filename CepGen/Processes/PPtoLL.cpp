#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"
#include <assert.h>

namespace CepGen
{
  namespace Process
  {
    PPtoLL::PPtoLL() :
      GenericKTProcess( "pptoll", "ɣɣ → l⁺l¯", { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } ),
      y1_( 0. ), y2_( 0. ), pt_diff_( 0. ), phi_pt_diff_( 0. )
    {}

    void
    PPtoLL::preparePhaseSpace()
    {
      registerVariable( y1_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "First outgoing lepton rapidity" );
      registerVariable( y2_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "Second outgoing lepton rapidity" );
      registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::pt_diff], { 0., 50. }, "Leptons transverse momentum difference" );
      registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::phi_pt_diff], { 0., 2.*M_PI }, "Leptons azimuthal angle difference" );
    }

    double
    PPtoLL::computeKTFactorisedMatrixElement()
    {
      const double ml = event_->getByRole( Particle::CentralSystem )[0].mass(), ml2 = ml*ml;

      const unsigned int iterm11 = 1, // Long-long
                         iterm22 = 1, // Trans-trans
                         iterm12 = 1, // Long-trans
                         itermtt = 1; // Trans-trans(')

      //=================================================================
      //     How matrix element is calculated
      //=================================================================

      const bool off_shell = true;

      //=================================================================
      //     two terms in Wolfgang's formula for
      //     off-shell gamma gamma --> l^+ l^-
      //=================================================================

      const unsigned int imat1 = 2, imat2 = 0;

      //=================================================================
      //     matrix element computation
      //=================================================================
      //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2_, 2))/s_*s_));

      // Inner photons
      const double q1tx = qt1_*cos( phi_qt1_ ), q1ty = qt1_*sin( phi_qt1_ ),
                   q2tx = qt2_*cos( phi_qt2_ ), q2ty = qt2_*sin( phi_qt2_ );
      CG_DEBUG_LOOP( "PPtoLL" )
        << "q1t(x/y) = " << q1tx << " / " << q1ty << "\n\t"
        << "q2t(x/y) = " << q2tx << " / " << q2ty;

      // Two-photon system
      const double ptsumx = q1tx+q2tx,
                   ptsumy = q1ty+q2ty,
                   ptsum = sqrt( ptsumx*ptsumx+ptsumy*ptsumy );

      const double ptdiffx = pt_diff_*cos( phi_pt_diff_ ),
                   ptdiffy = pt_diff_*sin( phi_pt_diff_ );

      // Outgoing leptons
      const double pt1x = ( ptsumx+ptdiffx )*0.5, pt1y = ( ptsumy+ptdiffy )*0.5, pt1 = std::hypot( pt1x, pt1y ),
                   pt2x = ( ptsumx-ptdiffx )*0.5, pt2y = ( ptsumy-ptdiffy )*0.5, pt2 = std::hypot( pt2x, pt2y );

      const Limits pt_limits = cuts_.cuts.central[Cuts::pt_single];
      if ( !pt_limits.passes( pt1 ) || !pt_limits.passes( pt2 ) )
        return 0.;

      // transverse mass for the two leptons
      const double amt1 = sqrt( pt1*pt1+ml2 ),
                   amt2 = sqrt( pt2*pt2+ml2 );

      //=================================================================
      //     a window in transverse momentum difference
      //=================================================================

      if ( cuts_.cuts.central.count( Cuts::pt_diff ) > 0
        && !cuts_.cuts.central.at( Cuts::pt_diff ).passes( fabs( pt1-pt2 ) ) )
        return 0.;

      //=================================================================
      //     a window in rapidity distance
      //=================================================================

      if ( cuts_.cuts.central.count( Cuts::rapidity_diff ) > 0
        && !cuts_.cuts.central[Cuts::rapidity_diff].passes( fabs( y1_-y2_ ) ) )
        return 0.;

      //=================================================================
      //     auxiliary quantities
      //=================================================================

      const double alpha1 = amt1/sqs_*exp( y1_ ), beta1  = amt1/sqs_*exp( -y1_ ),
                   alpha2 = amt2/sqs_*exp( y2_ ), beta2  = amt2/sqs_*exp( -y2_ );

      CG_DEBUG_LOOP( "PPtoLL" )
        << "Sudakov parameters:\n\t"
        << "  alpha1/2 = " << alpha1 << " / " << alpha2 << "\n\t"
        << "   beta1/2 = " << beta1 << " / " << beta2 << ".";

      const double q1t2 = q1tx*q1tx+q1ty*q1ty, q2t2 = q2tx*q2tx+q2ty*q2ty;

      const double x1 = alpha1+alpha2, x2 = beta1+beta2;

      const double z1p = alpha1/x1, z1m = alpha2/x1,
                   z2p = beta1 /x2, z2m = beta2 /x2;
      CG_DEBUG_LOOP( "PPtoLL" )
        << "z(1/2)p = " << z1p << " / " << z2p << "\n\t"
        << "z(1/2)m = " << z1m << " / " << z2m << ".";

      if ( x1 > 1. || x2 > 1. )
        return 0.; // sanity check

      // FIXME FIXME FIXME
      const double ak10 = event_->getOneByRole( Particle::IncomingBeam1 ).energy(),
                   ak1z = event_->getOneByRole( Particle::IncomingBeam1 ).momentum().pz(),
                   ak20 = event_->getOneByRole( Particle::IncomingBeam2 ).energy(),
                   ak2z = event_->getOneByRole( Particle::IncomingBeam2 ).momentum().pz();
      CG_DEBUG_LOOP( "PPtoLL" )
        << "incoming particles: p1: " << ak1z << " / " << ak10 << "\n\t"
        << "                    p2: " << ak2z << " / " << ak20;

      //=================================================================
      //     additional conditions for energy-momentum conservation
      //=================================================================

      const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
      const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh(y1_-y2_) - ptsum*ptsum );
      CG_DEBUG_LOOP( "PPtoLL" )
        << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
        << "dilepton invariant mass = " << invm << " GeV.";

      if ( ( cuts_.mode == Kinematics::Mode::ElasticInelastic
          || cuts_.mode == Kinematics::Mode::InelasticInelastic )
        && ( sqrt( s1_eff ) <= ( MY_+invm ) ) )
        return 0.;
      if ( ( cuts_.mode == Kinematics::Mode::InelasticElastic
          || cuts_.mode == Kinematics::Mode::InelasticInelastic )
        && ( sqrt( s2_eff ) <= ( MX_+invm ) ) )
        return 0.;

      //const double qcaptx = pcaptx, qcapty = pcapty;

      //=================================================================
      //     four-momenta of the outgoing protons (or remnants)
      //=================================================================

      const double px_plus  = ( 1.-x1 )*fabs( ak1z )*M_SQRT2,
                   px_minus = ( MX_*MX_ + q1tx*q1tx + q1ty*q1ty )*0.5/px_plus;

      const double py_minus = ( 1.-x2 )*fabs( ak2z )*M_SQRT2, // warning! sign of pz??
                   py_plus  = ( MY_*MY_ + q2tx*q2tx + q2ty*q2ty )*0.5/py_minus;

      CG_DEBUG_LOOP( "PPtoLL" )
        << "px± = " << px_plus << " / " << px_minus << "\n\t"
        << "py± = " << py_plus << " / " << py_minus << ".";

      PX_ = Particle::Momentum( -q1tx, -q1ty, ( px_plus-px_minus )*M_SQRT1_2, ( px_plus+px_minus )*M_SQRT1_2 );
      PY_ = Particle::Momentum( -q2tx, -q2ty, ( py_plus-py_minus )*M_SQRT1_2, ( py_plus+py_minus )*M_SQRT1_2 );

      CG_DEBUG_LOOP( "PPtoLL" )
        << "First remnant:  " << PX_ << ", mass = " << PX_.mass() << "\n\t"
        << "Second remnant: " << PY_ << ", mass = " << PY_.mass() << ".";

      assert( fabs( PX_.mass()-MX_ ) < 1.e-6 );
      assert( fabs( PY_.mass()-MY_ ) < 1.e-6 );

      //=================================================================
      //     four-momenta of the outgoing l^+ and l^-
      //=================================================================

      const Particle::Momentum p1( pt1x, pt1y, alpha1*ak1z + beta1*ak2z, alpha1*ak10 + beta1*ak20 );
      const Particle::Momentum p2( pt2x, pt2y, alpha2*ak1z + beta2*ak2z, alpha2*ak10 + beta2*ak20 );
      CG_DEBUG_LOOP( "PPtoLL" )
        << "unboosted first lepton:  " << p1 << ", mass = " << p1.mass() << "\n\t"
        << "          second lepton: " << p2 << ", mass = " << p2.mass() << ".";

      Pl1_ = Particle::Momentum( pt1x, pt1y, sqrt( pt1*pt1 + ml2 )*sinh( y1_ ), sqrt( pt1*pt1 + ml2 )*cosh( y1_ ) );
      Pl2_ = Particle::Momentum( pt2x, pt2y, sqrt( pt2*pt2 + ml2 )*sinh( y2_ ), sqrt( pt2*pt2 + ml2 )*cosh( y2_ ) );

      CG_DEBUG_LOOP( "PPtoLL" )
        << "First lepton:  " << Pl1_ << ", mass = " << Pl1_.mass() << "\n\t"
        << "Second lepton: " << Pl2_ << ", mass = " << Pl2_.mass() << ".";

      assert( fabs( Pl1_.mass()-event_->getByRole( Particle::CentralSystem )[0].mass() ) < 1.e-6 );
      assert( fabs( Pl2_.mass()-event_->getByRole( Particle::CentralSystem )[1].mass() ) < 1.e-6 );

      //=================================================================
      //     four-momenta squared of the virtual photons
      //=================================================================

      // FIXME FIXME FIXME /////////////////////
      const Particle::Momentum q1( q1tx, q1ty, 0., 0. );
      const Particle::Momentum q2( q2tx, q2ty, 0., 0. );
      //////////////////////////////////////////

      CG_DEBUG_LOOP( "PPtoLL" )
        << "First photon*:  " << q1 << ", mass2 = " << q1.mass2() << "\n\t"
        << "Second photon*: " << q2 << ", mass2 = " << q2.mass2() << ".";

      //=================================================================
      //     Mendelstam variables
      //=================================================================

      //const double shat = s_*x1*x2; // ishat = 1 (approximation)
      //const double shat = ( q1+q2 ).mass2(); // ishat = 2 (exact formula)

      const double that1 = ( q1-p1 ).mass2(), that2 = ( q2-p2 ).mass2();
      const double uhat1 = ( q1-p2 ).mass2(), uhat2 = ( q2-p1 ).mass2();
      CG_DEBUG_LOOP( "PPtoLL" )
        << "that(1/2) = " << that1 << " / " << that2 << "\n\t"
        << "uhat(1/2) = " << uhat1 << " / " << uhat2 << ".";

      //const double mll = sqrt( shat );

      const double that = 0.5*( that1+that2 ), uhat = 0.5*( uhat1+uhat2 );

      //=================================================================
      //     matrix elements
      //=================================================================
      double amat2 = 0.;
      if ( !off_shell ) {

        //=================================================================
        //     on-shell formula for M^2
        //=================================================================
        const double ml4 = ml2*ml2, ml8 = ml4*ml4;

        const double term1  =  6. *ml8,
                     term2  = -3. *ml4*that*that,
                     term3  = -14.*ml4*that*uhat,
                     term4  = -3. *ml4*uhat*uhat,
                     term5  =      ml2*that*that*that,
                     term6  =  7.* ml2*that*that*uhat,
                     term7  =  7.* ml2*that*uhat*uhat,
                     term8  =      ml2*uhat*uhat*uhat,
                     term9  =         -that*that*that*uhat,
                     term10 =         -that*uhat*uhat*uhat;

        const double auxil_gamgam = -2.*( term1+term2+term3+term4+term5+term6+term7+term8+term9+term10 )/( pow( ( ml2-that )*( ml2-uhat ), 2) );
        const double g_em_sq = 4.*M_PI*Constants::alphaEM;
        amat2 = g_em_sq*g_em_sq*auxil_gamgam;
      }
      else if ( off_shell ) {

        //=================================================================
        //     Wolfgang's formulae
        //=================================================================

        const double ak1_x = z1m*pt1x - z1p*pt2x, ak1_y = z1m*pt1y - z1p*pt2y,
                     ak2_x = z2m*pt1x - z2p*pt2x, ak2_y = z2m*pt1y - z2p*pt2y;

        const double t1abs = ( q1t2 + x1*( MX_*MX_-mp2_ )+x1*x1*mp2_ )/( 1.-x1 ),
                     t2abs = ( q2t2 + x2*( MY_*MY_-mp2_ )+x2*x2*mp2_ )/( 1.-x2 );

        const double eps12 = ml2 + z1p*z1m*t1abs,
                     eps22 = ml2 + z2p*z2m*t2abs;

        const double Phi10 = 1./( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                            -1./( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi11_x = ( ak1_x + z1p*q2tx )/( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                              -( ak1_x - z1m*q2tx )/( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi11_y = ( ak1_y + z1p*q2ty )/( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                              -( ak1_y - z1m*q2ty )/( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi102 = Phi10*Phi10;

        const double Phi20 = 1./( pow( ak2_x + z2p*q1tx, 2 )+pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                            -1./( pow( ak2_x - z2m*q1tx, 2 )+pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi21_x = ( ak2_x + z2p*q1tx )/( pow( ak2_x + z2p*q1tx, 2 ) + pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                              -( ak2_x - z2m*q1tx )/( pow( ak2_x - z2m*q1tx ,2 ) + pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi21_y = ( ak2_y + z2p*q1ty )/( pow( ak2_x + z2p*q1tx, 2 ) + pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                              -( ak2_y - z2m*q1ty )/( pow( ak2_x - z2m*q1tx, 2 ) + pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi202 = Phi20*Phi20;

        const double Phi11_dot_e = ( Phi11_x*q1tx + Phi11_y*q1ty )/qt1_, Phi11_cross_e = ( Phi11_x*q1ty-Phi11_y*q1tx )/qt1_;
        const double Phi21_dot_e = ( Phi21_x*q2tx + Phi21_y*q2ty )/qt2_, Phi21_cross_e = ( Phi21_x*q2ty-Phi21_y*q2tx )/qt2_;
        CG_DEBUG_LOOP( "PPtoLL" )
          << "Phi1: E, px, py = " << Phi10 << ", " << Phi11_x << ", " << Phi11_y << "\n\t"
          << "Phi2: E, px, py = " << Phi20 << ", " << Phi21_x << ", " << Phi21_y << "\n\t"
          << "(dot):   " << Phi11_dot_e << " / " << Phi21_dot_e << "\n\t"
          << "(cross): " << Phi11_cross_e << " / " << Phi21_cross_e << ".";

        const double aux2_1 = iterm11 * ( ml2 + 4.*z1p*z1p*z1m*z1m*t1abs ) * Phi102
                             +iterm22 * ( ( z1p*z1p + z1m*z1m )*( Phi11_dot_e*Phi11_dot_e + Phi11_cross_e*Phi11_cross_e ) )
                             +itermtt * ( Phi11_cross_e*Phi11_cross_e - Phi11_dot_e*Phi11_dot_e )
                             -iterm12 * 4.*z1p*z1m*( z1p-z1m ) * Phi10 * ( q1tx*Phi11_x + q1ty*Phi11_y );

        const double aux2_2 = iterm11 * ( ml2 + 4.*z2p*z2p*z2p*z2m*t2abs ) * Phi202
                             +iterm22 * ( ( z2p*z2p + z2m*z2m )*( Phi21_dot_e*Phi21_dot_e + Phi21_cross_e*Phi21_cross_e ) )
                             +itermtt * ( Phi21_cross_e*Phi21_cross_e - Phi21_dot_e*Phi21_dot_e )
                             -iterm12 * 4.*z2p*z2m*( z2p-z2m ) * Phi20 * ( q2tx*Phi21_x + q2ty*Phi21_y );

        //=================================================================
        //     convention of matrix element as in our kt-factorization
        //     for heavy flavours
        //=================================================================

        const double norm = 16.*M_PI*M_PI*Constants::alphaEM*Constants::alphaEM * pow( x1*x2*s_, 2 );

        double amat2_1 = norm * aux2_1*2.*z1p*z1m*t1abs/( q1t2*q2t2 ) * t2abs/q2t2;
        double amat2_2 = norm * aux2_2*2.*z2p*z2m*t2abs/( q1t2*q2t2 );

        //=================================================================
        //     symmetrization
        //=================================================================

        amat2 = 0.5*( imat1*amat2_1 + imat2*amat2_2 );

        CG_DEBUG_LOOP( "PPtoLL" )
          << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
          << "amat2(1/2), amat2 = " << amat2_1 << " / " << amat2_2 << " / " << amat2 << ".";
      }

      //============================================
      //     unintegrated photon distributions
      //============================================

      GenericKTProcess::computeIncomingFluxes( x1, q1t2, x2, q2t2 );

      //=================================================================
      //     factor 2.*pi from integration over phi_sum
      //     factor 1/4 from jacobian of transformations
      //     factors 1/pi and 1/pi due to integration over
      //       d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
      //=================================================================

      const double aintegral = amat2 / ( 16.*M_PI*M_PI*( x1*x2*s_ )*( x1*x2*s_ ) )
                             * flux1_/M_PI * flux2_/M_PI * 0.25
                             * Constants::GeV2toBarn;

      //=================================================================
      return aintegral*qt1_*qt2_*pt_diff_;
      //=================================================================
    }

    void
    PPtoLL::fillCentralParticlesKinematics()
    {
      // randomise the charge of the outgoing leptons
      short sign = ( drand() > 0.5 ) ? +1 : -1;

      //=================================================================
      //     first outgoing lepton
      //=================================================================
      Particle& ol1 = event_->getByRole( Particle::CentralSystem )[0];
      ol1.setPdgId( ol1.pdgId(), sign );
      ol1.setStatus( Particle::FinalState );
      ol1.setMomentum( Pl1_ );

      //=================================================================
      //     second outgoing lepton
      //=================================================================
      Particle& ol2 = event_->getByRole( Particle::CentralSystem )[1];
      ol2.setPdgId( ol2.pdgId(), -sign );
      ol2.setStatus( Particle::FinalState );
      ol2.setMomentum( Pl2_ );
    }
  }
}
