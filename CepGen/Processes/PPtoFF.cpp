#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"
#include <assert.h>

namespace CepGen
{
  namespace Process
  {
    PPtoFF::PPtoFF() :
      GenericKTProcess( "pptoff", "ɣɣ → f⁺f¯", { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } ),
      y1_( 0. ), y2_( 0. ), pt_diff_( 0. ), phi_pt_diff_( 0. )
    {}

    void
    PPtoFF::preparePhaseSpace()
    {
      registerVariable( y1_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "First outgoing fermion rapidity" );
      registerVariable( y2_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "Second outgoing fermion rapidity" );
      registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central.pt_diff, { 0., 50. }, "Fermions transverse momentum difference" );
      registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central.phi_pt_diff, { 0., 2.*M_PI }, "Fermions azimuthal angle difference" );

      const PDG& pdg_f = cuts_.central_system[0];
      mf_ = ParticleProperties::mass( pdg_f );
      mf2_ = mf_*mf_;
      qf_ = ParticleProperties::charge( pdg_f );
      colf_ = ParticleProperties::colours( pdg_f );
      CG_DEBUG( "PPtoFF:prepare" )
        << "Produced system (" << pdg_f << pdg_f << ") will have mass = " << mf_ << " GeV, charge = " << qf_ << " e";
    }

    double
    PPtoFF::computeKTFactorisedMatrixElement()
    {
      //=================================================================
      //     How matrix element is calculated
      //=================================================================

      const unsigned short method = 1;

      //=================================================================
      //     matrix element computation
      //=================================================================
      //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2_, 2))/s_*s_));

      // Inner photons
      const double q1tx = qt1_*cos( phi_qt1_ ), q1ty = qt1_*sin( phi_qt1_ ),
                   q2tx = qt2_*cos( phi_qt2_ ), q2ty = qt2_*sin( phi_qt2_ );
      CG_DEBUG_LOOP( "PPtoFF" )
        << "q1t(x/y) = " << q1tx << " / " << q1ty << "\n\t"
        << "q2t(x/y) = " << q2tx << " / " << q2ty;

      // Two-photon system
      const double ptsumx = q1tx+q2tx, ptsumy = q1ty+q2ty,
                   ptsum = std::hypot( ptsumx, ptsumy );

      const double ptdiffx = pt_diff_*cos( phi_pt_diff_ ),
                   ptdiffy = pt_diff_*sin( phi_pt_diff_ );

      // Outgoing fermions
      const double pt1x = ( ptsumx+ptdiffx )*0.5, pt1y = ( ptsumy+ptdiffy )*0.5, pt1 = std::hypot( pt1x, pt1y ),
                   pt2x = ( ptsumx-ptdiffx )*0.5, pt2y = ( ptsumy-ptdiffy )*0.5, pt2 = std::hypot( pt2x, pt2y );

      const Limits& pt_limits = cuts_.cuts.central.pt_single;
      if ( !pt_limits.passes( pt1 ) || !pt_limits.passes( pt2 ) )
        return 0.;

      // transverse mass for the two fermions
      const double amt1 = std::hypot( pt1, mf_ ), amt2 = std::hypot( pt2, mf_ );

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

      CG_DEBUG_LOOP( "PPtoFF" )
        << "Sudakov parameters:\n\t"
        << "  alpha1/2 = " << alpha1 << " / " << alpha2 << "\n\t"
        << "   beta1/2 = " << beta1 << " / " << beta2 << ".";

      const double q1t2 = q1tx*q1tx+q1ty*q1ty, q2t2 = q2tx*q2tx+q2ty*q2ty;

      const double x1 = alpha1+alpha2, x2 = beta1+beta2;

      const double z1p = alpha1/x1, z1m = alpha2/x1,
                   z2p = beta1 /x2, z2m = beta2 /x2;
      CG_DEBUG_LOOP( "PPtoFF" )
        << "z(1/2)p = " << z1p << " / " << z2p << "\n\t"
        << "z(1/2)m = " << z1m << " / " << z2m << ".";

      if ( x1 > 1. || x2 > 1. )
        return 0.; // sanity check

      // FIXME FIXME FIXME
      const double ak10 = event_->getOneByRole( Particle::IncomingBeam1 ).energy(),
                   ak1z = event_->getOneByRole( Particle::IncomingBeam1 ).momentum().pz(),
                   ak20 = event_->getOneByRole( Particle::IncomingBeam2 ).energy(),
                   ak2z = event_->getOneByRole( Particle::IncomingBeam2 ).momentum().pz();
      CG_DEBUG_LOOP( "PPtoFF" )
        << "incoming particles: p1: " << ak1z << " / " << ak10 << "\n\t"
        << "                    p2: " << ak2z << " / " << ak20;

      //=================================================================
      //     additional conditions for energy-momentum conservation
      //=================================================================

      const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
      const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh(y1_-y2_) - ptsum*ptsum );
      CG_DEBUG_LOOP( "PPtoFF" )
        << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
        << "central system's invariant mass = " << invm << " GeV.";

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

      CG_DEBUG_LOOP( "PPtoFF" )
        << "px± = " << px_plus << " / " << px_minus << "\n\t"
        << "py± = " << py_plus << " / " << py_minus << ".";

      PX_ = Particle::Momentum( -q1tx, -q1ty, ( px_plus-px_minus )*M_SQRT1_2, ( px_plus+px_minus )*M_SQRT1_2 );
      PY_ = Particle::Momentum( -q2tx, -q2ty, ( py_plus-py_minus )*M_SQRT1_2, ( py_plus+py_minus )*M_SQRT1_2 );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "First remnant:  " << PX_ << ", mass = " << PX_.mass() << "\n\t"
        << "Second remnant: " << PY_ << ", mass = " << PY_.mass() << ".";

      assert( fabs( PX_.mass()-MX_ ) < 1.e-6 );
      assert( fabs( PY_.mass()-MY_ ) < 1.e-6 );

      //=================================================================
      //     four-momenta of the outgoing l^+ and l^-
      //=================================================================

      const Particle::Momentum p1( pt1x, pt1y, alpha1*ak1z + beta1*ak2z, alpha1*ak10 + beta1*ak20 );
      const Particle::Momentum p2( pt2x, pt2y, alpha2*ak1z + beta2*ak2z, alpha2*ak10 + beta2*ak20 );
      CG_DEBUG_LOOP( "PPtoFF" )
        << "unboosted first fermion:  " << p1 << ", mass = " << p1.mass() << "\n\t"
        << "          second fermion: " << p2 << ", mass = " << p2.mass() << ".";

      p_f1_ = Particle::Momentum( pt1x, pt1y, sqrt( pt1*pt1 + mf2_ )*sinh( y1_ ), sqrt( pt1*pt1 + mf2_ )*cosh( y1_ ) );
      p_f2_ = Particle::Momentum( pt2x, pt2y, sqrt( pt2*pt2 + mf2_ )*sinh( y2_ ), sqrt( pt2*pt2 + mf2_ )*cosh( y2_ ) );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "First fermion:  " << p_f1_ << ", mass = " << p_f1_.mass() << "\n\t"
        << "Second fermion: " << p_f2_ << ", mass = " << p_f2_.mass() << ".";

      assert( fabs( p_f1_.mass()-mf_ ) < 1.e-6 );
      assert( fabs( p_f2_.mass()-mf_ ) < 1.e-6 );

      //=================================================================
      //     four-momenta squared of the virtual photons
      //=================================================================

      // FIXME FIXME FIXME /////////////////////
      const Particle::Momentum q1( q1tx, q1ty, 0., 0. );
      const Particle::Momentum q2( q2tx, q2ty, 0., 0. );
      //////////////////////////////////////////

      CG_DEBUG_LOOP( "PPtoFF" )
        << "First photon*:  " << q1 << ", mass2 = " << q1.mass2() << "\n\t"
        << "Second photon*: " << q2 << ", mass2 = " << q2.mass2() << ".";

      //=================================================================
      //     Mendelstam variables
      //=================================================================

      //const double shat = s_*x1*x2; // ishat = 1 (approximation)
      const double shat = ( q1+q2 ).mass2(); // ishat = 2 (exact formula)

      const double that1 = ( q1-p1 ).mass2(), that2 = ( q2-p2 ).mass2();
      const double uhat1 = ( q1-p2 ).mass2(), uhat2 = ( q2-p1 ).mass2();
      CG_DEBUG_LOOP( "PPtoFF" )
        << "that(1/2) = " << that1 << " / " << that2 << "\n\t"
        << "uhat(1/2) = " << uhat1 << " / " << uhat2 << ".";

      //const double mll = sqrt( shat );

      const double that = 0.5*( that1+that2 ), uhat = 0.5*( uhat1+uhat2 );

      //=================================================================
      //     matrix elements
      //=================================================================
      double amat2 = 0.;

      CG_DEBUG_LOOP( "PPtoFF" )
        << "matrix element mode: " << method << ".";

      if ( method == 0 )
        amat2 = onShellME( shat, that, uhat );

      else if ( method == 1 ) {
        Particle::Momentum ak1 = ( z1m*p_f1_-z1p*p_f2_ ), ak2 = ( z2m*p_f1_-z2p*p_f2_ );
        const double t1abs = ( q1.pt2() + x1*( MX_*MX_-mp2_ )+x1*x1*mp2_ )/( 1.-x1 ),
                     t2abs = ( q2.pt2() + x2*( MY_*MY_-mp2_ )+x2*x2*mp2_ )/( 1.-x2 );
        amat2 = offShellME( t1abs, t2abs, z1m, z1p, z2m, z2p, q1, q2, ak1, ak2 ) * pow( x1*x2*s_, 2 );
      }

      const double g_em = 4.*M_PI*Constants::alphaEM*qf_;

      //============================================
      //     unintegrated photon distributions
      //============================================

      GenericKTProcess::computeIncomingFluxes( x1, q1t2, x2, q2t2 );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "Incoming photon fluxes for (x/kt2) = "
        << "(" << x1 << "/" << q1t2 << "), "
        << "(" << x2 << "/" << q2t2 << "):\n\t"
        << flux1_ << ", " << flux2_ << ".";

      //=================================================================
      //     factor 2.*pi from integration over phi_sum
      //     factor 1/4 from jacobian of transformations
      //     factors 1/pi and 1/pi due to integration over
      //       d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
      //=================================================================

      const double aintegral = colf_ * g_em*g_em * amat2
                             * 1. / pow( 4.*M_PI*( x1*x2*s_ ), 2 )
                             * flux1_/M_PI * flux2_/M_PI * 0.25
                             * Constants::GeV2toBarn;

      //=================================================================
      return aintegral*qt1_*qt2_*pt_diff_;
      //=================================================================
    }

    void
    PPtoFF::fillCentralParticlesKinematics()
    {
      // randomise the charge of the outgoing fermions
      short sign = ( drand() > 0.5 ) ? +1 : -1;

      //=================================================================
      //     first outgoing fermion
      //=================================================================
      Particle& of1 = event_->getByRole( Particle::CentralSystem )[0];
      of1.setPdgId( of1.pdgId(), sign );
      of1.setStatus( Particle::FinalState );
      of1.setMomentum( p_f1_ );

      //=================================================================
      //     second outgoing fermion
      //=================================================================
      Particle& of2 = event_->getByRole( Particle::CentralSystem )[1];
      of2.setPdgId( of2.pdgId(), -sign );
      of2.setStatus( Particle::FinalState );
      of2.setMomentum( p_f2_ );
    }

    double
    PPtoFF::onShellME( double shat, double that, double uhat ) const
    {
      //=================================================================
      //     on-shell formula for M^2
      //=================================================================
      const double ml4 = mf2_*mf2_, ml8 = ml4*ml4;

      const double term1  =  6. *ml8,
                   term2  = -3. *ml4 *that*that,
                   term3  = -14.*ml4 *that*uhat,
                   term4  = -3. *ml4 *uhat*uhat,
                   term5  =      mf2_*that*that*that,
                   term6  =  7.* mf2_*that*that*uhat,
                   term7  =  7.* mf2_*that*uhat*uhat,
                   term8  =      mf2_*uhat*uhat*uhat,
                   term9  =          -that*that*that*uhat,
                   term10 =          -that*uhat*uhat*uhat;

      return -2.*( term1+term2+term3+term4+term5+term6+term7+term8+term9+term10 )/( pow( ( mf2_-that )*( mf2_-uhat ), 2) );
    }

    double
    PPtoFF::offShellME( double t1abs, double t2abs, double z1m, double z1p, double z2m, double z2p, const Particle::Momentum& q1, const Particle::Momentum& q2, const Particle::Momentum& ak1, const Particle::Momentum& ak2 ) const
    {
      //=================================================================
      //     Wolfgang's formulae
      //=================================================================

      const double z1 = z1p*z1m, z2 = z2p*z2m;
      const double eps12 = mf2_ + z1*t1abs,
                   eps22 = mf2_ + z2*t2abs;

      const Particle::Momentum phi1(
        ( ak1+z1p*q2 ).px()/( ( ak1+z1p*q2 ).pt2() + eps12 )-( ak1-z1m*q2 ).px()/( ( ak1-z1m*q2 ).pt2() + eps12 ),
        ( ak1+z1p*q2 ).py()/( ( ak1+z1p*q2 ).pt2() + eps12 )-( ak1-z1m*q2 ).py()/( ( ak1-z1m*q2 ).pt2() + eps12 ),
        0.,
        1./( ( ak1+z1p*q2 ).pt2() + eps12 )-1./( ( ak1-z1m*q2 ).pt2() + eps12 )
      );

      const Particle::Momentum phi2(
        ( ak2+z2p*q1 ).px()/( ( ak2+z2p*q1 ).pt2() + eps22 )-( ak2-z2m*q1 ).px()/( ( ak2-z2m*q1 ).pt2() + eps22 ),
        ( ak2+z2p*q1 ).py()/( ( ak2+z2p*q1 ).pt2() + eps22 )-( ak2-z2m*q1 ).py()/( ( ak2-z2m*q1 ).pt2() + eps22 ),
        0.,
        1./( ( ak2+z2p*q1 ).pt2() + eps22 )-1./( ( ak2-z2m*q1 ).pt2() + eps22 )
      );

      const double dot1 = phi1.threeProduct( q1 )/qt1_, cross1 = phi1.crossProduct( q1 )/qt1_;
      const double dot2 = phi2.threeProduct( q2 )/qt2_, cross2 = phi2.crossProduct( q2 )/qt2_;
      CG_DEBUG_LOOP( "PPtoFF" )
        << "phi1 = " << phi1 << "\n\t"
        << "phi2 = " << phi2 << "\n\t"
        << "(dot):   " << dot1 << " / " << dot2 << "\n\t"
        << "(cross): " << cross1 << " / " << cross2 << ".";

      //=================================================================
      //     six terms in Wolfgang's formula for
      //     off-shell gamma gamma --> l^+ l^-
      //=================================================================

      const unsigned short imat1 = 1, imat2 = 1;
      const unsigned short itermLL = 1, itermTT = 1, itermLT = 1, itermtt = 1;

      const double aux2_1 = itermLL * ( mf2_ + 4.*z1*z1*t1abs ) * phi1.energy2()
                           +itermTT * ( ( z1p*z1p + z1m*z1m )*( dot1*dot1 + cross1*cross1 ) )
                           +itermtt * ( cross1*cross1 - dot1*dot1 )
                           -itermLT * 4.*z1*( z1p-z1m ) * phi1.energy() * q1.threeProduct( phi1 );

      const double aux2_2 = itermLL * ( mf2_ + 4.*z2*z2*t2abs ) * phi2.energy2()
                           +itermTT * ( ( z2p*z2p + z2m*z2m )*( dot2*dot2 + cross2*cross2 ) )
                           +itermtt * ( cross2*cross2 - dot2*dot2 )
                           -itermLT * 4.*z2*( z2p-z2m ) * phi2.energy() * q2.threeProduct( phi2 );

      //=================================================================
      //     convention of matrix element as in our kt-factorization
      //     for heavy flavours
      //=================================================================

      const double amat2_1 = aux2_1*2.*z1p*z1m*t1abs/( q1.pt2()*q2.pt2() ) * t2abs/q2.pt2(),
                   amat2_2 = aux2_2*2.*z2p*z2m*t2abs/( q1.pt2()*q2.pt2() );

      //=================================================================
      //     symmetrization
      //=================================================================

      CG_DEBUG_LOOP( "PPtoFF" )
        << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
        << "amat2(1/2), amat2 = " << amat2_1 << " / " << amat2_2 << " / " << ( 0.5*( imat1*amat2_1 + imat2*amat2_2 ) ) << ".";

      return 0.5*( imat1*amat2_1 + imat2*amat2_2 );
    }
  }
}
