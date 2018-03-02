#include "PPtoLL.h"
#include "CepGen/Core/Exception.h"
#include <assert.h>

namespace CepGen
{
  namespace Process
  {
    PPtoLL::PPtoLL() :
      GenericKTProcess( "pptoll", "ɣɣ → l⁺l¯", 4, { { Photon, Photon } }, { Muon, Muon } ),
      y1_( 0. ), y2_( 0. ), pt_diff_( 0. ), phi_pt_diff_( 0. )
    {}

    void
    PPtoLL::preparePhaseSpace()
    {
      jacobian_ = GenericKTProcess::minimalJacobian();

      // Outgoing leptons
      if ( cuts_.cuts.central.count( Cuts::rapidity_single ) == 0
        || !cuts_.cuts.central.at( Cuts::rapidity_single ).valid()) {
        InWarning( "Failed to retrieve a rapidity range for the outgoing leptons from the user configuration!\n\t"
                   "Setting it to the default | y(l) | < 6 value." );
        cuts_.cuts.central[Cuts::rapidity_single] = { -6., 6. };
      }
      rap_limits_ = cuts_.cuts.central.at( Cuts::rapidity_single );
      jacobian_ *= pow( rap_limits_.range(), 2 );

      if ( cuts_.cuts.central.count( Cuts::pt_diff ) == 0
        || !cuts_.cuts.central.at( Cuts::pt_diff ).valid() ) {
        InWarning( "Failed to retrieve a leptons pT difference range from the user configuration!\n\t"
                   "Setting it to the default ΔpT < 50 GeV value." );
        cuts_.cuts.central[Cuts::pt_diff] = { 0., 50. };
      }
      ptdiff_limits_ = cuts_.cuts.central.at( Cuts::pt_diff );
      jacobian_ *= ptdiff_limits_.range();

      if ( cuts_.cuts.central.count( Cuts::phi_pt_diff ) == 0
        || !cuts_.cuts.central.at( Cuts::phi_pt_diff ).valid() ) {
        InWarning( "Failed to retrieve a leptons azimuthal angle difference range from the user configuration!\n\t"
                   "Setting it to the default 0 < Δɸ < 2π value." );
        cuts_.cuts.central[Cuts::phi_pt_diff] = { 0., 2*M_PI };
      }
      phi_pt_diff_limits_ = cuts_.cuts.central.at( Cuts::phi_pt_diff );
      jacobian_ *= phi_pt_diff_limits_.range();
    }

    void
    PPtoLL::prepareKTKinematics()
    {
      y1_ = rap_limits_.x( xkt( 0 ) );
      y2_ = rap_limits_.x( xkt( 1 ) );

      DebuggingInsideLoop( Form( "Leptons rapidities (%.2f < y < %.2f): %f / %f", rap_limits_.min(), rap_limits_.max(), y1_, y2_ ) );

      pt_diff_ = ptdiff_limits_.x( xkt( 2 ) );
      phi_pt_diff_ = phi_pt_diff_limits_.x( xkt( 3 ) );

      DebuggingInsideLoop( Form( "leptons pt difference:\n\t"
                                 "  mag = %f (%.2f < Dpt < %.2f)\n\t"
                                 "  phi = %f",
                                 pt_diff_, ptdiff_limits_.min(), ptdiff_limits_.max(), phi_pt_diff_ ) );
    }

    double
    PPtoLL::computeKTFactorisedMatrixElement()
    {
      const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;
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
      //const double stild = s_/2.*(1+sqrt(1.-(4*pow(mp2, 2))/s_*s_));

      // Inner photons
      const double q1tx = qt1_*cos( phi_qt1_ ), q1ty = qt1_*sin( phi_qt1_ ),
                   q2tx = qt2_*cos( phi_qt2_ ), q2ty = qt2_*sin( phi_qt2_ );
      DebuggingInsideLoop( Form( "q1t(x/y) = %e / %e\n\t"
                                 "q2t(x/y) = %e / %e", q1tx, q1ty, q2tx, q2ty ) );

      // Two-photon system
      const double ptsumx = q1tx+q2tx,
                   ptsumy = q1ty+q2ty,
                   ptsum = sqrt( ptsumx*ptsumx+ptsumy*ptsumy );

      const double ptdiffx = pt_diff_*cos( phi_pt_diff_ ),
                   ptdiffy = pt_diff_*sin( phi_pt_diff_ );

      // Outgoing leptons
      const double pt1x = ( ptsumx+ptdiffx )*0.5, pt1y = ( ptsumy+ptdiffy )*0.5, pt1 = sqrt( pt1x*pt1x+pt1y*pt1y ),
                   pt2x = ( ptsumx-ptdiffx )*0.5, pt2y = ( ptsumy-ptdiffy )*0.5, pt2 = sqrt( pt2x*pt2x+pt2y*pt2y );

      const Kinematics::Limits pt_limits = cuts_.cuts.central[Cuts::pt_single];
      if ( pt_limits.hasMin() && ( pt1 < pt_limits.min() || pt2 < pt_limits.min() ) ) return 0.;
      if ( pt_limits.hasMax() && ( pt1 > pt_limits.max() || pt2 > pt_limits.max() ) ) return 0.;

      // transverse mass for the two leptons
      const double amt1 = sqrt( pt1*pt1+ml2 ),
                   amt2 = sqrt( pt2*pt2+ml2 );

      //=================================================================
      //     a window in transverse momentum difference
      //=================================================================

      const Kinematics::Limits ptdiff_limits = cuts_.cuts.central[Cuts::pt_diff];
      if ( ptdiff_limits.hasMax() && fabs( pt1-pt2 ) > ptdiff_limits.max() ) return 0.;

      //=================================================================
      //     a window in rapidity distance
      //=================================================================

      const double dely = fabs( y1_-y2_ );
      const Kinematics::Limits dely_limits = cuts_.cuts.central[Cuts::rapidity_diff];
      if ( dely_limits.hasMin() && dely < dely_limits.min() ) return 0.;
      if ( dely_limits.hasMax() && dely > dely_limits.max() ) return 0.;

      //=================================================================
      //     auxiliary quantities
      //=================================================================

      const double alpha1 = amt1/sqs_*exp(  y1_ ),
                   alpha2 = amt2/sqs_*exp(  y2_ ),
                   beta1  = amt1/sqs_*exp( -y1_ ),
                   beta2  = amt2/sqs_*exp( -y2_ );
      DebuggingInsideLoop( Form( "Sudakov parameters:\n\t"
                                 "  alpha1/2 = %f / %f\n\t"
                                 "   beta1/2 = %f / %f", alpha1, alpha2, beta1, beta2 ) );

      const double q1t2 = q1tx*q1tx+q1ty*q1ty,
                   q2t2 = q2tx*q2tx+q2ty*q2ty;

      //const double old_x2 = 0.; //FIXME figure out where this comes from
      //const double delta_x1 = (MX_*MX_+q2t2)/((1.-old_x2)*s_);

      //x1 = alpha1+alpha2+delta_x1;
      const double x1 = alpha1+alpha2,
                   x2 = beta1 +beta2;

      /*const double xi_x1 = log10(x1);
      const double xi_x2 = log10(x2);*/

      const double z1p = alpha1/x1, z1m = alpha2/x1,
                   z2p = beta1 /x2, z2m = beta2 /x2;
      DebuggingInsideLoop( Form( "z(1/2)p = %f / %f\n\t"
                                 "z(1/2)m = %f / %f", z1p, z2p, z1m, z2m ) );

      if ( x1 > 1. || x2 > 1. ) return 0.; // sanity check

      // FIXME FIXME FIXME
      const double ak10 = event_->getOneByRole( Particle::IncomingBeam1 ).energy(),
                   ak1z = event_->getOneByRole( Particle::IncomingBeam1 ).momentum().pz(),
                   ak20 = event_->getOneByRole( Particle::IncomingBeam2 ).energy(),
                   ak2z = event_->getOneByRole( Particle::IncomingBeam2 ).momentum().pz();
      DebuggingInsideLoop( Form( "incoming particles: p1: %f / %f\n\t"
                                 "                    p2: %f / %f", ak1z, ak10, ak2z, ak20 ) );

      //=================================================================
      //     additional conditions for energy-momentum conservation
      //=================================================================

      const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
      const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh(y1_-y2_) - ptsum*ptsum );
      DebuggingInsideLoop( Form( "s(1/2)_eff = %f / %f GeV^2\n\t"
                                 "dilepton invariant mass = %f GeV", s1_eff, s2_eff, invm ) );

      switch ( cuts_.mode ) {
        case Kinematics::ElasticInelastic:   if ( sqrt( s1_eff ) <= ( MY_+invm ) ) return 0.;
        case Kinematics::InelasticElastic:   if ( sqrt( s2_eff ) <= ( MX_+invm ) ) return 0.;
        case Kinematics::InelasticInelastic: if ( sqrt( s1_eff ) <= ( MY_+invm ) ) return 0.;
                                             if ( sqrt( s2_eff ) <= ( MX_+invm ) ) return 0.;
        default: break;
      }

      //const double qcaptx = pcaptx, qcapty = pcapty;

      //=================================================================
      //     four-momenta of the outgoing protons (or remnants)
      //=================================================================

      const double px_plus  = ( 1.-x1 )*fabs( ak1z )*sqrt( 2. ),
                   px_minus = ( MX_*MX_ + q1tx*q1tx + q1ty*q1ty )*0.5/px_plus;

      const double py_minus = ( 1.-x2 )*fabs( ak2z )*sqrt( 2. ), // warning! sign of pz??
                   py_plus  = ( MY_*MY_ + q2tx*q2tx + q2ty*q2ty )*0.5/py_minus;

      DebuggingInsideLoop( Form( "px_(+/-) = %f / %f\n\t"
                                 "py_(+/-) = %f / %f", px_plus, px_minus, py_plus, py_minus ) );

      PX_ = Particle::Momentum( -q1tx, -q1ty, 0.5*( px_plus-px_minus )*sqrt( 2. ), 0.5*( px_plus+px_minus )*sqrt( 2. ) );
      PY_ = Particle::Momentum( -q2tx, -q2ty, 0.5*( py_plus-py_minus )*sqrt( 2. ), 0.5*( py_plus+py_minus )*sqrt( 2. ) );

      DebuggingInsideLoop( Form( "First remnant:  (E,p) = (%f, %f, %f, %f), mass = %f\n\t"
                                 "Second remnant: (E,p) = (%f, %f, %f, %f), mass = %f",
                                 PX_.px(), PX_.py(), PX_.pz(), PX_.energy(), PX_.mass(),
                                 PY_.px(), PY_.py(), PY_.pz(), PY_.energy(), PY_.mass() ) );

      assert( fabs( PX_.mass()-MX_ ) < 1.e-6 );
      assert( fabs( PY_.mass()-MY_ ) < 1.e-6 );

      //=================================================================
      //     four-momenta of the outgoing l^+ and l^-
      //=================================================================

      Particle::Momentum p1( pt1x, pt1y, alpha1*ak1z + beta1*ak2z, alpha1*ak10 + beta1*ak20 ),
                         p2( pt2x, pt2y, alpha2*ak1z + beta2*ak2z, alpha2*ak10 + beta2*ak20 );
      DebuggingInsideLoop( Form( "unboosted first lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                                 "          second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                                 p1.px(), p1.py(), p1.pz(), p1.energy(), p1.mass(),
                                 p2.px(), p2.py(), p2.pz(), p2.energy(), p2.mass() ) );

      Pl1_ = Particle::Momentum( pt1x, pt1y, sqrt( pt1*pt1 + ml2 )*sinh( y1_ ), sqrt( pt1*pt1 + ml2 )*cosh( y1_ ) );
      Pl2_ = Particle::Momentum( pt2x, pt2y, sqrt( pt2*pt2 + ml2 )*sinh( y2_ ), sqrt( pt2*pt2 + ml2 )*cosh( y2_ ) );

      DebuggingInsideLoop( Form( "First lepton:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                                 "Second lepton: (E,p), m = (%f, %f, %f, %f), %f",
                                 Pl1_.px(), Pl1_.py(), Pl1_.pz(), Pl1_.energy(), Pl1_.mass(),
                                 Pl2_.px(), Pl2_.py(), Pl2_.pz(), Pl2_.energy(), Pl2_.mass() ) );

      assert( fabs( Pl1_.mass()-event_->getByRole( Particle::CentralSystem )[0].mass() ) < 1.e-6 );
      assert( fabs( Pl2_.mass()-event_->getByRole( Particle::CentralSystem )[1].mass() ) < 1.e-6 );

      //=================================================================
      //     four-momenta squared of the virtual photons
      //=================================================================

      // FIXME FIXME FIXME /////////////////////
      Particle::Momentum q1( q1tx, q1ty, 0., 0. ),
                         q2( q2tx, q2ty, 0., 0. );
      //////////////////////////////////////////

      DebuggingInsideLoop( Form( "First photon*:  (E,p), m2 = (%f, %f, %f, %f), %e\n\t"
                                 "Second photon*: (E,p), m2 = (%f, %f, %f, %f), %e",
                                 q1.px(), q1.py(), q1.pz(), q1.energy(), q1.mass2(),
                                 q2.px(), q2.py(), q2.pz(), q2.energy(), q2.mass2() ) );
      //const double q12 = q1.mass2(), q22 = q2.mass2();

      //=================================================================
      //     Mendelstam variables
      //=================================================================

      //const double shat = s_*x1*x2; // ishat = 1 (approximation)
      //const double shat = ( q1+q2 ).mass2(); // ishat = 2 (exact formula)

      const double that1 = ( q1-p1 ).mass2(), that2 = ( q2-p2 ).mass2(),
                   uhat1 = ( q1-p2 ).mass2(), uhat2 = ( q2-p1 ).mass2();
      DebuggingInsideLoop( Form( "that(1/2) = %f / %f\n\t"
                                 "uhat(1/2) = %f / %f",
                                 that1, that2, uhat1, uhat2 ) );

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

        double aux2_1, aux2_2;

        const double ak1_x = z1m*pt1x - z1p*pt2x, ak1_y = z1m*pt1y - z1p*pt2y,
                     ak2_x = z2m*pt1x - z2p*pt2x, ak2_y = z2m*pt1y - z2p*pt2y;

        const double t1abs = ( q1t2 + x1*( MX_*MX_-mp2 )+x1*x1*mp2 )/( 1.-x1 ),
                     t2abs = ( q2t2 + x2*( MY_*MY_-mp2 )+x2*x2*mp2 )/( 1.-x2 );

        const double eps12 = ml2 + z1p*z1m*t1abs,
                     eps22 = ml2 + z2p*z2m*t2abs;

        const double Phi10 = 1./( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                            -1./( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi11_x = ( ak1_x + z1p*q2tx )/( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                              -( ak1_x - z1m*q2tx )/( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi11_y = ( ak1_y + z1p*q2ty )/( pow( ak1_x + z1p*q2tx, 2 ) + pow( ak1_y + z1p*q2ty, 2 ) + eps12 )
                              -( ak1_y - z1m*q2ty )/( pow( ak1_x - z1m*q2tx, 2 ) + pow( ak1_y - z1m*q2ty, 2 ) + eps12 ),
                     Phi102 = Phi10*Phi10;
        //const double Phi112 = Phi11_x*Phi11_x + Phi11_y*Phu11_y;

        const double Phi20 = 1./( pow( ak2_x + z2p*q1tx, 2 )+pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                            -1./( pow( ak2_x - z2m*q1tx, 2 )+pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi21_x = ( ak2_x + z2p*q1tx )/( pow( ak2_x + z2p*q1tx, 2 ) + pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                              -( ak2_x - z2m*q1tx )/( pow( ak2_x - z2m*q1tx ,2 ) + pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi21_y = ( ak2_y + z2p*q1ty )/( pow( ak2_x + z2p*q1tx, 2 ) + pow( ak2_y + z2p*q1ty, 2 ) + eps22 )
                              -( ak2_y - z2m*q1ty )/( pow( ak2_x - z2m*q1tx, 2 ) + pow( ak2_y - z2m*q1ty, 2 ) + eps22 ),
                     Phi202 = Phi20*Phi20;
        //const double Phi212 = Phi21_x*Phi21_x + Phi21_y*Phi21_y;

        /*aux2_1 = iterm11*( ml2 + 4.*z1p*z1m*t1abs )*Phi102
                +iterm22*( z1p*z1p + z1m*z1m )*Phi112
                -iterm12*4.*z1p*z1m*( z1p-z1m )*Phi10*( q1tx*Phi11_x + q1ty*Phi11_y );
        aux2_2 = iterm11*( ml2+4.*z2p*z2m*t2abs)*Phi202
                +iterm22*( z2p*z2p + z2m*z2m )*Phi212
                -iterm12*4.*z2p*z2m*( z2p-z2m )*Phi20*( q2tx*Phi21_x + q2ty*Phi21_y );*/

        const double Phi11_dot_e = ( Phi11_x*q1tx + Phi11_y*q1ty )/qt1_, Phi11_cross_e = ( Phi11_x*q1ty-Phi11_y*q1tx )/qt1_;
        const double Phi21_dot_e = ( Phi21_x*q2tx + Phi21_y*q2ty )/qt2_, Phi21_cross_e = ( Phi21_x*q2ty-Phi21_y*q2tx )/qt2_;
        DebuggingInsideLoop( Form( "Phi1: E, px, py = %e, %e, %e\n\t"
                                   "Phi2: E, px, py = %e, %e, %e\n\t"
                                   "(dot):   %e / %e\n\t"
                                   "(cross): %e / %e",
                                   Phi10, Phi11_x, Phi11_y, Phi20, Phi21_x, Phi21_y,
                                   Phi11_dot_e, Phi21_dot_e, Phi11_cross_e, Phi21_cross_e ) );

        aux2_1 = iterm11 * ( ml2 + 4.*z1p*z1p*z1m*z1m*t1abs ) * Phi102
                +iterm22 * ( ( z1p*z1p + z1m*z1m )*( Phi11_dot_e*Phi11_dot_e + Phi11_cross_e*Phi11_cross_e ) )
                +itermtt * ( Phi11_cross_e*Phi11_cross_e - Phi11_dot_e*Phi11_dot_e )
                -iterm12 * 4.*z1p*z1m*( z1p-z1m ) * Phi10 * ( q1tx*Phi11_x + q1ty*Phi11_y );

        aux2_2 = iterm11 * ( ml2 + 4.*z2p*z2p*z2p*z2m*t2abs ) * Phi202
                +iterm22 * ( ( z2p*z2p + z2m*z2m )*( Phi21_dot_e*Phi21_dot_e + Phi21_cross_e*Phi21_cross_e ) )
                +itermtt * ( Phi21_cross_e*Phi21_cross_e - Phi21_dot_e*Phi21_dot_e )
                -iterm12 * 4.*z2p*z2m*( z2p-z2m ) * Phi20 * ( q2tx*Phi21_x + q2ty*Phi21_y );

        //=================================================================
        //     convention of matrix element as in our kt-factorization
        //     for heavy flavours
        //=================================================================

        const double norm = 16.*M_PI*M_PI*Constants::alphaEM*Constants::alphaEM;

        double amat2_1 = norm * pow( x1*x2*s_, 2 ) * aux2_1*2.*z1p*z1m*t1abs/( q1t2*q2t2 ) * t2abs/q2t2;
        double amat2_2 = norm * pow( x1*x2*s_, 2 ) * aux2_2*2.*z2p*z2m*t2abs/( q1t2*q2t2 );

        //=================================================================
        //     symmetrization
        //=================================================================

        amat2 = 0.5*( imat1*amat2_1 + imat2*amat2_2 );

        DebuggingInsideLoop( Form( "aux2(1/2) = %e / %e\n\t"
                                   "amat2(1/2), amat2 = %e / %e / %e", aux2_1, aux2_2, amat2_1, amat2_2, amat2 ) );
        /*const double xx1 = alpha1+alpha2, xx2 = beta1+beta2;

        const double sudakov_2 = ( MX_*MX_ - mp2+q2t2+xx2*mp2 )/( ( 1.-xx2 )*s_ );
        const double sudakov_1 = ( q1t2 + xx1*mp2 )/( ( 1.-xx1 )*s_ );
        const double ratio1 = sudakov_1 / xx1,
                     ratio2 = sudakov_2 / xx2;*/

        //if ( ratio1 > 0.01 ) return 0.;
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

      const double aintegral = amat2 / ( 16.*M_PI*M_PI*x1*x1*x2*x2*s_*s_ )
                             * flux1_/M_PI * flux2_/M_PI
                             * Constants::GeV2toBarn * 0.25;

      //=================================================================
      return aintegral*qt1_*qt2_*pt_diff_;
      //=================================================================
    }

    void
    PPtoLL::fillCentralParticlesKinematics()
    {
      // randomise the charge of the outgoing leptons
      int sign = ( drand()>.5 ) ? +1 : -1;

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
