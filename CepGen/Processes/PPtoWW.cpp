#include "PPtoWW.h"
#include "CepGen/Core/Exception.h"
#include <assert.h>

using namespace CepGen::Process;

PPtoWW::PPtoWW() :
  GenericKTProcess( "pptoww", "ɣɣ → W⁺W¯", 4, { { Photon, Photon } }, { W, W } )
{}

void
PPtoWW::prepareKTKinematics()
{
  const Kinematics::Limits rap_limits = cuts_.cuts.central[Cuts::rapidity_single];

  // outgoing Ws
  y1_ = rap_limits.x( xkt( 0 ) );
  y2_ = rap_limits.x( xkt( 1 ) );
  DebuggingInsideLoop( Form( "W bosons rapidities (%.2f < y < %.2f): %f / %f", rap_limits.min(), rap_limits.max(), y1_, y2_ ) );

  Kinematics::Limits ptdiff_limits = cuts_.cuts.central[Cuts::pt_diff];
  if ( !ptdiff_limits.hasMax() ) ptdiff_limits.max() = 500.; //FIXME
  pt_diff_ = ptdiff_limits.x( xkt( 2 ) );

  phi_pt_diff_ = 2.*M_PI*xkt( 3 );

  DebuggingInsideLoop( Form( "W bosons pt difference:\n\t"
                             "  mag = %f (%.2f < Dpt < %.2f)\n\t"
                             "  phi = %f",
                             pt_diff_, ptdiff_limits.min(), ptdiff_limits.max(), phi_pt_diff_ ) );
}

double
PPtoWW::computeJacobian()
{
  double jac = GenericKTProcess::minimalJacobian();
  jac *= cuts_.cuts.central[Cuts::rapidity_single].range(); // d(y1)
  jac *= cuts_.cuts.central[Cuts::rapidity_single].range(); // d(y2)
  jac *= cuts_.cuts.central[Cuts::pt_diff].range(); // d(Dpt)
  jac *= 2.*M_PI; // d(phiDpt)

  return jac;
}

double
PPtoWW::computeKTFactorisedMatrixElement()
{
  const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;
  const double mw = ParticleProperties::mass( W ), mw2 = mw*mw;

  //=================================================================
  //     How matrix element is calculated
  //=================================================================

  const unsigned short method = 1;

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
  const double pt1x = 0.5 * ( ptsumx+ptdiffx ), pt1y = 0.5 * ( ptsumy+ptdiffy ), pt1 = sqrt( pt1x*pt1x+pt1y*pt1y ),
               pt2x = 0.5 * ( ptsumx-ptdiffx ), pt2y = 0.5 * ( ptsumy-ptdiffy ), pt2 = sqrt( pt2x*pt2x+pt2y*pt2y );

  const Kinematics::Limits pt_limits = cuts_.cuts.central[Cuts::pt_single];
  if ( pt_limits.hasMin() && ( pt1 < pt_limits.min() || pt2 < pt_limits.min() ) ) return 0.;
  if ( pt_limits.hasMax() && ( pt1 > pt_limits.max() || pt2 > pt_limits.max() ) ) return 0.;

  // transverse mass for the two leptons
  const double amt1 = sqrt( pt1*pt1+mw2 ),
               amt2 = sqrt( pt2*pt2+mw2 );

  //=================================================================
  //     a window in two-boson invariant mass
  //=================================================================

  const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh( y1_-y2_ ) - ptsum*ptsum );
  const Kinematics::Limits invm_limits = cuts_.cuts.central[Cuts::mass_sum];
  if ( invm_limits.hasMin() && invm < invm_limits.min() ) return 0.;
  if ( invm_limits.hasMax() && invm > invm_limits.max() ) return 0.;

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

  const double alpha1 = amt1/sqs_*exp( +y1_ ),
               alpha2 = amt2/sqs_*exp( +y2_ ),
               beta1  = amt1/sqs_*exp( -y1_ ),
               beta2  = amt2/sqs_*exp( -y2_ );
  DebuggingInsideLoop( Form( "Sudakov parameters:\n\t"
                             "  alpha1/2 = %f / %f\n\t"
                             "   beta1/2 = %f / %f", alpha1, alpha2, beta1, beta2 ) );

  const double q1t2 = q1tx*q1tx + q1ty*q1ty,
               q2t2 = q2tx*q2tx + q2ty*q2ty;

  //const double old_x2 = 0.; //FIXME figure out where this comes from
  //const double delta_x1 = (MX_*MX_+q2t2)/((1.-old_x2)*s_);

  //x1 = alpha1+alpha2+delta_x1;
  const double x1 = alpha1 + alpha2,
               x2 = beta1  + beta2;

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
               px_minus = ( MX_*MX_ + q1t2 )*0.5/px_plus;

  const double py_minus = ( 1.-x2 )*fabs( ak2z )*sqrt( 2. ), // warning! sign of pz??
               py_plus  = ( MY_*MY_ + q2t2 )*0.5/py_minus;

  DebuggingInsideLoop( Form( "px_(+/-) = %f / %f\n\t"
                             "py_(+/-) = %f / %f", px_plus, px_minus, py_plus, py_minus ) );

  PX_ = Particle::Momentum( -q1tx, -q1ty, ( px_plus-px_minus )/sqrt( 2. ), ( px_plus+px_minus )/sqrt( 2. ) );
  PY_ = Particle::Momentum( -q2tx, -q2ty, ( py_plus-py_minus )/sqrt( 2. ), ( py_plus+py_minus )/sqrt( 2. ) );

  DebuggingInsideLoop( Form( "First remnant:  (E,p) = (%f, %f, %f, %f)\n\t"
                             "Second remnant: (E,p) = (%f, %f, %f, %f)",
                             PX_.px(), PX_.py(), PX_.pz(), PX_.energy(),
                             PY_.px(), PY_.py(), PY_.pz(), PY_.energy() ) );

  /*assert( fabs( PX_.mass()-MX_ ) < 1.e-6 );
  assert( fabs( PY_.mass()-MY_ ) < 1.e-6 );*/

  //=================================================================
  //     four-momenta squared of the virtual photons
  //=================================================================

  const double ww = 0.5 * ( 1.+sqrt( 1.-4.*mp2/s_ ) );

  // FIXME FIXME FIXME /////////////////////
  Particle::Momentum q1( q1tx, q1ty, +0.5 * x1*ww*sqs_*( 1.-q1t2/x1/x1/ww/ww/s_ ), 0.5 * x1*ww*sqs_*( 1.+q1t2/x1/x1/ww/ww/s_ ) ),
                     q2( q2tx, q2ty, -0.5 * x2*ww*sqs_*( 1.-q2t2/x2/x2/ww/ww/s_ ), 0.5 * x2*ww*sqs_*( 1.+q2t2/x2/x2/ww/ww/s_ ) );
  //////////////////////////////////////////

  DebuggingInsideLoop( Form( "First photon*:  (E,p), m2 = (%f, %f, %f, %f), %e\n\t"
                             "Second photon*: (E,p), m2 = (%f, %f, %f, %f), %e",
                             q1.px(), q1.py(), q1.pz(), q1.energy(), q1.mass2(),
                             q2.px(), q2.py(), q2.pz(), q2.energy(), q2.mass2() ) );
  //const double q12 = q1.mass2(), q22 = q2.mass2();

  //=================================================================
  //     four-momenta of the outgoing W^+ and W^-
  //=================================================================

  p_w1_ = Particle::Momentum( pt1x, pt1y, alpha1*ak1z + beta1*ak2z, alpha1*ak10 + beta1*ak20 );
  p_w2_ = Particle::Momentum( pt2x, pt2y, alpha2*ak1z + beta2*ak2z, alpha2*ak10 + beta2*ak20 );

  DebuggingInsideLoop( Form( "First W:  (E,p), m = (%f, %f, %f, %f), %f\n\t"
                             "Second W: (E,p), m = (%f, %f, %f, %f), %f",
                             p_w1_.px(), p_w2_.py(), p_w1_.pz(), p_w1_.energy(), p_w1_.mass(),
                             p_w2_.px(), p_w2_.py(), p_w2_.pz(), p_w2_.energy(), p_w2_.mass() ) );

  //assert( fabs( p_w1_.mass()-event_->getByRole( Particle::CentralSystem )[0].mass() ) < 1.e-6 );
  //assert( fabs( p_w2_.mass()-event_->getByRole( Particle::CentralSystem )[1].mass() ) < 1.e-6 );

  //=================================================================
  //     Mendelstam variables
  //=================================================================

  //const double shat = s_*x1*x2; // ishat = 1 (approximation)
  const double shat = ( q1+q2 ).mass2(); // ishat = 2 (exact formula)

  const double that1 = ( q1-p_w1_ ).mass2(), that2 = ( q2-p_w2_ ).mass2(),
               uhat1 = ( q1-p_w2_ ).mass2(), uhat2 = ( q2-p_w1_ ).mass2();
  DebuggingInsideLoop( Form( "that(1/2) = %f / %f\n\t"
                             "uhat(1/2) = %f / %f",
                             that1, that2, uhat1, uhat2 ) );

  //const double mll = sqrt( shat );

  const double that = 0.5*( that1+that2 ), uhat = 0.5*( uhat1+uhat2 );

  //=================================================================
  //     matrix elements
  //=================================================================
  double amat2 = 0.;
  if ( method == 0 ) {

    //=================================================================
    //     matrix element for gamma gamma --> W^+ W^-
    //     (Denner+Dittmaier+Schuster)
    //     (work in collaboration with C. Royon)
    //=================================================================

    const double mw4 = mw2*mw2;

    const double term1  = 2.*shat * ( 2.*shat+3.*mw2 ) / ( 3.*( mw2-that )*( mw2-uhat ) );
    const double term2  = 2.*shat*shat * ( shat*shat + 3.*mw4 ) / ( 3.*pow( mw2-that, 2 )*pow( mw2-uhat, 2 ) );

    const double auxil_gamgam = 1. - term1 + term2;
    const double beta = sqrt( 1.-4.*mw2/shat );

    amat2 = 3.*Constants::alphaEM*Constants::alphaEM*beta / ( 2.*shat ) * auxil_gamgam / ( beta/( 64.*M_PI*M_PI*shat ) );
  }
  else if ( method == 1 ) {

    //=================================================================
    //     off-shell Nachtmann formulae
    //=================================================================

    const double e2 = 4.*M_PI*Constants::alphaEM;

    const double phi_diff = phi_qt1_-phi_qt2_, phi_sum = phi_qt1_+phi_qt2_;
    double amat2_0 = 0., amat2_1 = 0., amat2_interf = 0.;
    for ( const auto lam3 : { -1, 0, 1 } ) {
      for ( const auto lam4 : { -1, 0, 1 } ) {
        double ampli_pp = WWamplitude( shat, that, uhat, +1, +1, lam3, lam4 );
        double ampli_mm = WWamplitude( shat, that, uhat, -1, -1, lam3, lam4 );
        double ampli_pm = WWamplitude( shat, that, uhat, +1, -1, lam3, lam4 );
        double ampli_mp = WWamplitude( shat, that, uhat, -1, +1, lam3, lam4 );
        amat2_0 += ampli_pp*ampli_pp + ampli_mm*ampli_mm + 2.*cos( 2.*phi_diff )*ampli_pp*ampli_mm;
        amat2_1 += ampli_pm*ampli_pm + ampli_mp*ampli_mp + 2.*cos( 2.*phi_sum  )*ampli_pm*ampli_mp;
        amat2_interf -= 2.*( cos( phi_sum+phi_diff )*( ampli_pp*ampli_pm+ampli_mm*ampli_mp ) + cos( phi_sum-phi_diff )*( ampli_pp*ampli_mp+ampli_mm*ampli_pm ) );
      }
    }
    amat2 = e2*e2 * ( amat2_0 + amat2_1 + amat2_interf );
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

  const double aintegral = amat2 * ( 2.*M_PI )/ ( 16.*M_PI*M_PI*( x1*x2*s_ )*( x1*x2*s_ ) )
                         * flux1_/M_PI * flux2_/M_PI * 0.25
                         * Constants::GeV2toBarn * 0.5 / M_PI;
  /*const double aintegral = amat2 / ( 16.*M_PI*M_PI*x1*x1*x2*x2*s_*s_ )
                         * flux1_/M_PI * flux2_/M_PI
                         * Constants::GeV2toBarn * 0.25;*/

  //=================================================================
  return aintegral*qt1_*qt2_*pt_diff_;
  //=================================================================
}

void
PPtoWW::fillCentralParticlesKinematics()
{
  // randomise the charge of the outgoing leptons
  short sign = ( drand()>.5 ) ? +1 : -1;

  //=================================================================
  //     first outgoing lepton
  //=================================================================
  Particle& ow1 = event_->getByRole( Particle::CentralSystem )[0];
  ow1.setPdgId( ow1.pdgId(), sign );
  ow1.setStatus( Particle::Undecayed );
  ow1.setMomentum( p_w1_ );

  //=================================================================
  //     second outgoing lepton
  //=================================================================
  Particle& ow2 = event_->getByRole( Particle::CentralSystem )[1];
  ow2.setPdgId( ow2.pdgId(), -sign );
  ow2.setStatus( Particle::Undecayed );
  ow2.setMomentum( p_w2_ );
}

double
PPtoWW::WWamplitude( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 ) const
{
  const double mw = ParticleProperties::mass( W ), mw2 = mw*mw;
  const double sqrt2 = sqrt( 2. );

  // then compute some kinematic variables
  const double cos_theta = ( that-uhat ) / shat / sqrt( 1.+1.e-10-4.*mw2/shat ), cos_theta2 = cos_theta*cos_theta;
  const double sin_theta2 = 1.-cos_theta2, sin_theta = sqrt( sin_theta2 );
  const double beta = sqrt( 1.-4.*mw2/shat ), beta2 = beta*beta;
  const double gamma = 1./sqrt( 1.-beta2 ), gamma2 = gamma*gamma;
  const double invA = 1./( 1.-beta2*cos_theta2 );

  const double term1 = 1./gamma2*( ( gamma2+1. )*( 1.-lam1*lam2 )*sin_theta2 - ( 1.+lam1*lam2 ) );
  const double term2 = -sqrt2/gamma*( lam1-lam2 ) * ( 1.+lam1*lam3*cos_theta )*sin_theta;
  const double term3 = -0.5*( 2.*beta*( lam1+lam2 )*( lam3+lam4 ) - ( 1./gamma2 )*( 1.+lam3*lam4 )*( 2.*lam1*lam2+( 1.-lam1*lam2 ) * cos_theta2 )+( 1.+lam1*lam2*lam3*lam4 )*( 3.+lam1*lam2 ) + 2.*( lam1-lam2 )*( lam3-lam4 )*cos_theta + ( 1.-lam1*lam2 )*( 1.-lam3*lam4 )*cos_theta2 );
  const double term4 = -sqrt2/gamma*( lam2-lam1 )*( 1.+lam2*lam4*cos_theta )*sin_theta;

  if ( lam3 == 0 && lam4 == 0 ) return invA*term1;
  if ( lam4 == 0 )              return invA*term2;
  if ( lam3 == 0 )              return invA*term4;
  if ( lam3 != 0 && lam4 != 0 ) return invA*term3;
  return 0.;
}
