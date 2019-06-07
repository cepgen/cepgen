#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Processes/ProcessesHandler.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include <iomanip>

namespace cepgen
{
  namespace proc
  {
    PPtoFF::PPtoFF( const ParametersList& params ) :
      GenericKTProcess( params, "pptoff", "ɣɣ → f⁺f¯", { { PDG::photon, PDG::photon } }, { PDG::muon, PDG::muon } ),
      pair_  ( params.get<ParametersList>( "pair" ).get<int>( "pdgid" ) ),
      method_( (ME)params.get<int>( "method", (int)ME::offShell ) ),
      p_mat1_( 0 ), p_mat2_( 0 ), p_term_ll_( 0 ), p_term_lt_( 0 ), p_term_tt1_( 0 ), p_term_tt2_( 0 ),
      y1_( 0. ), y2_( 0. ), pt_diff_( 0. ), phi_pt_diff_( 0. )
    {
      if ( method_ == ME::offShell ) { // off-shell matrix element
        const auto& ofp = params.get<ParametersList>( "offShellParameters" );
        p_mat1_ = ofp.get<int>( "mat1", 1 );
        p_mat2_ = ofp.get<int>( "mat2", 1 );
        p_term_ll_ = ofp.get<int>( "termLL", 1 );
        p_term_lt_ = ofp.get<int>( "termLT", 1 );
        p_term_tt1_ = ofp.get<int>( "termTT", 1 );
        p_term_tt2_ = ofp.get<int>( "termtt", 1 );
      }
    }

    void
    PPtoFF::preparePhaseSpace()
    {
      registerVariable( y1_, Mapping::linear, kin_.cuts.central.rapidity_single, { -6., 6. }, "First outgoing fermion rapidity" );
      registerVariable( y2_, Mapping::linear, kin_.cuts.central.rapidity_single, { -6., 6. }, "Second outgoing fermion rapidity" );
      registerVariable( pt_diff_, Mapping::linear, kin_.cuts.central.pt_diff, { 0., 50. }, "Fermions transverse momentum difference" );
      registerVariable( phi_pt_diff_, Mapping::linear, kin_.cuts.central.phi_pt_diff, { 0., 2.*M_PI }, "Fermions azimuthal angle difference" );

      const auto& pair_info = PDG::get()( pair_ ); // all properties on the fermion pair
      if ( !pair_info.fermion || pair_info.charge == 0. )
        throw CG_FATAL( "PPtoFF:prepare" )
          << "Invalid fermion pair selected: " << PDG::get()( pair_ ).description
          << " (" << (int)pair_ << ")!";

      mf_ = pair_info.mass; mf2_ = mf_*mf_;
      qf_ = pair_info.charge;
      colf_ = pair_info.colours;
      CG_DEBUG( "PPtoFF:prepare" )
        << "Produced particles: " << pair_ << " ("
        << "mass = " << mf_ << " GeV, "
        << "charge = " << std::setprecision( 2 ) << qf_ << " e)\n"
        << "matrix element computation method: " << (int)method_ << ".";
    }

    double
    PPtoFF::computeKTFactorisedMatrixElement()
    {
      //=================================================================
      //     matrix element computation
      //=================================================================

      //--- incoming photons (in two-photon frame, hence fully transverse)
      const auto q1t = Particle::Momentum::fromPThetaPhi( qt1_, M_PI_2, phi_qt1_ );
      const auto q2t = Particle::Momentum::fromPThetaPhi( qt2_, M_PI_2, phi_qt2_ );

      CG_DEBUG_LOOP( "PPtoFF" ) << "q(1/2)t = " << q1t << ", " << q2t << ".";

      //--- two-photon system
      const Particle::Momentum ptsum = q1t+q2t;
      const auto ptdiff = Particle::Momentum::fromPThetaPhi( pt_diff_, M_PI_2, phi_pt_diff_ );

      //--- outgoing fermions
      const Particle::Momentum p1_cm = 0.5*( ptsum+ptdiff ), p2_cm = 0.5*( ptsum-ptdiff );

      //=================================================================
      //     a window in single particle transverse momentum
      //=================================================================

      const auto& pt_limits = kin_.cuts.central.pt_single;
      if ( !pt_limits.passes( p1_cm.pt() ) || !pt_limits.passes( p2_cm.pt() ) )
        return 0.;

      //=================================================================
      //     a window in transverse momentum difference
      //=================================================================

      if ( !kin_.cuts.central.pt_diff.passes( fabs( p1_cm.pt()-p2_cm.pt() ) ) )
        return 0.;

      //=================================================================
      //     a window in rapidity distance
      //=================================================================

      if ( !kin_.cuts.central.rapidity_diff.passes( fabs( y1_-y2_ ) ) )
        return 0.;

      //=================================================================
      //     auxiliary quantities
      //=================================================================

      // transverse mass for the two fermions
      const double amt1 = std::hypot( p1_cm.pt(), mf_ ), amt2 = std::hypot( p2_cm.pt(), mf_ );
      const double alpha1 = amt1/sqs_*exp( y1_ ), beta1  = amt1/sqs_*exp( -y1_ ),
                   alpha2 = amt2/sqs_*exp( y2_ ), beta2  = amt2/sqs_*exp( -y2_ );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "Sudakov parameters:\n\t"
        << "  alpha(1/2) = " << alpha1 << ", " << alpha2 << "\n\t"
        << "   beta(1/2) = " << beta1 << ", " << beta2 << ".";

      const double x1 = alpha1+alpha2, x2 = beta1+beta2;
      { // sanity check for x_i values
        const Limits x_limits{ 0., 1. };
        if ( !x_limits.passes( x1 ) || !x_limits.passes( x2 ) )
          return 0.;
      }

      //=================================================================
      //     additional conditions for energy-momentum conservation
      //=================================================================

      const double s1_eff = x1*s_-qt1_*qt1_, s2_eff = x2*s_-qt2_*qt2_;
      const double invm = sqrt( amt1*amt1 + amt2*amt2 + 2.*amt1*amt2*cosh(y1_-y2_) - ptsum.pt2() );
      CG_DEBUG_LOOP( "PPtoFF" )
        << "s(1/2)eff = " << s1_eff << ", " << s2_eff << " GeV²\n\t"
        << "central system's invariant mass = " << invm << " GeV.";

      if ( ( kin_.mode == KinematicsMode::ElasticInelastic
          || kin_.mode == KinematicsMode::InelasticInelastic )
        && ( sqrt( s1_eff ) <= ( MY_+invm ) ) )
        return 0.;
      if ( ( kin_.mode == KinematicsMode::InelasticElastic
          || kin_.mode == KinematicsMode::InelasticInelastic )
        && ( sqrt( s2_eff ) <= ( MX_+invm ) ) )
        return 0.;

      //=================================================================
      //     four-momenta of the outgoing protons (or remnants)
      //=================================================================

      const Particle::Momentum& ak1 = event_->getOneByRole( Particle::IncomingBeam1 ).momentum(),
                               &ak2 = event_->getOneByRole( Particle::IncomingBeam2 ).momentum();
      CG_DEBUG_LOOP( "PPtoFF" )
        << "incoming particles: p(1/2) = " << ak1 << ", " << ak2 << ".";

      const double px_plus  = ( 1.-x1 )*M_SQRT2*ak1.p(),
                   px_minus = ( MX_*MX_ + q1t.pt2() )*0.5/px_plus;

      const double py_minus = ( 1.-x2 )*M_SQRT2*ak2.p(),
                   py_plus  = ( MY_*MY_ + q2t.pt2() )*0.5/py_minus;

      CG_DEBUG_LOOP( "PPtoFF" )
        << "px± = " << px_plus << ", " << px_minus << "\n\t"
        << "py± = " << py_plus << ", " << py_minus << ".";

      PX_ = Particle::Momentum( -q1t.px(), -q1t.py(), ( px_plus-px_minus )*M_SQRT1_2, ( px_plus+px_minus )*M_SQRT1_2 );
      PY_ = Particle::Momentum( -q2t.px(), -q2t.py(), ( py_plus-py_minus )*M_SQRT1_2, ( py_plus+py_minus )*M_SQRT1_2 );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "First remnant:  " << PX_ << ", mass = " << PX_.mass() << "\n\t"
        << "Second remnant: " << PY_ << ", mass = " << PY_.mass() << ".";

      if ( fabs( PX_.mass()-MX_ ) > 1.e-4 )
        throw CG_FATAL( "PPtoFF" ) << "Invalid X system mass: " << PX_.mass() << "/" << MX_ << ".";
      if ( fabs( PY_.mass()-MY_ ) > 1.e-4 )
        throw CG_FATAL( "PPtoFF" ) << "Invalid Y system mass: " << PY_.mass() << "/" << MY_ << ".";

      //=================================================================
      //     four-momenta of the outgoing l^+ and l^-
      //=================================================================

      const Particle::Momentum p1 = p1_cm + alpha1*ak1 + beta1*ak2;
      const Particle::Momentum p2 = p2_cm + alpha2*ak1 + beta2*ak2;
      CG_DEBUG_LOOP( "PPtoFF" )
        << "unboosted first fermion:  " << p1 << ", mass = " << p1.mass() << "\n\t"
        << "          second fermion: " << p2 << ", mass = " << p2.mass() << ".";

      p_f1_ = Particle::Momentum::fromPxPyYM( p1_cm.px(), p1_cm.py(), y2_, mf_ );
      p_f2_ = Particle::Momentum::fromPxPyYM( p2_cm.px(), p2_cm.py(), y1_, mf_ );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "First fermion:  " << p_f1_ << ", mass = " << p_f1_.mass() << "\n\t"
        << "Second fermion: " << p_f2_ << ", mass = " << p_f2_.mass() << ".";

      if ( fabs( p_f1_.mass()-mf_ ) > 1.e-4 )
        throw CG_FATAL( "PPtoFF" ) << "Invalid fermion 1 mass: "
          << p_f1_.mass() << "/" << mf_ << ".";
      if ( fabs( p_f2_.mass()-mf_ ) > 1.e-4 )
        throw CG_FATAL( "PPtoFF" ) << "Invalid fermion 2 mass: "
          << p_f2_.mass() << "/" << mf_ << ".";

      //=================================================================
      //     matrix elements
      //=================================================================
      double amat2 = 0.;

      //=================================================================
      //     How matrix element is calculated
      //=================================================================

      switch ( method_ ) {
        case ME::onShell: {
          //--- first compute Mendelstam variables
          //const double shat = s_*x1*x2; // approximation
          const double shat = ( q1t+q2t ).mass2(); // exact formula
          const double that1 = ( q1t-p1 ).mass2(), that2 = ( q2t-p2 ).mass2();
          const double uhat1 = ( q1t-p2 ).mass2(), uhat2 = ( q2t-p1 ).mass2();
          const double that = 0.5*( that1+that2 ), uhat = 0.5*( uhat1+uhat2 );

          amat2 = onShellME( shat, that, uhat );

          CG_DEBUG_LOOP( "PPtoFF:onShell" )
            << "that(1/2) = " << that1 << " / " << that2 << "\n\t"
            << "uhat(1/2) = " << uhat1 << " / " << uhat2 << "\n\t"
            << "squared matrix element: " << amat2 << ".";
        } break;
        case ME::offShell: {
          const double t1abs = ( q1t.pt2() + x1*( MX_*MX_-mp2_ )+x1*x1*mp2_ )/( 1.-x1 ),
                       t2abs = ( q2t.pt2() + x2*( MY_*MY_-mp2_ )+x2*x2*mp2_ )/( 1.-x2 );
          const double z1p = alpha1/x1, z1m = alpha2/x1,
                       z2p =  beta1/x2, z2m =  beta2/x2;
          CG_DEBUG_LOOP( "PPtoFF:offShell" )
            << "z(1/2)p = " << z1p << ", " << z2p << "\n\t"
            << "z(1/2)m = " << z1m << ", " << z2m << ".";
          amat2 = offShellME( t1abs, t2abs, z1m, z1p, z2m, z2p, q1t, q2t ) * pow( x1*x2*s_, 2 );
        } break;
      }

      //============================================
      //     unintegrated photon distributions
      //============================================

      const std::pair<double,double> fluxes
        = GenericKTProcess::incomingFluxes( x1, q1t.pt2(), x2, q2t.pt2() );

      CG_DEBUG_LOOP( "PPtoFF" )
        << "Incoming photon fluxes for (x/kt2) = "
        << "(" << x1 << "/" << q1t.pt2() << "), "
        << "(" << x2 << "/" << q2t.pt2() << "):\n\t"
        << fluxes.first << ", " << fluxes.second << ".";

      //=================================================================
      //     factor 2.*pi from integration over phi_sum
      //     factor 1/4 from jacobian of transformations
      //     factors 1/pi and 1/pi due to integration over
      //       d^2 kappa_1 d^2 kappa_2 instead d kappa_1^2 d kappa_2^2
      //=================================================================

      const double g_em = 4.*M_PI*constants::ALPHA_EM*qf_*qf_;
      const double aintegral = amat2 * colf_ * ( g_em*g_em )
                             * 1. / pow( 4.*M_PI*( x1*x2*s_ ), 2 )
                             * fluxes.first*M_1_PI * fluxes.second*M_1_PI * 0.25
                             * constants::GEVM2_TO_PB;

      //=================================================================
      return aintegral*qt1_*qt2_*pt_diff_;
      //=================================================================
    }

    void
    PPtoFF::fillCentralParticlesKinematics()
    {
      // randomise the charge of the outgoing fermions
      const short sign = ( drand() > 0.5 ) ? +1 : -1;

      //=================================================================
      //     first outgoing fermion
      //=================================================================
      Particle& of1 = (*event_)[Particle::CentralSystem][0];
      of1.setPdgId( pair_, sign );
      of1.setStatus( Particle::Status::FinalState );
      of1.setMomentum( p_f1_ );

      //=================================================================
      //     second outgoing fermion
      //=================================================================
      Particle& of2 = (*event_)[Particle::CentralSystem][1];
      of2.setPdgId( pair_, -sign );
      of2.setStatus( Particle::Status::FinalState );
      of2.setMomentum( p_f2_ );
    }

    double
    PPtoFF::onShellME( double shat, double that, double uhat ) const
    {
      CG_DEBUG_LOOP( "PPtoFF:onShell" )
        << "shat: " << shat << ", that: " << that << ", uhat: " << uhat << ".";

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

      return -2.*( term1+term2+term3+term4+term5
                  +term6+term7+term8+term9+term10 )/( pow( ( mf2_-that )*( mf2_-uhat ), 2) );
    }

    double
    PPtoFF::offShellME( double t1abs, double t2abs,
                        double z1m, double z1p, double z2m, double z2p,
                        const Particle::Momentum& q1, const Particle::Momentum& q2 ) const
    {
      const double z1 = z1p*z1m, z2 = z2p*z2m;
      const double eps12 = mf2_+z1*t1abs, eps22 = mf2_+z2*t2abs;

      const Particle::Momentum ak1 = ( z1m*p_f1_-z1p*p_f2_ ), ak2 = ( z2m*p_f1_-z2p*p_f2_ );
      const Particle::Momentum ph_p1 = ak1+z1p*q2, ph_m1 = ak1-z1m*q2;
      const Particle::Momentum ph_p2 = ak2+z2p*q1, ph_m2 = ak2-z2m*q1;

      const Particle::Momentum phi1(
        ph_p1.px()/( ph_p1.pt2()+eps12 )-ph_m1.px()/( ph_m1.pt2()+eps12 ),
        ph_p1.py()/( ph_p1.pt2()+eps12 )-ph_m1.py()/( ph_m1.pt2()+eps12 ),
        0.,
        1./( ph_p1.pt2()+eps12 )-1./( ph_m1.pt2()+eps12 )
      );

      const Particle::Momentum phi2(
        ph_p2.px()/( ph_p2.pt2()+eps22 )-ph_m2.px()/( ph_m2.pt2()+eps22 ),
        ph_p2.py()/( ph_p2.pt2()+eps22 )-ph_m2.py()/( ph_m2.pt2()+eps22 ),
        0.,
        1./( ph_p2.pt2()+eps22 )-1./( ph_m2.pt2()+eps22 )
      );

      const double dot1 = phi1.threeProduct( q1 )/qt1_, cross1 = phi1.crossProduct( q1 )/qt1_;
      const double dot2 = phi2.threeProduct( q2 )/qt2_, cross2 = phi2.crossProduct( q2 )/qt2_;
      CG_DEBUG_LOOP( "PPtoFF:offShell" )
        << "phi1 = " << phi1 << "\n\t"
        << "phi2 = " << phi2 << "\n\t"
        << "(dot):   " << dot1 << " / " << dot2 << "\n\t"
        << "(cross): " << cross1 << " / " << cross2 << ".";

      const double aux2_1 = p_term_ll_ * ( mf2_ + 4.*z1*z1*t1abs ) * phi1.energy2()
                           +p_term_tt1_* ( ( z1p*z1p + z1m*z1m )*( dot1*dot1 + cross1*cross1 ) )
                           +p_term_tt2_* ( cross1*cross1 - dot1*dot1 )
                           -p_term_lt_ * 4.*z1*( z1p-z1m ) * phi1.energy() * q1.threeProduct( phi1 );

      const double aux2_2 = p_term_ll_ * ( mf2_ + 4.*z2*z2*t2abs ) * phi2.energy2()
                           +p_term_tt1_* ( ( z2p*z2p + z2m*z2m )*( dot2*dot2 + cross2*cross2 ) )
                           +p_term_tt2_* ( cross2*cross2 - dot2*dot2 )
                           -p_term_lt_ * 4.*z2*( z2p-z2m ) * phi2.energy() * q2.threeProduct( phi2 );

      //=================================================================
      //     convention of matrix element as in our kt-factorization
      //     for heavy flavours
      //=================================================================

      const double amat2_1 = aux2_1*2.*z1*q1.pt2()/( q1.pt2()*q2.pt2() ),
                   amat2_2 = aux2_2*2.*z2*q2.pt2()/( q1.pt2()*q2.pt2() );

      //=================================================================
      //     symmetrization
      //=================================================================

      const double amat2 = 0.5*( p_mat1_*amat2_1 + p_mat2_*amat2_2 );
      CG_DEBUG_LOOP( "PPtoFF:offShell" )
        << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
        << "amat2(1/2), amat2 = " << amat2_1 << " / " << amat2_2 << " / " << amat2 << ".";

      return amat2;
    }
  }
}
// register process and define aliases
REGISTER_PROCESS( pptoll, PPtoFF )
REGISTER_PROCESS( pptoff, PPtoFF )
REGISTER_PROCESS( pptoqq, PPtoFF )
