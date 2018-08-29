#include "CepGen/Processes/GamGamLL.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

namespace CepGen
{
  namespace Process
  {
    //---------------------------------------------------------------------------------------------

    GamGamLL::Masses::Masses() :
      MX2_( 0. ), MY2_( 0. ), Ml2_( 0. ),
      w12_( 0. ), w31_( 0. ), dw31_( 0. ), w52_( 0. ), dw52_( 0. )
    {}

    //---------------------------------------------------------------------------------------------

    GamGamLL::GamGamLL( const ParametersList& params ) :
      GenericProcess( "lpair", "pp → p(*) ( ɣɣ → l⁺l¯ ) p(*)" ),
      n_opt_( params.get<int>( "nopt", 0 ) ),
      pair_( params.get<int>( "pair", 0 ) ),
      ep1_( 0. ), ep2_( 0. ), p_cm_( 0. ),
      ec4_( 0. ), pc4_( 0. ), mc4_( 0. ), w4_( 0. ),
      p12_( 0. ), p1k2_( 0. ), p2k1_( 0. ),
      p13_( 0. ), p14_( 0. ), p25_( 0. ),
      q1dq_( 0. ), q1dq2_( 0. ),
      s1_( 0. ), s2_( 0. ),
      epsi_( 0. ),
      g5_( 0. ), g6_( 0. ), a5_( 0. ), a6_( 0. ), bb_( 0. ),
      gram_( 0. ),
      dd1_( 0. ), dd2_( 0. ), dd3_( 0. ), dd4_( 0. ), dd5_( 0. ),
      delta_( 0. ),
      g4_( 0. ), sa1_( 0. ), sa2_( 0. ),
      sl1_( 0. ),
      cos_theta4_( 0. ), sin_theta4_( 0. ),
      al4_( 0. ), be4_( 0. ), de3_( 0. ), de5_( 0. ),
      pt4_( 0. ),
      jacobian_( 0. )
    {}

    //---------------------------------------------------------------------------------------------

    void
    GamGamLL::addEventContent()
    {
      GenericProcess::setEventContent( {
        { Particle::IncomingBeam1, PDG::Proton },
        { Particle::IncomingBeam2, PDG::Proton },
        { Particle::Parton1, PDG::Photon },
        { Particle::Parton2, PDG::Photon }
      }, {
        { Particle::OutgoingBeam1, { PDG::Proton } },
        { Particle::OutgoingBeam2, { PDG::Proton } },
        { Particle::CentralSystem, { (PDG)pair_, (PDG)pair_ } }
      } );
    }

    unsigned int
    GamGamLL::numDimensions() const
    {
      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElectronProton: default:
          throw CG_FATAL( "GamGamLL" )
            << "Process mode " << cuts_.mode << " not (yet) supported! "
            << "Please contact the developers to consider an implementation.";
        case Kinematics::Mode::ElasticElastic:
          return 7;
        case Kinematics::Mode::ElasticInelastic:
        case Kinematics::Mode::InelasticElastic:
          return 8;
        case Kinematics::Mode::InelasticInelastic:
          return 9;
      }
    }

    //---------------------------------------------------------------------------------------------

    void
    GamGamLL::setKinematics( const Kinematics& kin )
    {
      GenericProcess::setKinematics( kin );

      masses_.Ml2_ = event_->getByRole( Particle::CentralSystem )[0].mass2();

      w_limits_ = cuts_.cuts.central.mass_single;
      if ( !w_limits_.hasMax() )
        w_limits_.max() = s_;
      // The minimal energy for the central system is its outgoing leptons' mass energy (or wmin_ if specified)
      if ( !w_limits_.hasMin() )
        w_limits_.min() = 4.*masses_.Ml2_;
      // The maximal energy for the central system is its CM energy with the outgoing particles' mass energy substracted (or wmax if specified)
      w_limits_.max() = std::min( pow( sqs_-MX_-MY_, 2 ), w_limits_.max() );

      CG_DEBUG_LOOP( "GamGamLL:setKinematics" )
        << "w limits = " << w_limits_ << "\n\t"
        << "wmax/wmin = " << w_limits_.max()/w_limits_.min();

      q2_limits_ = cuts_.cuts.initial.q2;
      mx_limits_ = cuts_.cuts.remnants.mass_single;
    }

    //---------------------------------------------------------------------------------------------

    bool
    GamGamLL::pickin()
    {
      CG_DEBUG_LOOP( "GamGamLL" )
        << "Optimised mode? " << n_opt_;

      jacobian_ = 0.;

      w4_ = mc4_*mc4_;

      // sig1 = sigma and sig2 = sigma' in [1]
      const double sig = mc4_+MY_;
      double sig1 = sig*sig,
             sig2 = sig1;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "mc4 = " << mc4_ << "\n\t"
        << "sig1 = " << sig1 << "\n\t"
        << "sig2 = " << sig2;

      // Mass difference between the first outgoing particle
      // and the first incoming particle
      masses_.w31_ = masses_.MX2_-w1_;
      // Mass difference between the second outgoing particle
      // and the second incoming particle
      masses_.w52_ = masses_.MY2_-w2_;
      // Mass difference between the two incoming particles
      masses_.w12_ = w1_-w2_;
      // Mass difference between the central two-photons system
      // and the second outgoing particle
      const double d6 = w4_-masses_.MY2_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "w1 = " << w1_ << "\n\t"
        << "w2 = " << w2_ << "\n\t"
        << "w3 = " << masses_.MX2_ << "\n\t"
        << "w4 = " << w4_ << "\n\t"
        << "w5 = " << masses_.MY2_;;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "w31 = " << masses_.w31_ << "\n\t"
        << "w52 = " << masses_.w52_ << "\n\t"
        << "w12 = " << masses_.w12_;;

      const double ss = s_+masses_.w12_;

      const double rl1 = ss*ss-4.*w1_*s_; // lambda(s, m1**2, m2**2)
      if ( rl1 <= 0. ) {
        CG_WARNING( "GamGamLL" ) << "rl1 = " << rl1 << " <= 0";
        return false;
      }
      sl1_ = sqrt( rl1 );

      s2_ = 0.;
      double ds2 = 0.;
      if ( n_opt_ == 0 ) {
        const double smax = s_+masses_.MX2_-2.*MX_*sqs_;
        map( x(2), Limits( sig1, smax ), s2_, ds2, "s2" );
        sig1 = s2_; //FIXME!!!!!!!!!!!!!!!!!!!!
      }

      CG_DEBUG_LOOP( "GamGamLL" )
        << "s2 = " << s2_;

      const double sp = s_+masses_.MX2_-sig1,
                   d3 = sig1-w2_;

      const double rl2 = sp*sp-4.*s_*masses_.MX2_; // lambda(s, m3**2, sigma)
      if ( rl2 <= 0. ) {
        CG_DEBUG( "GamGamLL" ) << "rl2 = " << rl2 << " <= 0";
        return false;
      }
      const double sl2 = sqrt( rl2 );

      double t1_max = w1_+masses_.MX2_-( ss*sp+sl1_*sl2 )/( 2.*s_ ), // definition from eq. (A.4) in [1]
             t1_min = ( masses_.w31_*d3+( d3-masses_.w31_ )*( d3*w1_-masses_.w31_*w2_ )/s_ )/t1_max; // definition from eq. (A.5) in [1]

      // FIXME dropped in CDF version
      if ( t1_max > -q2_limits_.min() ) {
        CG_WARNING( "GamGamLL" ) << "t1max = " << t1_max << " > -q2min = " << ( -q2_limits_.min() );
        return false;
      }
      if ( t1_min < -q2_limits_.max() && q2_limits_.hasMax() ) {
        CG_DEBUG( "GamGamLL" ) << "t1min = " << t1_min << " < -q2max = " << -q2_limits_.max();
        return false;
      }
      if ( t1_max < -q2_limits_.max() && q2_limits_.hasMax() )
        t1_max = -q2_limits_.max();
      if ( t1_min > -q2_limits_.min() && q2_limits_.hasMin() )
        t1_min = -q2_limits_.min();
      /////

      // t1, the first photon propagator, is defined here
      t1_ = 0.;
      double dt1 = 0.;
      map( x(0), Limits( t1_min, t1_max ), t1_, dt1, "t1" );
      // changes wrt mapt1 : dx->-dx
      dt1 *= -1.;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Definition of t1 = " << t1_ << " according to\n\t"
        << "(t1min, t1max) = (" << t1_min << ", " << t1_max << ")";

      dd4_ = w4_-t1_;

      const double d8 = t1_-w2_, t13 = t1_-w1_-masses_.MX2_;

      sa1_ = -pow( t1_-masses_.w31_, 2 )/4.+w1_*t1_;
      if ( sa1_ >= 0. ) {
        CG_WARNING( "GamGamLL" ) << "sa1_ = " << sa1_ << " >= 0";
        return false;
      }

      const double sl3 = sqrt( -sa1_ );

      // one computes splus and (s2x=s2max)
      double splus, s2max;
      if ( w1_ != 0. ) {
        const double inv_w1 = 1./w1_;
        const double sb = masses_.MX2_ + 0.5 * ( s_*( t1_-masses_.w31_ )+masses_.w12_*t13 )*inv_w1,
                     sd = sl1_*sl3*inv_w1,
                     se =( s_*( t1_*( s_+t13-w2_ )-w2_*masses_.w31_ )+masses_.MX2_*( masses_.w12_*d8+w2_*masses_.MX2_ ) )*inv_w1;

        if ( fabs( ( sb-sd )/sd ) >= 1. ) {
          splus = sb-sd;
          s2max = se/splus;
        }
        else {
          s2max = sb+sd;
          splus = se/s2max;
        }
      }
      else { // 3
        s2max = ( s_*( t1_*( s_+d8-masses_.MX2_ )-w2_*masses_.MX2_ )+w2_*masses_.MX2_*( w2_+masses_.MX2_-t1_ ) )/( ss*t13 );
        splus = sig2;
      }
      // 4
      double s2x = s2max;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "s2x = s2max = " << s2x;

      if ( n_opt_ < 0 ) { // 5
        if ( splus > sig2 ) {
          sig2 = splus;
          CG_DEBUG_LOOP( "GamGamLL" )
            << "sig2 truncated to splus = " << splus;
        }
        if ( n_opt_ < -1 )
          map( x(2), Limits( sig2, s2max ), s2_, ds2, "s2" );
        else
          mapla( t1_, w2_, x(2), sig2, s2max, s2_, ds2 ); // n_opt_==-1
        s2x = s2_;
      }
      else if ( n_opt_ == 0 )
        s2x = s2_; // 6

      CG_DEBUG_LOOP( "GamGamLL" )
        << "s2x = " << s2x;

      // 7
      const double r1 = s2x-d8, r2 = s2x-d6;

      const double rl4 = ( r1*r1-4.*w2_*s2x )*( r2*r2-4.*masses_.MY2_*s2x );
      if ( rl4 <= 0. ) {
        CG_DEBUG_LOOP( "GamGamLL" )
          << "rl4 = " << rl4 << " <= 0";
        return false;
      }
      const double sl4 = sqrt( rl4 );

      // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
      const double t2_max = w2_+masses_.MY2_-( r1*r2+sl4 )/s2x * 0.5,
                   t2_min = ( masses_.w52_*dd4_+( dd4_-masses_.w52_ )*( dd4_*w2_-masses_.w52_*t1_ )/s2x )/t2_max;

      // t2, the second photon propagator, is defined here
      t2_ = 0.;
      double dt2 = 0.;
      map( x(1), Limits( t2_min, t2_max ), t2_, dt2, "t2" );
      // changes wrt mapt2 : dx->-dx

      dt2 *= -1.;

      // \f$\delta_6=m_4^2-m_5^2\f$ as defined in Vermaseren's paper
      const double tau = t1_-t2_,
                   r3 = dd4_-t2_,
                   r4 = masses_.w52_-t2_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "r1 = " << r1 << "\n\t"
        << "r2 = " << r2 << "\n\t"
        << "r3 = " << r3 << "\n\t"
        << "r4 = " << r4;

      const double b = r3*r4-2.*( t1_+w2_ )*t2_;
      const double c = t2_*d6*d8+( d6-d8 )*( d6*w2_-d8*masses_.MY2_ );

      const double t25 = t2_-w2_-masses_.MY2_;

      sa2_ = -0.25 * r4*r4 + w2_*t2_;
      if ( sa2_ >= 0. ) {
        CG_WARNING( "GamGamLL" )
          <<  "sa2_ = " << sa2_ << " >= 0";
        return false;
      }

      const double sl6 = 2.*sqrt( -sa2_ );

      g4_ = -r3*r3/4.+t1_*t2_;
      if ( g4_ >= 0. ) {
        CG_WARNING( "GamGamLL" ) << "g4_ = " << g4_ << " >= 0";
        return false;
      }

      const double sl7 = 2.*sqrt( -g4_ ),
                   sl5 = sl6*sl7;

      double s2p, s2min;
      if ( fabs( ( sl5-b )/sl5 ) >= 1. ) {
        s2p = 0.5 * ( sl5-b )/t2_;
        s2min = c/( t2_*s2p );
      }
      else { // 8
        s2min = 0.5 * ( -sl5-b )/t2_;
        s2p = c/( t2_*s2min );
      }
      // 9
      if ( n_opt_ > 1 )
        map( x( 2 ), Limits( s2min, s2max ), s2_, ds2, "s2" );
      else if ( n_opt_ == 1 )
        mapla( t1_, w2_, x( 2 ), s2min, s2max, s2_, ds2 );

      const double ap = -0.25*pow( s2_+d8, 2 )+s2_*t1_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "s2 = " << s2_ << ", s2max = " << s2max << ", splus = " << splus;

      if ( w1_ != 0. )
        dd1_ = -0.25 * ( s2_-s2max ) * ( s2_-splus ) * w1_; // 10
      else
        dd1_ =  0.25 * ( s2_-s2max ) * ss * t13;
      // 11
      dd2_ = -0.25 * t2_*( s2_-s2p )*( s2_-s2min );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "t2    = " << t2_ << "\n\t"
        << "s2    = " << s2_ << "\n\t"
        << "s2p   = " << s2p << "\n\t"
        << "s2min = " << s2min << "\n\t"
        << "dd2   = " << dd2_;;

      const double yy4 = cos( M_PI*x( 3 ) );
      const double dd = dd1_*dd2_;
      p12_ = 0.5 * ( s_-w1_-w2_ );
      const double st = s2_-t1_-w2_;
      const double delb = ( 2.*w2_*r3+r4*st )*( 4.*p12_*t1_-( t1_-masses_.w31_ )*st )/( 16.*ap );

      if ( dd <= 0. ) {
        CG_DEBUG_LOOP( "GamGamLL" )
          << std::scientific
          << "dd = " << dd << " <= 0\n\t"
          << "dd1 = " << dd1_ << "\t"
          << "dd2 = " << dd2_ << std::fixed;
        return false;
      }

      delta_ = delb - yy4*st*sqrt( dd )/ ap * 0.5;
      s1_ = t2_+w1_+( 2.*p12_*r3-4.*delta_ )/st;

      if ( ap >= 0. ) {
        CG_DEBUG_LOOP( "GamGamLL" )
          <<  "ap = " << ap << " >= 0";
        return false;
      }

      jacobian_ = ds2 * dt1 * dt2 * 0.125 * M_PI*M_PI/( sl1_*sqrt( -ap ) );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Jacobian = " << std::scientific << jacobian_ << std::fixed;

      gram_ = ( 1.-yy4*yy4 )*dd/ap;

      p13_ = -0.5 * t13;
      p14_ =  0.5 * ( tau+s1_-masses_.MX2_ );
      p25_ = -0.5 * t25;

      p1k2_ = 0.5 * ( s1_-t2_-w1_ );
      p2k1_ = 0.5 * st;

      if ( w2_ != 0. ) {
        const double inv_w2 = 1./w2_;
        const double sbb = 0.5 * ( s_*( t2_-masses_.w52_ )-masses_.w12_*t25 )*inv_w2 + masses_.MY2_,
                     sdd = 0.5 * sl1_*sl6*inv_w2,
                     see = ( s_*( t2_*( s_+t25-w1_ )-w1_*masses_.w52_ )+masses_.MY2_*( w1_*masses_.MY2_-masses_.w12_*( t2_-w1_ ) ) )*inv_w2;
        double s1m = 0., s1p = 0.;
        if ( sbb/sdd >= 0. ) {
          s1p = sbb+sdd;
          s1m = see/s1p;
        }
        else {
          s1m = sbb-sdd;
          s1p = see/s1m;
        } // 12
        dd3_ = -0.25 * w2_*( s1p-s1_ )*( s1m-s1_ ); // 13
      }
      else { // 14
        const double s1p = ( s_*( t2_*( s_-masses_.MY2_+t2_-w1_ )-w1_*masses_.MY2_ )+w1_*masses_.MY2_*( w1_+masses_.MY2_-t2_ ) )/( t25*( s_-masses_.w12_ ) );
        dd3_ = -0.25 * t25*( s_-masses_.w12_ )*( s1p-s1_ );
      }
      // 15
      //const double acc3 = (s1p-s1_)/(s1p+s1_);

      const double ssb = t2_+0.5 * w1_-r3*( masses_.w31_-t1_ )/t1_,
                   ssd = sl3*sl7/t1_,
                   sse = ( t2_-w1_ )*( w4_-masses_.MX2_ )+( t2_-w4_+masses_.w31_ )*( ( t2_-w1_)*masses_.MX2_-(w4_-masses_.MX2_)*w1_)/t1_;

      double s1pp, s1pm;
      if ( ssb/ssd >= 0. ) {
        s1pp = ssb+ssd;
        s1pm = sse/s1pp;
      }
      else { // 16
        s1pm = ssb-ssd;
        s1pp = sse/s1pm;
      }
      // 17
      dd4_ = -0.25 * t1_*( s1_-s1pp )*( s1_-s1pm );
      //const double acc4 = ( s1_-s1pm )/( s1_+s1pm );
      dd5_ = dd1_+dd3_+( ( p12_*( t1_-masses_.w31_ )*0.5-w1_*p2k1_ )*( p2k1_*( t2_-masses_.w52_ )-w2_*r3 )
                        -delta_*( 2.*p12_*p2k1_-w2_*( t1_-masses_.w31_ ) ) ) / p2k1_;

      return true;
    }

    //---------------------------------------------------------------------------------------------

    bool
    GamGamLL::orient()
    {
      if ( !pickin() || jacobian_ == 0. ) {
        CG_DEBUG_LOOP( "GamGamLL" )
          << "Pickin failed! Jacobian = " << jacobian_;
        return false;
      }

      const double re = 0.5 / sqs_;
      ep1_ = re*( s_+masses_.w12_ );
      ep2_ = re*( s_-masses_.w12_ );

      CG_DEBUG_LOOP( "GamGamLL" )
        << std::scientific
        << " re = " << re << "\n\t"
        << "w12_ = " << masses_.w12_
        << std::fixed;
      CG_DEBUG_LOOP( "GamGamLL" )
        << "Incoming particles' energy = " << ep1_ << ", " << ep2_;

      p_cm_ = re*sl1_;

      de3_ = re*(s2_-masses_.MX2_+masses_.w12_);
      de5_ = re*(s1_-masses_.MY2_-masses_.w12_);

      // Final state energies
      const double ep3 = ep1_-de3_,
                   ep5 = ep2_-de5_;
      ec4_ = de3_+de5_;

      if ( ec4_ < mc4_ ) {
        CG_WARNING( "GamGamLL" )
          << "ec4_ = " << ec4_ << " < mc4_ = " << mc4_ << "\n\t"
          << "==> de3 = " << de3_ << ", de5 = " << de5_;
        return false;
      }

      // What if the protons' momenta are not along the z-axis?
      pc4_ = sqrt( ec4_*ec4_-mc4_*mc4_ );

      if ( pc4_ == 0. ) {
        CG_WARNING( "GamGamLL" ) << "pzc4 is null and should not be...";
        return false;
      }

      const double pp3 = sqrt( ep3*ep3-masses_.MX2_ ), pt3 = sqrt( dd1_/s_ )/p_cm_,
                   pp5 = sqrt( ep5*ep5-masses_.MY2_ ), pt5 = sqrt( dd3_/s_ )/p_cm_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Central system's energy: E4 = " << ec4_ << "\n\t"
        << "               momentum: p4 = " << pc4_ << "\n\t"
        << "         invariant mass: m4 = " << mc4_ << "\n\t"
        << "Outgoing particles' energy: E3 = " << ep3 << "\n\t"
        << "                            E5 = " << ep5;

      const double sin_theta3 = pt3/pp3,
                   sin_theta5 = pt5/pp5;

      CG_DEBUG_LOOP( "GamGamLL" )
        << std::scientific
        << "sin_theta3 = " << sin_theta3 << "\n\t"
        << "sin_theta5 = " << sin_theta5
        << std::fixed;

      if ( sin_theta3 > 1. ) {
        CG_WARNING( "GamGamLL" )
          << "sin_theta3 = " << sin_theta3 << " > 1";
        return false;
      }
      if ( sin_theta5 > 1. ) {
        CG_WARNING( "GamGamLL" )
          << "sin_theta5 = " << sin_theta5 << " > 1";
        return false;
      }

      double ct3 = sqrt( 1.-sin_theta3*sin_theta3 ),
             ct5 = sqrt( 1.-sin_theta5*sin_theta5 );

      if ( ep1_*ep3 < p13_ )
        ct3 *= -1.;
      if ( ep2_*ep5 > p25_ )
        ct5 *= -1.;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "ct3 = " << ct3 << "\n\t"
        << "ct5 = " << ct5;

      if ( dd5_ < 0. ) {
        CG_WARNING( "GamGamLL" )
          <<  "dd5 = " << dd5_ << " < 0";
        return false;
      }

      // Centre of mass system kinematics (theta4 and phi4)
      pt4_ = sqrt( dd5_/s_ )/p_cm_;
      sin_theta4_ = pt4_/pc4_;

      if ( sin_theta4_ > 1. ) {
        CG_WARNING( "GamGamLL" )
          << "st4 = " << sin_theta4_ << " > 1";
        return false;
      }

      cos_theta4_ = sqrt( 1.-sin_theta4_*sin_theta4_ );
      if ( ep1_*ec4_ < p14_ )
        cos_theta4_ *= -1.;

      al4_ = 1.-cos_theta4_;
      be4_ = 1.+cos_theta4_;

      if ( cos_theta4_ < 0. ) be4_ = sin_theta4_*sin_theta4_/al4_;
      else                    al4_ = sin_theta4_*sin_theta4_/be4_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "ct4 = " << cos_theta4_ << "\n\t"
        << "al4 = " << al4_ << ", be4 = " << be4_;

      const double rr  = sqrt( -gram_/s_ )/( p_cm_*pt4_ );
      const double sin_phi3 =  rr / pt3,
                   sin_phi5 = -rr / pt5;

      if ( fabs( sin_phi3 ) > 1. ) {
        CG_WARNING( "GamGamLL" )
          << "sin(phi_3) = " << sin_phi3 << " while it must be in [-1 ; 1]";
        return false;
      }
      if ( fabs( sin_phi5 ) > 1. ) {
        CG_WARNING( "GamGamLL" )
          << "sin(phi_5) = " << sin_phi5 << " while it must be in [-1 ; 1]";
        return false;
      }

      const double cos_phi3 = -sqrt( 1.-sin_phi3*sin_phi3 ),
                   cos_phi5 = -sqrt( 1.-sin_phi5*sin_phi5 );

      p3_lab_ = Particle::Momentum( pp3*sin_theta3*cos_phi3, pp3*sin_theta3*sin_phi3, pp3*ct3, ep3 );
      p5_lab_ = Particle::Momentum( pp5*sin_theta5*cos_phi5, pp5*sin_theta5*sin_phi5, pp5*ct5, ep5 );

      const double a1 = p3_lab_.px()-p5_lab_.px();

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Kinematic quantities\n\t"
        << "cos(theta3) = " << ct3 << "\t" << "sin(theta3) = " << sin_theta3 << "\n\t"
        << "cos( phi3 ) = " << cos_phi3 << "\t" << "sin( phi3 ) = " << sin_phi3 << "\n\t"
        << "cos(theta4) = " << cos_theta4_ << "\t" << "sin(theta4) = " << sin_theta4_ << "\n\t"
        << "cos(theta5) = " << ct5 << "\t" << "sin(theta5) = " << sin_theta5 << "\n\t"
        << "cos( phi5 ) = " << cos_phi5 << "\t" << "sin( phi5 ) = " << sin_phi5 << "\n\t"
        << "a1 = " << a1;

      if ( fabs( pt4_+p3_lab_.px()+p5_lab_.px() ) < fabs( fabs( a1 )-pt4_ ) ) {
        CG_DEBUG_LOOP( "GamGamLL" )
          << "|pt4+pt3*cos(phi3)+pt5*cos(phi5)| < | |a1|-pt4 |\n\t"
          << "pt4 = " << pt4_ << "\t"
          << "pt5 = " << pt5 << "\n\t"
          << "cos(phi3) = " << cos_phi3 << "\t"
          << "cos(phi5) = " << cos_phi5 << "\n\t"
          << "a1 = " << a1;
        return true;
      }
      if ( a1 < 0. ) p5_lab_[0] = -p5_lab_.px();
      else           p3_lab_[0] = -p3_lab_.px();
      return true;
    }

    //---------------------------------------------------------------------------------------------

    double
    GamGamLL::computeOutgoingPrimaryParticlesMasses( double x, double outmass, double lepmass, double& dw )
    {
      const double mx0 = mp_+ParticleProperties::mass( PDG::PiZero ); // 1.07
      const double wx2min = pow( std::max( mx0, mx_limits_.min() ), 2 ),
                   wx2max = pow( std::min( sqs_-outmass-2.*lepmass, mx_limits_.max() ), 2 );

      double mx2 = 0., dmx2 = 0.;
      map( x, Limits( wx2min, wx2max ), mx2, dmx2, "mx2" );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "mX^2 in range (" << wx2min << ", " << wx2max << "), x = " << x << "\n\t"
        << "mX^2 = " << mx2 << ", d(mX^2) = " << dmx2 << "\n\t"
        << "mX = " << sqrt( mx2 ) << ", d(mX) = " << sqrt( dmx2 );

      dw = sqrt( dmx2 );
      return sqrt( mx2 );
    }

    //---------------------------------------------------------------------------------------------

    void
    GamGamLL::beforeComputeWeight()
    {
      if ( !GenericProcess::is_point_set_ ) return;

      const Particle& p1 = event_->getOneByRole( Particle::IncomingBeam1 ),
                     &p2 = event_->getOneByRole( Particle::IncomingBeam2 );

      ep1_ = p1.energy();
      ep2_ = p2.energy();

      //std::cout << __PRETTY_FUNCTION__ << ":" << w_limits_ << "|" << q2_limits_ << "|" << mx_limits_ << std::endl;
      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElectronProton: default:
          CG_ERROR( "GamGamLL" ) << "Case not yet supported!"; break;
        case Kinematics::Mode::ElasticElastic:
          masses_.dw31_ = masses_.dw52_ = 0.; break;
        case Kinematics::Mode::InelasticElastic: {
          const double m = computeOutgoingPrimaryParticlesMasses( x( 7 ), p1.mass(), sqrt( masses_.Ml2_ ), masses_.dw31_ );
          event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( m );
          event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( ParticleProperties::mass( p2.pdgId() ) );
        } break;
        case Kinematics::Mode::ElasticInelastic: {
          const double m = computeOutgoingPrimaryParticlesMasses( x( 7 ), p2.mass(), sqrt( masses_.Ml2_ ), masses_.dw52_ );
          event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( ParticleProperties::mass( p1.pdgId() ) );
          event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( m );
        } break;
        case Kinematics::Mode::InelasticInelastic: {
          const double mx = computeOutgoingPrimaryParticlesMasses( x( 7 ), p2.mass(), sqrt( masses_.Ml2_ ), masses_.dw31_ );
          event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( mx );
          const double my = computeOutgoingPrimaryParticlesMasses( x( 8 ), p1.mass(), sqrt( masses_.Ml2_ ), masses_.dw52_ );
          event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( my );
        } break;
      }
      MX_ = event_->getOneByRole( Particle::OutgoingBeam1 ).mass();
      MY_ = event_->getOneByRole( Particle::OutgoingBeam2 ).mass();
      masses_.MX2_ = MX_*MX_;
      masses_.MY2_ = MY_*MY_;
    }

    //---------------------------------------------------------------------------------------------

    double
    GamGamLL::computeWeight()
    {
      CG_DEBUG_LOOP( "GamGamLL" )
        << "sqrt(s) = " << sqs_ << " GeV\n\t"
        << "m(X1) = " << MX_ << " GeV\t"
        << "m(X2) = " << MY_ << " GeV";

      // compute the two-photon energy for this point
      w4_ = 0.;
      double dw4 = 0.;
      map( x( 4 ), w_limits_, w4_, dw4, "w4" );
      mc4_ = sqrt( w4_ );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Computed value for w4 = " << w4_ << " → mc4 = " << mc4_;

      if ( !orient() )
        return 0.;

      if ( jacobian_ == 0. ) {
        CG_WARNING( "GamGamLL" ) << "dj = " << jacobian_;
        return 0.;
      }

      if ( t1_ > 0. ) {
        CG_WARNING( "GamGamLL" ) << "t1 = " << t1_ << " > 0";
        return 0.;
      }
      if ( t2_ > 0. ) {
        CG_WARNING( "GamGamLL" ) << "t2 = " << t2_ << " > 0";
        return 0.;
      }

      const double ecm6 = w4_ / ( 2.*mc4_ ),
                   pp6cm = sqrt( ecm6*ecm6-masses_.Ml2_ );

      jacobian_ *= dw4*pp6cm/( mc4_*Constants::sconstb*s_ );

      // Let the most obscure part of this code begin...

      const double e1mp1 = w1_ / ( ep1_+p_cm_ ),
                   e3mp3 = masses_.MX2_ / ( p3_lab_.energy()+p3_lab_.p() );

      const double al3 = pow( sin( p3_lab_.theta() ), 2 )/( 1.+( p3_lab_.theta() ) );

      // 2-photon system kinematics ?!
      const double eg = ( w4_+t1_-t2_ )/( 2.*mc4_ );
      double pg = sqrt( eg*eg-t1_ );

      const double pgx = -p3_lab_.px()*cos_theta4_-sin_theta4_*( de3_-e1mp1 + e3mp3 + p3_lab_.p()*al3 ),
                   pgy = -p3_lab_.py(),
                   pgz = mc4_*de3_/( ec4_+pc4_ )-ec4_*de3_*al4_/mc4_-p3_lab_.px()*ec4_*sin_theta4_/mc4_+ec4_*cos_theta4_/mc4_*( p3_lab_.p()*al3+e3mp3-e1mp1 );

      CG_DEBUG_LOOP( "GamGamLL" ) << "pg = " << Particle::Momentum( pgx, pgy, pgz );

      const double pgp = sqrt( pgx*pgx + pgy*pgy ), // outgoing proton (3)'s transverse momentum
                   pgg = sqrt( pgp*pgp + pgz*pgz ); // outgoing proton (3)'s momentum
      if ( pgg > pgp*0.9 && pgg > pg )
        pg = pgg; //FIXME ???

      // angles for the 2-photon system ?!
      const double cpg = pgx/pgp, spg = pgy/pgp;
      const double stg = pgp/pg;

      const int theta_sign = ( pgz>0. ) ? 1 : -1;
      const double ctg = theta_sign*sqrt( 1.-stg*stg );

      double xx6 = x( 5 );

      const double amap = 0.5 * ( w4_-t1_-t2_ ),
                   bmap = 0.5 * sqrt( ( pow( w4_-t1_-t2_, 2 )-4.*t1_*t2_ )*( 1.-4.*masses_.Ml2_/w4_ ) ),
                   ymap = ( amap+bmap )/( amap-bmap ),
                   beta = pow( ymap, 2.*xx6-1. );
      xx6 = 0.5 * ( 1. + amap/bmap*( beta-1. )/( beta+1. ) );
      xx6 = std::max( 0., std::min( xx6, 1. ) ); // xx6 in [0., 1.]

      CG_DEBUG_LOOP( "GamGamLL" )
        << "amap = " << amap << "\n\t"
        << "bmap = " << bmap << "\n\t"
        << "ymap = " << ymap << "\n\t"
        << "beta = " << beta;

      // 3D rotation of the first outgoing lepton wrt the CM system
      const double theta6cm = acos( 1.-2.*xx6 );

      // match the Jacobian
      jacobian_ *= ( amap+bmap*cos( theta6cm ) );
      jacobian_ *= ( amap-bmap*cos( theta6cm ) );
      jacobian_ /= amap;
      jacobian_ /= bmap;
      jacobian_ *= log( ymap );
      jacobian_ *= 0.5;

      CG_DEBUG_LOOP( "GamGamLL" ) << "Jacobian = " << jacobian_;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "ctcm6 = " << cos( theta6cm ) << "\n\t"
        << "stcm6 = " << sin( theta6cm );

      const double phi6cm = 2.*M_PI*x( 6 );

      // First outgoing lepton's 3-momentum in the centre of mass system
      Particle::Momentum p6cm = Particle::Momentum::fromPThetaPhi( pp6cm, theta6cm, phi6cm );

      CG_DEBUG_LOOP( "GamGamLL" ) << "p3cm6 = " << p6cm;

      const double h1 = stg*p6cm.pz()+ctg*p6cm.px();
      const double pc6z = ctg*p6cm.pz()-stg*p6cm.px(), pc6x = cpg*h1-spg*p6cm.py();

      const double qcx = 2.*pc6x, qcz = 2.*pc6z;
      // qcy == QCY is never defined

      const double el6 = ( ec4_*ecm6+pc4_*pc6z ) / mc4_,
                   h2  = ( ec4_*pc6z+pc4_*ecm6 ) / mc4_;

      CG_DEBUG_LOOP( "GamGamLL" ) << "h1 = " << h1 << "\n\th2 = " << h2;

      // first outgoing lepton's 3-momentum
      const double p6x = cos_theta4_*pc6x+sin_theta4_*h2,
                   p6y = cpg*p6cm.py()+spg*h1,
                   p6z = cos_theta4_*h2-sin_theta4_*pc6x;

      // first outgoing lepton's kinematics
      p6_cm_ = Particle::Momentum( p6x, p6y, p6z, el6 );
      CG_DEBUG_LOOP( "GamGamLL" ) << "p6(cm) = " << p6_cm_;

      const double hq = ec4_*qcz/mc4_;

      const Particle::Momentum qve(
        cos_theta4_*qcx+sin_theta4_*hq,
        2.*p6y,
        cos_theta4_*hq-sin_theta4_*qcx,
        pc4_*qcz/mc4_ // energy
      );

      // Available energy for the second lepton is the 2-photon system's energy with the first lepton's energy removed
      const double el7 = ec4_-el6;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "Outgoing kinematics\n\t"
        << " first outgoing lepton: p = " << p6_cm_.p() << ", E = " << p6_cm_.energy() << "\n\t"
        << "second outgoing lepton: p = " << p7_cm_.p() << ", E = " << p7_cm_.energy();;

      // Second outgoing lepton's 3-momentum
      const double p7x = -p6x + pt4_,
                   p7y = -p6y,
                   p7z = -p6z + pc4_*cos_theta4_;

      // second outgoing lepton's kinematics
      p7_cm_ = Particle::Momentum( p7x, p7y, p7z, el7 );

      //p6_cm_ = Particle::Momentum(pl6*st6*cp6, pl6*st6*sp6, pl6*ct6, el6);
      //p7_cm_ = Particle::Momentum(pl7*st7*cp7, pl7*st7*sp7, pl7*ct7, el7);

      q1dq_ = eg*( 2.*ecm6-mc4_ )-2.*pg*p6cm.pz();
      q1dq2_ = ( w4_-t1_-t2_ ) * 0.5;

      const double phi3 = p3_lab_.phi(), cos_phi3 = cos( phi3 ), sin_phi3 = sin( phi3 ),
                   phi5 = p5_lab_.phi(), cos_phi5 = cos( phi5 ), sin_phi5 = sin( phi5 );

      bb_ = t1_*t2_+( w4_*pow( sin( theta6cm ), 2 ) + 4.*masses_.Ml2_*pow( cos( theta6cm ), 2 ) )*pg*pg;

      const double c1 = p3_lab_.pt() * ( qve.px()*sin_phi3  - qve.py()*cos_phi3   ),
                   c2 = p3_lab_.pt() * ( qve.pz()*ep1_ - qve.energy() *p_cm_ ),
                   c3 = ( masses_.w31_*ep1_*ep1_ + 2.*w1_*de3_*ep1_ - w1_*de3_*de3_ + p3_lab_.pt2()*ep1_*ep1_ ) / ( p3_lab_.energy()*p_cm_ + p3_lab_.pz()*ep1_ );

      const double b1 = p5_lab_.pt() * ( qve.px()*sin_phi5  - qve.py()*cos_phi5   ),
                   b2 = p5_lab_.pt() * ( qve.pz()*ep2_ + qve.energy() *p_cm_ ),
                   b3 = ( masses_.w52_*ep2_*ep2_ + 2.*w2_*de5_*ep2_ - w2_*de5_*de5_ + p5_lab_.pt2()*ep2_*ep2_ ) / ( ep2_*p5_lab_.pz() - p5_lab_.energy()*p_cm_ );

      const double r12 =  c2*sin_phi3 + qve.py()*c3,
                   r13 = -c2*cos_phi3 - qve.px()*c3;

      const double r22 =  b2*sin_phi5 + qve.py()*b3,
                   r23 = -b2*cos_phi5 - qve.px()*b3;

      epsi_ = p12_*c1*b1 + r12*r22 + r13*r23;

      g5_ = w1_*c1*c1 + r12*r12 + r13*r13;
      g6_ = w2_*b1*b1 + r22*r22 + r23*r23;

      const double pt3 = p3_lab_.pt(), pt5 = p5_lab_.pt();
      a5_ = -( qve.px()*cos_phi3 + qve.py()*sin_phi3 )*pt3*p1k2_
            -( ep1_*qve.energy()-p_cm_*qve.pz() )*( cos_phi3*cos_phi5 + sin_phi3*sin_phi5 )*pt3*pt5
            +( de5_*qve.pz()+qve.energy()*( p_cm_+p5_lab_.pz() ) )*c3;
      a6_ = -( qve.px()*cos_phi5 + qve.py()*sin_phi5 )*pt5*p2k1_
            -( ep2_*qve.energy()+p_cm_*qve.pz() )*( cos_phi3*cos_phi5 + sin_phi3*sin_phi5 )*pt3*pt5
            +( de3_*qve.pz()-qve.energy()*( p_cm_-p3_lab_.pz() ) )*b3;

      CG_DEBUG_LOOP( "GamGamLL" )
        << "a5 = " << a5_ << "\n\t"
        << "a6 = " << a6_;

      ////////////////////////////////////////////////////////////////
      // END of GAMGAMLL subroutine in the FORTRAN version
      ////////////////////////////////////////////////////////////////

      const Particle::Momentum cm = event_->getOneByRole( Particle::IncomingBeam1 ).momentum()
                                  + event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

      ////////////////////////////////////////////////////////////////
      // INFO from f.f
      ////////////////////////////////////////////////////////////////

      const double gamma = cm.energy() / sqs_, betgam = cm.pz() / sqs_;

      //--- kinematics computation for both leptons

      p6_cm_.betaGammaBoost( gamma, betgam );
      p7_cm_.betaGammaBoost( gamma, betgam );

      //--- cut on mass of final hadronic system (MX/Y)

      if ( mx_limits_.valid() ) {
        if ( ( cuts_.mode == Kinematics::Mode::InelasticElastic
            || cuts_.mode == Kinematics::Mode::InelasticInelastic )
          && !mx_limits_.passes( MX_ ) )
          return 0.;
        if ( ( cuts_.mode == Kinematics::Mode::ElasticInelastic
            || cuts_.mode == Kinematics::Mode::InelasticInelastic )
          && !mx_limits_.passes( MY_ ) )
          return 0.;
      }

      //--- cut on the proton's Q2 (first photon propagator T1)

      if ( !cuts_.cuts.initial.q2.passes( -t1_ ) )
        return 0.;

      //--- cuts on outgoing leptons' kinematics

      if ( !cuts_.cuts.central.mass_sum.passes( ( p6_cm_+p7_cm_ ).mass() ) )
        return 0.;

      //----- cuts on the individual leptons

      if ( cuts_.cuts.central.pt_single.valid() ) {
        const Limits& pt_limits = cuts_.cuts.central.pt_single;
        if ( !pt_limits.passes( p6_cm_.pt() ) || !pt_limits.passes( p7_cm_.pt() ) )
          return 0.;
      }

      if ( cuts_.cuts.central.energy_single.valid() ) {
        const Limits& energy_limits = cuts_.cuts.central.energy_single;
        if ( !energy_limits.passes( p6_cm_.energy() ) || !energy_limits.passes( p7_cm_.energy() ) )
          return 0.;
      }

      if ( cuts_.cuts.central.eta_single.valid() ) {
        const Limits& eta_limits = cuts_.cuts.central.eta_single;
        if ( !eta_limits.passes( p6_cm_.eta() ) || !eta_limits.passes( p7_cm_.eta() ) )
          return 0.;
      }

      //--- compute the structure functions factors

      switch ( cuts_.mode ) { // inherited from CDF version
        case Kinematics::Mode::ElectronProton: default: jacobian_ *= periPP( 1, 2 ); break; // ep case
        case Kinematics::Mode::ElasticElastic:          jacobian_ *= periPP( 2, 2 ); break; // elastic case
        case Kinematics::Mode::InelasticElastic:        jacobian_ *= periPP( 3, 2 )*( masses_.dw31_*masses_.dw31_ ); break;
        case Kinematics::Mode::ElasticInelastic:        jacobian_ *= periPP( 3, 2 )*( masses_.dw52_*masses_.dw52_ ); break; // single-dissociative case
        case Kinematics::Mode::InelasticInelastic:      jacobian_ *= periPP( 3, 3 )*( masses_.dw31_*masses_.dw31_ )*( masses_.dw52_*masses_.dw52_ ); break; // double-dissociative case
      }

      //--- compute the event weight using the Jacobian

      return Constants::GeV2toBarn*jacobian_;
    }

    //---------------------------------------------------------------------------------------------

    void
    GamGamLL::fillKinematics( bool )
    {
      const Particle::Momentum cm = event_->getOneByRole( Particle::IncomingBeam1 ).momentum() + event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

      const double gamma  = cm.energy()/sqs_, betgam = cm.pz()/sqs_;

      Particle::Momentum plab_ip1 = Particle::Momentum( 0., 0.,  p_cm_, ep1_ ).betaGammaBoost( gamma, betgam );
      Particle::Momentum plab_ip2 = Particle::Momentum( 0., 0., -p_cm_, ep2_ ).betaGammaBoost( gamma, betgam );
      p3_lab_.betaGammaBoost( gamma, betgam );
      p5_lab_.betaGammaBoost( gamma, betgam );

      //----- parameterise a random rotation around z-axis
      const int rany = ( rand() >= 0.5 * RAND_MAX ) ? 1 : -1,
                ransign = ( rand() >= 0.5 * RAND_MAX ) ? 1 : -1;
      const double ranphi = ( 0.5 * rand()/RAND_MAX )*2.*M_PI;

      Particle::Momentum plab_ph1 = ( plab_ip1-p3_lab_ ).rotatePhi( ranphi, rany );
      Particle::Momentum plab_ph2 = ( plab_ip2-p5_lab_ ).rotatePhi( ranphi, rany );

      p3_lab_.rotatePhi( ranphi, rany );
      p5_lab_.rotatePhi( ranphi, rany );
      p6_cm_.rotatePhi( ranphi, rany );
      p7_cm_.rotatePhi( ranphi, rany );

      /*if ( symmetrise_ && rand() >= .5*RAND_MAX ) {
        p6_cm_.mirrorZ();
        p7_cm_.mirrorZ();
      }*/

      //----- incoming protons
      event_->getOneByRole( Particle::IncomingBeam1 ).setMomentum( plab_ip1 );
      event_->getOneByRole( Particle::IncomingBeam2 ).setMomentum( plab_ip2 );

      //----- first outgoing proton
      Particle& op1 = event_->getOneByRole( Particle::OutgoingBeam1 );

      op1.setMomentum( p3_lab_ );
      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElasticElastic:
        case Kinematics::Mode::ElasticInelastic:
        default:
          op1.setStatus( Particle::Status::FinalState ); // stable proton
          break;
        case Kinematics::Mode::InelasticElastic:
        case Kinematics::Mode::InelasticInelastic:
          op1.setStatus( Particle::Status::Unfragmented ); // fragmenting remnants
          op1.setMass( MX_ );
          break;
      }

      //----- second outgoing proton
      Particle& op2 = event_->getOneByRole( Particle::OutgoingBeam2 );
      op2.setMomentum( p5_lab_ );
      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElasticElastic:
        case Kinematics::Mode::InelasticElastic:
        default:
          op2.setStatus( Particle::Status::FinalState ); // stable proton
          break;
        case Kinematics::Mode::ElasticInelastic:
        case Kinematics::Mode::InelasticInelastic:
          op2.setStatus( Particle::Status::Unfragmented ); // fragmenting remnants
          op2.setMass( MY_ );
          break;
      }

      //----- first incoming photon
      Particle& ph1 = event_->getOneByRole( Particle::Parton1 );
      ph1.setMomentum( plab_ph1 );

      //----- second incoming photon
      Particle& ph2 = event_->getOneByRole( Particle::Parton2 );
      ph2.setMomentum( plab_ph2 );

      Particles& central_system = event_->getByRole( Particle::CentralSystem );

      //----- first outgoing lepton
      Particle& ol1 = central_system[0];
      ol1.setPdgId( ol1.pdgId(), ransign );
      ol1.setMomentum( p6_cm_ );
      ol1.setStatus( Particle::Status::FinalState );

      //----- second outgoing lepton
      Particle& ol2 = central_system[1];
      ol2.setPdgId( ol2.pdgId(), -ransign );
      ol2.setMomentum( p7_cm_ );
      ol2.setStatus( Particle::Status::FinalState );

      //----- intermediate two-lepton system
      event_->getOneByRole( Particle::Intermediate ).setMomentum( p6_cm_+p7_cm_ );
    }

    //---------------------------------------------------------------------------------------------

    double
    GamGamLL::periPP( int nup_, int ndown_ )
    {
      CG_DEBUG_LOOP( "GamGamLL" )
        << " Nup  = " << nup_ << "\n\t"
        << "Ndown = " << ndown_;

      FormFactors fp1, fp2;
      formFactors( -t1_, -t2_, fp1, fp2 );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "u1 = " << fp1.FM << "\n\t"
        << "u2 = " << fp1.FE << "\n\t"
        << "v1 = " << fp2.FM << "\n\t"
        << "v2 = " << fp2.FE;

      const double qqq = q1dq_*q1dq_,
                   qdq = 4.*masses_.Ml2_-w4_;
      const double t11 = 64. *(  bb_*( qqq-g4_-qdq*( t1_+t2_+2.*masses_.Ml2_ ) )-2.*( t1_+2.*masses_.Ml2_ )*( t2_+2.*masses_.Ml2_ )*qqq ) * t1_*t2_, // magnetic-magnetic
                   t12 = 128.*( -bb_*( dd2_+g6_ )-2.*( t1_+2.*masses_.Ml2_ )*( sa2_*qqq+a6_*a6_ ) ) * t1_, // electric-magnetic
                   t21 = 128.*( -bb_*( dd4_+g5_ )-2.*( t2_+2.*masses_.Ml2_ )*( sa1_*qqq+a5_*a5_ ) ) * t2_, // magnetic-electric
                   t22 = 512.*(  bb_*( delta_*delta_-gram_ )-pow( epsi_-delta_*( qdq+q1dq2_ ), 2 )-sa1_*a6_*a6_-sa2_*a5_*a5_-sa1_*sa2_*qqq ); // electric-electric

      const double peripp = ( fp1.FM*fp2.FM*t11 + fp1.FE*fp2.FM*t21 + fp1.FM*fp2.FE*t12 + fp1.FE*fp2.FE*t22 ) / pow( 2.*t1_*t2_*bb_, 2 );

      CG_DEBUG_LOOP( "GamGamLL" )
        << "t11 = " << t11 << "\t"
        << "t12 = " << t12 << "\n\t"
        << "t21 = " << t21 << "\t"
        << "t22 = " << t22 << "\n\t"
        << "=> PeriPP = " << peripp;

      return peripp;
    }

    void
    GamGamLL::map( double expo, const Limits& lim, double& out, double& dout, const std::string& var_name_ )
    {
      const double y = lim.max()/lim.min();
      out = lim.min()*pow( y, expo );
      dout = out*log( y );
      CG_DEBUG_LOOP( "GamGamLL:map" )
        << "Mapping variable \"" << var_name_ << "\"\n\t"
        << "limits = " << lim << "\n\t"
        << "max/min = " << y << "\n\t"
        << "exponent = " << expo << "\n\t"
        << "output = " << out << "\n\t"
        << "d(output) = " << dout;
    }

    void
    GamGamLL::mapla( double y, double z, int u, double xm, double xp, double& x, double& d )
    {
      const double xmb = xm-y-z, xpb = xp-y-z;
      const double c = -4.*y*z;
      const double alp = sqrt( xpb*xpb + c ), alm = sqrt( xmb*xmb + c );
      const double am = xmb+alm, ap = xpb+alp;
      const double yy = ap/am, zz = pow( yy, u );

      x = y + z + ( am*zz - c / ( am*zz ) ) / 2.;
      const double ax = sqrt( pow( x-y-z, 2 ) + c );
      d = ax*log( yy );
    }

    void
    GamGamLL::formFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const
    {
      const double mx2 = MX_*MX_, my2 = MY_*MY_;

      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElectronElectron: default: {
          fp1 = FormFactors::trivial(); // electron (trivial) form factor
          fp2 = FormFactors::trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::Mode::ProtonElectron: {
          fp1 = FormFactors::protonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::Mode::ElectronProton: {
          fp1 = FormFactors::trivial(); // electron (trivial) form factor
          fp2 = FormFactors::protonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::Mode::ElasticElastic: {
          fp1 = FormFactors::protonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::protonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::Mode::ElasticInelastic: {
          fp1 = FormFactors::protonElastic( -t1_ );
          fp2 = FormFactors::protonInelastic( -t2_, w2_, my2, *cuts_.structure_functions );
        } break;
        case Kinematics::Mode::InelasticElastic: {
          fp1 = FormFactors::protonInelastic( -t1_, w1_, mx2, *cuts_.structure_functions );
          fp2 = FormFactors::protonElastic( -t2_ );
        } break;
        case Kinematics::Mode::InelasticInelastic: {
          fp1 = FormFactors::protonInelastic( -t1_, w1_, mx2, *cuts_.structure_functions );
          fp2 = FormFactors::protonInelastic( -t2_, w2_, my2, *cuts_.structure_functions );
        } break;
      }
    }
  }
}
