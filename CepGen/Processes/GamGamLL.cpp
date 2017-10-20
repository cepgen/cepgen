#include "GamGamLL.h"
#include "CepGen/Core/Exception.h"

using namespace CepGen::Process;

GamGamLL::GamGamLL( int nopt ) :
  GenericProcess( "lpair", "pp -> p(*) (gamma gamma -> l+ l-) p(*)" ),
  n_opt_( nopt ),
  MX2_( 0. ), MY2_( 0. ), Ml2_( 0. ),
  ep1_( 0. ), ep2_( 0. ), p_cm_( 0. ),
  w12_( 0. ), w31_( 0. ), dw31_( 0. ), w52_( 0. ), dw52_( 0. ),
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

void
GamGamLL::addEventContent()
{
  GenericProcess::setEventContent( {
    { Particle::IncomingBeam1, Proton },
    { Particle::IncomingBeam2, Proton },
    { Particle::Parton1, Photon },
    { Particle::Parton2, Photon }
  }, {
    { Particle::OutgoingBeam1, { Proton } },
    { Particle::OutgoingBeam2, { Proton } },
    { Particle::CentralSystem, { Muon, Muon } }
  } );
}

unsigned int
GamGamLL::numDimensions( const Kinematics::ProcessMode& process_mode ) const
{
  switch ( process_mode ) {
    case Kinematics::ElectronProton:     { InError( "Not supported yet!" ); }
    case Kinematics::ElasticElastic:
    default:                             return 7;
    case Kinematics::ElasticInelastic:
    case Kinematics::InelasticElastic:   return 8;
    case Kinematics::InelasticInelastic: return 9;
  }
}

bool
GamGamLL::pickin()
{
  DebuggingInsideLoop( Form( "Optimised mode? %i", n_opt_ ) );

  jacobian_ = 0.;

  w4_ = mc4_*mc4_;

  // sig1 = sigma and sig2 = sigma' in [1]
  const double sig = mc4_+MY_;
  double sig1 = sig*sig,
         sig2 = sig1;

  DebuggingInsideLoop( Form( "mc4 = %f\n\t"
                             "sig1 = %f\n\t"
                             "sig2 = %f", mc4_, sig1, sig2 ) );

  // Mass difference between the first outgoing particle and the first incoming
  // particle
  w31_ = MX2_-w1_;
  // Mass difference between the second outgoing particle and the second
  // incoming particle
  w52_ = MY2_-w2_;
  // Mass difference between the two incoming particles
  w12_ = w1_-w2_;
  // Mass difference between the central two-photons system and the second
  // outgoing particle
  const double d6 = w4_-MY2_;

  DebuggingInsideLoop( Form( "w1 = %f\n\t"
                             "w2 = %f\n\t"
                             "w3 = %f\n\t"
                             "w4 = %f\n\t"
                             "w5 = %f",
                             w1_, w2_, MX2_, w4_, MY2_ ) );

  DebuggingInsideLoop( Form( "w31 = %f\n\tw52 = %f\n\tw12 = %f", w31_, w52_, w12_ ) );

  const double ss = s_+w12_;

  const double rl1 = ss*ss-4.*w1_*s_; // lambda(s, m1**2, m2**2)
  if ( rl1 <= 0. ) { InWarning( Form( "rl1 = %f <= 0", rl1 ) ); return false; }
  sl1_ = sqrt( rl1 );

  s2_ = 0.;
  double ds2 = 0.;
  if ( n_opt_ == 0 ) {
    const double smax = s_+MX2_-2.*MX_*sqs_;
    Map( x(2), sig1, smax, s2_, ds2, "s2" );
    sig1 = s2_; //FIXME!!!!!!!!!!!!!!!!!!!!
  }

  DebuggingInsideLoop( Form( "s2 = %f", s2_ ) );

  //std::cout << "s=" << _s << ", w3=" << MX2_ << ", sig1=" << sig1 << std::endl;
  const double sp = s_+MX2_-sig1,
               d3 = sig1-w2_;

  const double rl2 = sp*sp-4.*s_*MX2_; // lambda(s, m3**2, sigma)
  if ( rl2 <= 0. ) { InWarning( Form( "rl2 = %f <= 0", rl2 ) ); return false; }
  const double sl2 = sqrt( rl2 );

  //std::cout << "ss=" << ss << ", sp=" << sp << ", sl1=" << sl1_ << ", sl2=" << sl2 << std::endl;
  double t1_max = w1_+MX2_-( ss*sp+sl1_*sl2 )/( 2.*s_ ), // definition from eq. (A.4) in [1]
         t1_min = ( w31_*d3+( d3-w31_ )*( d3*w1_-w31_*w2_ )/s_ )/t1_max; // definition from eq. (A.5) in [1]

  // FIXME dropped in CDF version
  const Kinematics::Limits q2_limits = cuts_.cuts.initial[Cuts::q2];
  if ( t1_max > -q2_limits.min() ) { InWarning( Form( "t1max = %f > -q2min = %f", t1_max, -q2_limits.min() ) ); return false; }
  if ( t1_min < -q2_limits.max() && q2_limits.hasMax() ) { Debugging( Form( "t1min = %f < -q2max = %f", t1_min, -q2_limits.max() ) ); return false; }
  if ( t1_max < -q2_limits.max() && q2_limits.hasMax() ) t1_max = -q2_limits.max();
  if ( t1_min > -q2_limits.min() )                       t1_min = -q2_limits.min();
  /////

  // t1, the first photon propagator, is defined here
  t1_ = 0.;
  double dt1 = 0.;
  Map( x(0), t1_min, t1_max, t1_, dt1, "t1" );
  // changes wrt mapt1 : dx->-dx
  dt1 = -dt1;

  DebuggingInsideLoop( Form( "Definition of t1 = %f according to\n\t"
                             "(t1min, t1max) = (%f, %f)", t1_, t1_min, t1_max ) );

  dd4_ = w4_-t1_;

  const double d8 = t1_-w2_,
	       t13 = t1_-w1_-MX2_;

  sa1_ = -pow( t1_-w31_, 2 )/4.+w1_*t1_;
  if ( sa1_ >= 0. ) { InWarning( Form( "sa1_ = %f >= 0", sa1_ ) ); return false; }

  const double sl3 = sqrt( -sa1_ );

  // one computes splus and (s2x=s2max)
  double splus, s2max;
  if ( w1_ != 0. ) {
    const double sb =( s_*( t1_-w31_ )+w12_*t13 )/( 2.*w1_ )+MX2_,
                 sd = sl1_*sl3/w1_,
                 se =( s_*( t1_*( s_+t13-w2_ )-w2_*w31_ )+MX2_*( w12_*d8+w2_*MX2_ ) )/w1_;

    if ( fabs( ( sb-sd )/sd )>=1. ) { splus = sb-sd; s2max = se/splus; }
    else                            { s2max = sb+sd; splus = se/s2max; }
  }
  else { // 3
    s2max = ( s_*(t1_*(s_+d8-MX2_)-w2_*MX2_ )+w2_*MX2_*( w2_+MX2_-t1_ ) )/( ss*t13 );
    splus = sig2;
  }
  // 4
  double s2x = s2max;

  DebuggingInsideLoop( Form( "s2x = s2max = %f", s2x ) );

  if ( n_opt_ < 0 ) { // 5
    if ( splus > sig2 ) {
      sig2 = splus;
      DebuggingInsideLoop( Form( "sig2 truncated to splus = %f", splus ) );
    }
    if ( n_opt_ < -1 ) { Map( x(2), sig2, s2max, s2_, ds2, "s2" ); }
    else               { Mapla( t1_, w2_, x(2), sig2, s2max, s2_, ds2 ); } // n_opt_==-1
    s2x = s2_;
  }
  else if ( n_opt_ == 0 ) { s2x = s2_; } // 6

  DebuggingInsideLoop( Form( "s2x = %f", s2x ) );

  // 7
  const double r1 = s2x-d8,
               r2 = s2x-d6;

  const double rl4 = ( r1*r1-4.*w2_*s2x )*( r2*r2-4.*MY2_*s2x );
  if ( rl4 <= 0. ) {
    DebuggingInsideLoop( Form( "rl4 = %f <= 0", rl4 ) );
    return false;
  }
  const double sl4 = sqrt( rl4 );

  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  const double t2_max = w2_+MY2_-( r1*r2+sl4 )/s2x * 0.5,
               t2_min = ( w52_*dd4_+( dd4_-w52_ )*( dd4_*w2_-w52_*t1_ )/s2x )/t2_max;

  // t2, the second photon propagator, is defined here
  t2_ = 0.;
  double dt2 = 0.;
  Map( x(1), t2_min, t2_max, t2_, dt2, "t2" );
  // changes wrt mapt2 : dx->-dx

  dt2 = -dt2;

  // \f$\delta_6=m_4^2-m_5^2\f$ as defined in Vermaseren's paper
  const double tau = t1_-t2_,
               r3 = dd4_-t2_,
               r4 = w52_-t2_;

  DebuggingInsideLoop( Form( "r1 = %f\n\tr2 = %f\n\tr3 = %f\n\tr4 = %f", r1, r2, r3, r4 ) );

  const double b = r3*r4-2.*( t1_+w2_ )*t2_,
               c = t2_*d6*d8+( d6-d8 )*( d6*w2_-d8*MY2_ );

  const double t25 = t2_-w2_-MY2_;

  sa2_ = -r4*r4/4.+w2_*t2_;
  if ( sa2_ >= 0. ) { InWarning( Form( "sa2_ = %f >= 0", sa2_ ) ); return false; }

  const double sl6 = 2.*sqrt( -sa2_ );

  g4_ = -r3*r3/4.+t1_*t2_;
  if ( g4_ >= 0. )  { InWarning( Form( "g4_ = %f >= 0", g4_ ) );   return false; }

  const double sl7 = 2.*sqrt( -g4_ ),
               sl5 = sl6*sl7;

  double s2p, s2min;
  if ( fabs( ( sl5-b )/sl5 ) >= 1. ) {
    s2p = ( sl5-b )/t2_ * 0.5;
    s2min = c/( t2_*s2p );
  }
  else { // 8
    s2min = ( -sl5-b )/t2_ * 0.5;
    s2p = c/( t2_*s2min );
  }
  // 9
  if ( n_opt_ > 1 )       Map( x( 2 ), s2min, s2max, s2_, ds2, "s2" );
  else if ( n_opt_ == 1 ) Mapla( t1_, w2_, x( 2 ), s2min, s2max, s2_, ds2 );

  const double ap = -0.25*pow( s2_+d8, 2 )+s2_*t1_;

  DebuggingInsideLoop( Form( "s2 = %f, s2max = %f, splus = %f", s2_, s2max, splus ) );

  if ( w1_!=0. ) dd1_ = -0.25 * ( s2_-s2max ) * ( s2_-splus ) * w1_; // 10
  else           dd1_ =  0.25 * ( s2_-s2max ) * ss * t13;
  // 11
  dd2_ = -t2_*( s2_-s2p )*( s2_-s2min ) * 0.25;

  DebuggingInsideLoop( Form( "t2 =%f\n\ts2 =%f\n\ts2p=%f\n\ts2min=%f\n\tdd2=%f", t2_, s2_, s2p, s2min, dd2_ ) );

  const double yy4 = cos( M_PI*x( 3 ) );
  const double dd = dd1_*dd2_;
  p12_ = ( s_-w1_-w2_ )*0.5;
  const double st = s2_-t1_-w2_;
  const double delb = ( 2.*w2_*r3+r4*st )*( 4.*p12_*t1_-( t1_-w31_ )*st )/( 16.*ap );

  if ( dd <= 0. ) {
    DebuggingInsideLoop( Form( "dd = %e <= 0\n\tdd1 = %e\tdd2 = %e", dd, dd1_, dd2_ ) );
    return false;
  }

  delta_ = delb - yy4*st*sqrt( dd )/ ap * 0.5;
  s1_ = t2_+w1_+( 2.*p12_*r3-4.*delta_ )/st;

  if ( ap >= 0. ) {
    DebuggingInsideLoop( Form( "ap = %f >= 0", ap ) );
    return false;
  }

  jacobian_ = ds2
            * dt1
            * dt2
            * M_PI*M_PI/( 8.*sl1_*sqrt( -ap ) );

  DebuggingInsideLoop( Form( "Jacobian = %e", jacobian_ ) );

  gram_ = ( 1.-yy4*yy4 )*dd/ap;

  p13_ = -t13 * 0.5;
  p14_ = (tau+s1_-MX2_) * 0.5;
  p25_ = -t25 * 0.5;

  /*const double //p15 = (s_+t2_-s1_-w2_)/2.,
               //p23 = (s_+t1_-s2_-w1_)/2.,
               //p24 = (s2_-tau-MY2_)/2.,
               //p34 = (s1_-MX2_-w4_)/2.,
               //p35 = (s_+w4_-s1_-s2_)/2.,
               //p45 = (s2_-w4_-MY2_)/2.;*/

  p1k2_ = ( s1_-t2_-w1_ ) * 0.5;
  p2k1_ = st * 0.5;

  double s1p, s1m;
  if ( w2_ != 0. ) {
    const double sbb = ( s_*( t2_-w52_ )-w12_*t25 )/w2_ * 0.5 + MY2_,
                 sdd = sl1_*sl6/w2_ * 0.5,
                 see = ( s_*( t2_*( s_+t25-w1_ )-w1_*w52_ )+MY2_*( w1_*MY2_-w12_*( t2_-w1_ ) ) )/w2_;
    if ( sbb/sdd >= 0. ) { s1p = sbb+sdd; s1m = see/s1p; }
    else                 { s1m = sbb-sdd; s1p = see/s1m; } // 12
    dd3_ = -w2_*( s1p-s1_ )*( s1m-s1_ ) * 0.25; // 13
  }
  else { // 14
    s1p = ( s_*( t2_*( s_-MY2_+t2_-w1_ )-w1_*MY2_ )+w1_*MY2_*( w1_+MY2_-t2_ ) )/( t25*( s_-w12_ ) );
    dd3_ = -t25*( s_-w12_ )*( s1p-s1_ ) * 0.25;
  }
  // 15
  //const double acc3 = (s1p-s1_)/(s1p+s1_);

  const double ssb = t2_+w1_-r3*( w31_-t1_ )/t1_ * 0.5,
               ssd = sl3*sl7/t1_,
               sse = ( t2_-w1_ )*( w4_-MX2_ )+( t2_-w4_+w31_ )*( ( t2_-w1_)*MX2_-(w4_-MX2_)*w1_)/t1_;

  double s1pp, s1pm;
  if (ssb/ssd>=0.) { s1pp = ssb+ssd; s1pm = sse/s1pp; }
  else             { s1pm = ssb-ssd; s1pp = sse/s1pm; } // 16
  // 17
  dd4_ = -t1_*( s1_-s1pp )*( s1_-s1pm ) * 0.25;
  //const double acc4 = ( s1_-s1pm )/( s1_+s1pm );
  dd5_ = dd1_+dd3_+( ( p12_*( t1_-w31_ )*0.5-w1_*p2k1_ )*( p2k1_*( t2_-w52_ )-w2_*r3 )-delta_*( 2.*p12_*p2k1_-w2_*( t1_-w31_ ) ) ) / p2k1_;

  return true;
}

bool
GamGamLL::orient()
{
  if ( !pickin() or jacobian_ == 0. ) { DebuggingInsideLoop( Form( "Pickin failed! Jacobian = %f", jacobian_ ) ); return false; }

  const double re = 0.5 / sqs_;
  ep1_ = re*( s_+w12_ );
  ep2_ = re*( s_-w12_ );

  DebuggingInsideLoop( Form( " re = %e\n\tw12_ = %e", re, w12_ ) );
  DebuggingInsideLoop( Form( "Incoming particles' energy = %f, %f", ep1_, ep2_ ) );

  p_cm_ = re*sl1_;

  de3_ = re*(s2_-MX2_+w12_);
  de5_ = re*(s1_-MY2_-w12_);

  // Final state energies
  const double ep3 = ep1_-de3_,
               ep5 = ep2_-de5_;
  ec4_ = de3_+de5_;

  if ( ec4_ < mc4_ ) { InWarning( Form( "ec4_ = %f < mc4_ = %f\n\t==> de3 = %f, de5 = %f", ec4_, mc4_, de3_, de5_ ) ); return false; }

  // What if the protons' momenta are not along the z-axis?
  pc4_ = sqrt( ec4_*ec4_-mc4_*mc4_ );

  if ( pc4_==0. ) {
    InWarning( "pzc4==0" );
    return false;
  }

  const double pp3 = sqrt( ep3*ep3-MX2_ ), pt3 = sqrt( dd1_/s_ )/p_cm_,
               pp5 = sqrt( ep5*ep5-MY2_ ), pt5 = sqrt( dd3_/s_ )/p_cm_;

  DebuggingInsideLoop( Form( "Central system's energy: E4 = %f\n\t"
                             "                 momentum: p4 = %f\n\t"
                             "                 invariant mass: m4 = %f\n\t"
                             "Outgoing particles' energy: E3 = %f\n\t"
                             "                            E5 = %f",
                             ec4_, pc4_, mc4_, ep3, ep5 ) );

  const double sin_theta3 = pt3/pp3,
               sin_theta5 = pt5/pp5;

  DebuggingInsideLoop( Form( "sin_theta3 = %e\n\tsin_theta5 = %e", sin_theta3, sin_theta5 ) );

  if (sin_theta3>1.) { InWarning( Form( "sin_theta3 = %e > 1", sin_theta3 ) ); return false; }
  if (sin_theta5>1.) { InWarning( Form( "sin_theta5 = %e > 1", sin_theta5 ) ); return false; }

  double ct3 = sqrt( 1.-sin_theta3*sin_theta3 ),
         ct5 = sqrt( 1.-sin_theta5*sin_theta5 );

  if ( ep1_*ep3 < p13_ ) ct3 *= -1.;
  if ( ep2_*ep5 > p25_ ) ct5 *= -1.;

  DebuggingInsideLoop( Form( "ct3 = %e\n\tct5 = %e", ct3, ct5 ) );

  if ( dd5_ < 0. ) { InWarning( Form( "dd5 = %f < 0", dd5_ ) ); return false; }

  // Centre of mass system kinematics (theta4 and phi4)
  pt4_ = sqrt( dd5_/s_ )/p_cm_;
  sin_theta4_ = pt4_/pc4_;

  if ( sin_theta4_ > 1. ) { InWarning( Form( "st4 = %f > 1", sin_theta4_ ) ); return false; }

  cos_theta4_ = sqrt( 1.-sin_theta4_*sin_theta4_ );
  if ( ep1_*ec4_ < p14_ ) cos_theta4_ *= -1.;

  al4_ = 1.-cos_theta4_;
  be4_ = 1.+cos_theta4_;

  if (cos_theta4_<0.) be4_ = sin_theta4_*sin_theta4_/al4_;
  else                al4_ = sin_theta4_*sin_theta4_/be4_;

  DebuggingInsideLoop( Form( "ct4 = %f\n\tal4 = %f, be4 = %f", cos_theta4_, al4_, be4_ ) );

  const double rr  = sqrt( -gram_/s_ )/( p_cm_*pt4_ );
  const double sin_phi3 =  rr / pt3,
               sin_phi5 = -rr / pt5;

  if ( fabs( sin_phi3 ) > 1. ) { InWarning( Form( "sin(phi_3) = %e while it must be in [-1 ; 1]", sin_phi3 ) ); return false; }
  if ( fabs( sin_phi5 ) > 1. ) { InWarning( Form( "sin(phi_5) = %e while it must be in [-1 ; 1]", sin_phi5 ) ); return false; }

  const double cos_phi3 = -sqrt( 1.-sin_phi3*sin_phi3 ),
               cos_phi5 = -sqrt( 1.-sin_phi5*sin_phi5 );

  p3_lab_ = Particle::Momentum( pp3*sin_theta3*cos_phi3, pp3*sin_theta3*sin_phi3, pp3*ct3, ep3 );
  p5_lab_ = Particle::Momentum( pp5*sin_theta5*cos_phi5, pp5*sin_theta5*sin_phi5, pp5*ct5, ep5 );

  const double a1 = p3_lab_.px()-p5_lab_.px();

  DebuggingInsideLoop( Form( "Kinematic quantities\n\t"
                             "cos(theta3) = %1.4f\tsin(theta3) = %1.4f\tcos( phi3 ) = %1.4f\tsin( phi3 ) = %1.4f\n\t"
                             "cos(theta4) = %1.4f\tsin(theta4) = %1.4f\n\t"
                             "cos(theta5) = %1.4f\tsin(theta5) = %1.4f\tcos( phi5 ) = %1.4f\tsin( phi5 ) = %1.4f\n\t"
                             "a1 = %f",
                             ct3, sin_theta3, cos_phi3, sin_phi3,
                             cos_theta4_, sin_theta4_,
                             ct5, sin_theta5, cos_phi5, sin_phi5, a1 ) );

  if ( fabs( pt4_+p3_lab_.px()+p5_lab_.px() ) < fabs( fabs( a1 )-pt4_ ) ) {
    DebuggingInsideLoop( Form( "|pp4+pp3*cos(phi3)+pp5*cos(phi5)| < | |a1|-pp4 |\n\t"
                               "pp4 = %f\tpp5 = %f\n\t"
                               "cos(phi3) = %f\tcos(phi5) = %f"
                               "a1 = %f",
                               pt4_, pt5, cos_phi3, cos_phi5, a1 ) );
    return true;
  }
  if ( a1 < 0. ) p5_lab_.setP( 0, -p5_lab_.px() );
  else           p3_lab_.setP( 0, -p3_lab_.px() );
  return true;
}

double
GamGamLL::computeOutgoingPrimaryParticlesMasses( double x, double outmass, double lepmass, double& dw )
{
  const double mx0 = ParticleProperties::mass( Proton )+ParticleProperties::mass( PiPlus ); // 1.07
  const Kinematics::Limits mx_limits = cuts_.cuts.remnants[Cuts::mass];
  const double wx2min = pow( std::max( mx0, mx_limits.min() ), 2 ),
               wx2max = pow( std::min( sqs_-outmass-2.*lepmass, mx_limits.max() ), 2 );

  double mx2 = 0., dmx2 = 0.;
  Map( x, wx2min, wx2max, mx2, dmx2, "mx2" );

  DebuggingInsideLoop( Form( "mX^2 in range (%f, %f), x = %f\n\t"
                             "mX^2 = %f, d(mX^2) = %f\n\t"
                             "mX = %f, d(mX) = %f", wx2min, wx2max, x, mx2, dmx2, sqrt( mx2 ), sqrt( dmx2 ) ) );

  dw = sqrt( dmx2 );
  return sqrt( mx2 );
}

void
GamGamLL::beforeComputeWeight()
{
  if ( !GenericProcess::is_point_set_ ) return;

  const Particle& p1 = event_->getOneByRole( Particle::IncomingBeam1 ),
                 &p2 = event_->getOneByRole( Particle::IncomingBeam2 );

  ep1_ = p1.energy();
  ep2_ = p2.energy();

  Ml2_ = event_->getByRole( Particle::CentralSystem )[0].mass2();

  switch ( cuts_.mode ) {
    case Kinematics::ElectronProton: default:
      { InError( "Case not yet supported!" ); } break;
    case Kinematics::ElasticElastic:
      dw31_ = dw52_ = 0.; break;
    case Kinematics::InelasticElastic: {
      const double m = computeOutgoingPrimaryParticlesMasses( x( 7 ), p1.mass(), sqrt( Ml2_ ), dw31_ );
      event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( m );
      event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( ParticleProperties::mass( p2.pdgId() ) );
    } break;
    case Kinematics::ElasticInelastic: {
      const double m = computeOutgoingPrimaryParticlesMasses( x( 7 ), p2.mass(), sqrt( Ml2_ ), dw52_ );
      event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( ParticleProperties::mass( p1.pdgId() ) );
      event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( m );
    } break;
    case Kinematics::InelasticInelastic: {
      const double mx = computeOutgoingPrimaryParticlesMasses( x( 7 ), p2.mass(), sqrt( Ml2_ ), dw31_ );
      event_->getOneByRole( Particle::OutgoingBeam1 ).setMass( mx );
      const double my = computeOutgoingPrimaryParticlesMasses( x( 8 ), p1.mass(), sqrt( Ml2_ ), dw52_ );
      event_->getOneByRole( Particle::OutgoingBeam2 ).setMass( my );
    } break;
  }
  MX_ = event_->getOneByRole( Particle::OutgoingBeam1 ).mass();
  MY_ = event_->getOneByRole( Particle::OutgoingBeam2 ).mass();
  MX2_ = MX_*MX_;
  MY2_ = MY_*MY_;
}

double
GamGamLL::computeWeight()
{
  if ( !is_outgoing_state_set_ ) { InWarning( "Output state not set!" ); return 0.; }

  DebuggingInsideLoop( Form( "sqrt(s)=%f\n\tm(X1)=%f\tm(X2)=%f", sqs_, MX_, MY_ ) );

  Kinematics::Limits& w_limits = cuts_.cuts.initial[Cuts::w];
  if ( !w_limits.hasMax() ) w_limits.max() = s_;

  // The minimal energy for the central system is its outgoing leptons' mass energy (or wmin_ if specified)
  const double wmin = std::max( 4.*Ml2_, w_limits.min() );

  // The maximal energy for the central system is its CM energy with the outgoing particles' mass energy substracted (or _wmax if specified)
  const double wmax = std::min( pow( sqs_-MX_-MY_, 2 ), w_limits.max() );

  DebuggingInsideLoop( Form( "wmin = %f\n\twmax = %f\n\twmax/wmin = %f", wmin, wmax, wmax/wmin ) );

  // compute the two-photon energy for this point
  w4_ = 0.;
  double dw4 = 0.;
  Map( x( 4 ), wmin, wmax, w4_, dw4, "w4" );
  Map( x( 4 ), wmin, wmax, w4_, dw4, "w4" );
  mc4_ = sqrt( w4_ );

  DebuggingInsideLoop( Form( "Computed value for w4 = %f -> mc4 = %f", w4_, mc4_ ) );

  if ( !orient() ) return 0.;

  if ( jacobian_ == 0. ) { InWarning( Form( "dj = %f", jacobian_ ) ); return 0.; }

  if ( t1_>0. )  { InWarning( Form( "t1 = %f > 0", t1_ ) ); return 0.; }
  if ( t2_>0. )  { InWarning( Form( "t2 = %f > 0", t2_ ) ); return 0.; }

  const double ecm6 = w4_ / ( 2.*mc4_ ),
               pp6cm = sqrt( ecm6*ecm6-Ml2_ );

  jacobian_ *= dw4*pp6cm/( mc4_*Constants::sconstb*s_ );

  // Let the most obscure part of this code begin...

  const double e1mp1 = w1_ / ( ep1_+p_cm_ ),
               e3mp3 = MX2_ / ( p3_lab_.energy()+p3_lab_.p() );

  const double al3 = pow( sin( p3_lab_.theta() ), 2 )/( 1.+( p3_lab_.theta() ) );

  // 2-photon system kinematics ?!
  const double eg = ( w4_+t1_-t2_ )/( 2.*mc4_ );
  double pg = sqrt( eg*eg-t1_ );

  const double pgx = -p3_lab_.px()*cos_theta4_-sin_theta4_*( de3_-e1mp1 + e3mp3 + p3_lab_.p()*al3 ),
               pgy = -p3_lab_.py(),
               pgz = mc4_*de3_/( ec4_+pc4_ )-ec4_*de3_*al4_/mc4_-p3_lab_.px()*ec4_*sin_theta4_/mc4_+ec4_*cos_theta4_/mc4_*( p3_lab_.p()*al3+e3mp3-e1mp1 );

  DebuggingInsideLoop( Form( "pg3 = (%f, %f, %f)\n\t"
                             "pg3^2 = %f",
                             pgx, pgy, pgz,
                             sqrt( pgx*pgx+pgy*pgy+pgz*pgz ) ) );

  const double pgp = sqrt( pgx*pgx + pgy*pgy ), // outgoing proton (3)'s transverse momentum
               pgg = sqrt( pgp*pgp + pgz*pgz ); // outgoing proton (3)'s momentum
  if ( pgg > pgp*0.9 && pgg > pg ) { pg = pgg; } //FIXME ???

  // Phi angle for the 2-photon system ?!
  const double cpg = pgx/pgp,
               spg = pgy/pgp;

  // Theta angle for the 2-photon system ?!
  const double stg = pgp/pg;

  const int theta_sign = ( pgz>0. ) ? 1 : -1;
  const double ctg = theta_sign*sqrt( 1.-stg*stg );

  double xx6 = x( 5 );

  const double amap = 0.5 * (w4_-t1_-t2_),
               bmap = 0.5 * sqrt( ( pow( w4_-t1_-t2_, 2 )-4.*t1_*t2_ )*( 1.-4.*Ml2_/w4_ ) ),
               ymap = ( amap+bmap )/( amap-bmap ),
               beta = pow( ymap, 2.*xx6-1. );
  xx6 = 0.5 * ( 1. + amap/bmap*( beta-1. )/( beta+1. ) );
  xx6 = std::max( 0., std::min( xx6, 1. ) ); // xx6 in [0., 1.]

  DebuggingInsideLoop( Form( "amap = %f\n\tbmap = %f\n\tymap = %f\n\tbeta = %f", amap, bmap, ymap, beta ) );

  // 3D rotation of the first outgoing lepton wrt the CM system
  const double theta6cm = acos( 1.-2.*xx6 );

  // match the Jacobian
  jacobian_ *= ( amap+bmap*cos( theta6cm ) );
  jacobian_ *= ( amap-bmap*cos( theta6cm ) );
  jacobian_ /= amap;
  jacobian_ /= bmap;
  jacobian_ *= log( ymap );
  jacobian_ *= 0.5;

  DebuggingInsideLoop( Form( "Jacobian = %e", jacobian_ ) );

  DebuggingInsideLoop( Form( "ctcm6 = %f\n\tstcm6 = %f", cos( theta6cm ), sin( theta6cm ) ) );

  const double phi6cm = 2.*M_PI*x( 6 );

  // First outgoing lepton's 3-momentum in the centre of mass system
  Particle::Momentum p6cm = Particle::Momentum::fromPThetaPhi( pp6cm, theta6cm, phi6cm );

  DebuggingInsideLoop( Form( "p3cm6 = (%f, %f, %f)", p6cm.px(), p6cm.py(), p6cm.pz() ) );

  const double h1 = stg*p6cm.pz()+ctg*p6cm.px();
  const double pc6z = ctg*p6cm.pz()-stg*p6cm.px(), pc6x = cpg*h1-spg*p6cm.py();

  const double qcx = 2.*pc6x, qcz = 2.*pc6z;
  // qcy == QCY is never defined

  const double el6 = ( ec4_*ecm6+pc4_*pc6z ) / mc4_,
               h2  = ( ec4_*pc6z+pc4_*ecm6 ) / mc4_;

  DebuggingInsideLoop( Form( "h1 = %f\n\th2 = %f", h1, h2 ) );

  // first outgoing lepton's 3-momentum
  const double p6x = cos_theta4_*pc6x+sin_theta4_*h2,
               p6y = cpg*p6cm.py()+spg*h1,
               p6z = cos_theta4_*h2-sin_theta4_*pc6x;

  // first outgoing lepton's kinematics
  p6_cm_ = Particle::Momentum( p6x, p6y, p6z, el6 );
  DebuggingInsideLoop( Form( "E6(cm) = %f\n\tP6(cm) = (%f, %f, %f)", el6, p6x, p6y, p6z ) );

  const double hq = ec4_*qcz/mc4_;

  const Particle::Momentum qve(
    cos_theta4_*qcx+sin_theta4_*hq,
    2.*p6y,
    cos_theta4_*hq-sin_theta4_*qcx,
    pc4_*qcz/mc4_ // energy
  );

  // Available energy for the second lepton is the 2-photon system's energy with the first lepton's energy removed
  const double el7 = ec4_-el6;

  DebuggingInsideLoop( Form( "Outgoing kinematics\n\t"
                             " first outgoing lepton: p = %f, E = %f\n\t"
                             "second outgoing lepton: p = %f, E = %f",
                             p6_cm_.p(), p6_cm_.energy(), p7_cm_.p(), p6_cm_.energy() ) );

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

  bb_ = t1_*t2_+( w4_*pow( sin( theta6cm ), 2 ) + 4.*Ml2_*pow( cos( theta6cm ), 2 ) )*pg*pg;

  const double c1 = p3_lab_.pt() * ( qve.px()*sin_phi3  - qve.py()*cos_phi3   ),
               c2 = p3_lab_.pt() * ( qve.pz()*ep1_ - qve.energy() *p_cm_ ),
               c3 = ( w31_*ep1_*ep1_ + 2.*w1_*de3_*ep1_ - w1_*de3_*de3_ + p3_lab_.pt2()*ep1_*ep1_ ) / ( p3_lab_.energy()*p_cm_ + p3_lab_.pz()*ep1_ );

  const double b1 = p5_lab_.pt() * ( qve.px()*sin_phi5  - qve.py()*cos_phi5   ),
               b2 = p5_lab_.pt() * ( qve.pz()*ep2_ + qve.energy() *p_cm_ ),
               b3 = ( w52_*ep2_*ep2_ + 2.*w2_*de5_*ep2_ - w2_*de5_*de5_ + p5_lab_.pt2()*ep2_*ep2_ ) / ( ep2_*p5_lab_.pz() - p5_lab_.energy()*p_cm_ );

  const double r12 =  c2*sin_phi3 + qve.py()*c3,
               r13 = -c2*cos_phi3 - qve.px()*c3;

  const double r22 =  b2*sin_phi5 + qve.py()*b3,
               r23 = -b2*cos_phi5 - qve.px()*b3;

  epsi_ = p12_*c1*b1 + r12*r22 + r13*r23;

  g5_ = w1_*c1*c1 + r12*r12 + r13*r13;
  g6_ = w2_*b1*b1 + r22*r22 + r23*r23;

  a5_ = -( qve.px()*cos_phi3 + qve.py()*sin_phi3 )*p3_lab_.pt()*p1k2_-( ep1_*qve.energy()-p_cm_*qve.pz() )*( cos_phi3*cos_phi5 + sin_phi3*sin_phi5 )*p3_lab_.pt()*p5_lab_.pt() + ( de5_*qve.pz()+qve.energy()*( p_cm_+p5_lab_.pz() ) )*c3;
  a6_ = -( qve.px()*cos_phi5 + qve.py()*sin_phi5 )*p5_lab_.pt()*p2k1_-( ep2_*qve.energy()+p_cm_*qve.pz() )*( cos_phi3*cos_phi5 + sin_phi3*sin_phi5 )*p3_lab_.pt()*p5_lab_.pt() + ( de3_*qve.pz()-qve.energy()*( p_cm_-p3_lab_.pz() ) )*b3;

  DebuggingInsideLoop( Form( "a5 = %f\n\ta6 = %f", a5_, a6_ ) );

  ////////////////////////////////////////////////////////////////
  // END of GAMGAMLL subroutine in the FORTRAN version
  ////////////////////////////////////////////////////////////////

  const Particle::Momentum cm = event_->getOneByRole( Particle::IncomingBeam1 ).momentum() + event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

  ////////////////////////////////////////////////////////////////
  // INFO from f.f
  ////////////////////////////////////////////////////////////////

  const double gamma = cm.energy() / sqs_, betgam = cm.pz() / sqs_;

  //--- kinematics computation for both leptons

  p6_cm_.betaGammaBoost( gamma, betgam );
  p7_cm_.betaGammaBoost( gamma, betgam );

  //--- cut on mass of final hadronic system (MX/Y)

  const Kinematics::Limits mx_limits = cuts_.cuts.remnants[Cuts::mass];
  if ( cuts_.mode == Kinematics::InelasticElastic || cuts_.mode == Kinematics::InelasticInelastic ) {
    if ( mx_limits.hasMin() && MX_ < mx_limits.min() ) return 0.;
    if ( mx_limits.hasMax() && MX_ > mx_limits.max() ) return 0.;
  }
  if ( cuts_.mode == Kinematics::ElasticInelastic || cuts_.mode == Kinematics::InelasticInelastic ) {
    if ( mx_limits.hasMin() && MY_ < mx_limits.min() ) return 0.;
    if ( mx_limits.hasMax() && MY_ > mx_limits.max() ) return 0.;
  }

  //--- cut on the proton's Q2 (first photon propagator T1)

  const Kinematics::Limits q2_limits = cuts_.cuts.initial[Cuts::q2];
  if ( q2_limits.hasMax() && t1_<-q2_limits.max() ) return 0;
  if ( q2_limits.hasMin() && t1_>-q2_limits.min() ) return 0.;

  //--- cuts on outgoing leptons' kinematics

  const Kinematics::Limits m_limits = cuts_.cuts.central[Cuts::mass_sum];
  if ( m_limits.hasMin() && ( p6_cm_+p7_cm_ ).mass() < m_limits.min() ) return 0.;
  if ( m_limits.hasMax() && ( p6_cm_+p7_cm_ ).mass() > m_limits.max() ) return 0.;

  //----- cuts on the individual leptons

  const Kinematics::Limits pt_limits = cuts_.cuts.central[Cuts::pt_single];
  if ( pt_limits.hasMin() && ( p6_cm_.pt() < pt_limits.min() || p7_cm_.pt() < pt_limits.min() ) ) return 0.;
  if ( pt_limits.hasMax() && ( p6_cm_.pt() > pt_limits.max() || p7_cm_.pt() > pt_limits.max() ) ) return 0.;

  const Kinematics::Limits energy_limits = cuts_.cuts.central[Cuts::energy_single];
  if ( energy_limits.hasMin() && ( p6_cm_.energy() < energy_limits.min() || p7_cm_.energy() < energy_limits.min() ) ) return 0.;
  if ( energy_limits.hasMax() && ( p6_cm_.energy() > energy_limits.max() || p7_cm_.energy() > energy_limits.max() ) ) return 0.;

  const Kinematics::Limits eta_limits = cuts_.cuts.central[Cuts::eta_single];
  if ( eta_limits.hasMin() && ( p6_cm_.eta() < eta_limits.min() || p7_cm_.eta() < eta_limits.min() ) ) return 0.;
  if ( eta_limits.hasMax() && ( p6_cm_.eta() > eta_limits.max() || p7_cm_.eta() > eta_limits.max() ) ) return 0.;

  //--- compute the structure functions factors

  switch ( cuts_.mode ) { // inherited from CDF version
    case Kinematics::ElectronProton: default: jacobian_ *= periPP( 1, 2 ); break; // ep case
    case Kinematics::ElasticElastic:          jacobian_ *= periPP( 2, 2 ); break; // elastic case
    case Kinematics::InelasticElastic:        jacobian_ *= periPP( 3, 2 )*( dw31_*dw31_ ); break;
    case Kinematics::ElasticInelastic:        jacobian_ *= periPP( 3, 2 )*( dw52_*dw52_ ); break; // single-dissociative case
    case Kinematics::InelasticInelastic:      jacobian_ *= periPP( 3, 3 )*( dw31_*dw31_ )*( dw52_*dw52_ ); break; // double-dissociative case
  }

  //--- compute the event weight using the Jacobian

  return Constants::GeV2toBarn*jacobian_;
}

void
GamGamLL::fillKinematics( bool )
{
  const Particle::Momentum cm = event_->getOneByRole( Particle::IncomingBeam1 ).momentum() + event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

  const double gamma  = cm.energy()/sqs_, betgam = cm.pz()/sqs_;

  Particle::Momentum plab_ip1 = Particle::Momentum( 0., 0.,  p_cm_, ep1_ ).betaGammaBoost( gamma, betgam );
  Particle::Momentum plab_ip2 = Particle::Momentum( 0., 0., -p_cm_, ep2_ ).betaGammaBoost( gamma, betgam );
  p3_lab_.betaGammaBoost( gamma, betgam );
  p5_lab_.betaGammaBoost( gamma, betgam );

  //----- needed to parametrise a random rotation around z-axis
  const int rany = ( rand() >= RAND_MAX*0.5 ) ? 1 : -1,
            ransign = ( rand() >= RAND_MAX*0.5 ) ? 1 : -1;
  const double ranphi = ( rand()/RAND_MAX*0.5 )*2.*M_PI;

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
    case Kinematics::ElasticElastic:
    case Kinematics::ElasticInelastic:
    default:
      op1.setStatus( Particle::FinalState ); // stable proton
      break;
    case Kinematics::InelasticElastic:
    case Kinematics::InelasticInelastic:
      op1.setStatus( Particle::Undecayed ); // fragmenting remnants
      op1.setMass( MX_ );
      break;
  }

  //----- second outgoing proton
  Particle& op2 = event_->getOneByRole( Particle::OutgoingBeam2 );
  op2.setMomentum( p5_lab_ );
  switch ( cuts_.mode ) {
    case Kinematics::ElasticElastic:
    case Kinematics::InelasticElastic:
    default:
      op2.setStatus( Particle::FinalState ); // stable proton
      break;
    case Kinematics::ElasticInelastic:
    case Kinematics::InelasticInelastic:
      op2.setStatus( Particle::Undecayed ); // fragmenting remnants
      op2.setMass( MY_ );
      break;
  }

  //----- first incoming photon
  Particle& ph1 = event_->getOneByRole( Particle::Parton1 );
  ph1.setMomentum( plab_ph1 );
  ph1.setStatus( Particle::Incoming );

  //----- second incoming photon
  Particle& ph2 = event_->getOneByRole( Particle::Parton2 );
  ph2.setMomentum( plab_ph2 );
  ph2.setStatus( Particle::Incoming );

  Particles& central_system = event_->getByRole( Particle::CentralSystem );

  //----- first outgoing lepton
  Particle& ol1 = central_system[0];
  ol1.setPdgId( ol1.pdgId(), ransign );
  ol1.setMomentum( p6_cm_ );
  ol1.setStatus( Particle::FinalState );

  //----- second outgoing lepton
  Particle& ol2 = central_system[1];
  ol2.setPdgId( ol2.pdgId(), -ransign );
  ol2.setMomentum( p7_cm_ );
  ol2.setStatus( Particle::FinalState );

  //----- intermediate two-lepton system
  event_->getOneByRole( Particle::Intermediate ).setMomentum( p6_cm_+p7_cm_ );
}

double
GamGamLL::periPP( int nup_, int ndown_ )
{
  DebuggingInsideLoop( Form( " Nup  = %d\n\tNdown = %d", nup_, ndown_ ) );

  FormFactors fp1, fp2;
  GenericProcess::formFactors( -t1_, -t2_, fp1, fp2 );

  DebuggingInsideLoop( Form( "u1 = %f\n\tu2 = %f\n\tv1 = %f\n\tv2 = %f", fp1.FM, fp1.FE, fp2.FM, fp2.FE ) );

  const double qqq = q1dq_*q1dq_,
               qdq = 4.*Ml2_-w4_;
  const double t11 = 64. *(  bb_*( qqq-g4_-qdq*( t1_+t2_+2.*Ml2_ ) )-2.*( t1_+2.*Ml2_ )*( t2_+2.*Ml2_ )*qqq ) * t1_*t2_, // magnetic-magnetic
               t12 = 128.*( -bb_*( dd2_+g6_ )-2.*( t1_+2.*Ml2_ )*( sa2_*qqq+a6_*a6_ ) ) * t1_, // electric-magnetic
               t21 = 128.*( -bb_*( dd4_+g5_ )-2.*( t2_+2.*Ml2_ )*( sa1_*qqq+a5_*a5_ ) ) * t2_, // magnetic-electric
               t22 = 512.*(  bb_*( delta_*delta_-gram_ )-pow( epsi_-delta_*( qdq+q1dq2_ ), 2 )-sa1_*a6_*a6_-sa2_*a5_*a5_-sa1_*sa2_*qqq ); // electric-electric

  const double peripp = ( fp1.FM*fp2.FM*t11 + fp1.FE*fp2.FM*t21 + fp1.FM*fp2.FE*t12 + fp1.FE*fp2.FE*t22 ) / pow( 2.*t1_*t2_*bb_, 2 );

  DebuggingInsideLoop( Form( "t11 = %5.2f\tt12 = %5.2f\n\t"
                             "t21 = %5.2f\tt22 = %5.2f\n\t"
                             "=> PeriPP = %e", t11, t12, t21, t22, peripp ) );

  return peripp;
}
