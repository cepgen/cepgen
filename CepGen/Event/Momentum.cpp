#include "CepGen/Event/Particle.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <math.h>
#include <iomanip>

namespace cepgen
{
  using Momentum = Particle::Momentum;

  //----- Particle momentum methods

  Momentum::Momentum() :
    px_( 0. ), py_( 0. ), pz_( 0. ), p_( 0. ), energy_( 0. )
  {}

  Momentum::Momentum( double x, double y, double z, double t ) :
    px_( x ), py_( y ), pz_( z ), energy_( t )
  {
    computeP();
  }

  Momentum::Momentum( double* p ) :
    px_( p[0] ), py_( p[1] ), pz_( p[2] ), energy_( p[3] )
  {
    computeP();
  }

  //--- static constructors

  Momentum
  Momentum::fromPtEtaPhi( double pt, double eta, double phi, double e )
  {
    const double px = pt*cos( phi ),
                 py = pt*sin( phi ),
                 pz = pt*sinh( eta );
    return Momentum( px, py, pz, e );
  }

  Momentum
  Momentum::fromPThetaPhi( double p, double theta, double phi, double e )
  {
    const double px = p*sin( theta )*cos( phi ),
                 py = p*sin( theta )*sin( phi ),
                 pz = p*cos( theta );
    return Momentum( px, py, pz, e );
  }

  Momentum
  Momentum::fromPxPyPzE( double px, double py, double pz, double e )
  {
    return Momentum( px, py, pz, e );
  }

  Momentum
  Momentum::fromPxPyPzM( double px, double py, double pz, double m )
  {
    Momentum mom( px, py, pz );
    mom.setMass( m );
    return mom;
  }

  Momentum
  Momentum::fromPxPyYM( double px, double py, double rap, double m )
  {
    const double pt = std::hypot( px, py ), et = std::hypot( pt, m );
    return Momentum( px, py, et*sinh( rap ), et*cosh( rap ) );
  }

  //--- arithmetic operators

  Momentum&
  Momentum::operator+=( const Momentum& mom )
  {
    px_ += mom.px_;
    py_ += mom.py_;
    pz_ += mom.pz_;
    energy_ += mom.energy_;
    computeP();
    return *this;
  }

  Momentum&
  Momentum::operator-=( const Momentum& mom )
  {
    px_ -= mom.px_;
    py_ -= mom.py_;
    pz_ -= mom.pz_;
    energy_ -= mom.energy_;
    computeP();
    return *this;
  }

  bool
  Momentum::operator==( const Momentum& mom ) const
  {
    return ( px_ == mom.px_
          && py_ == mom.py_
          && pz_ == mom.pz_
          && energy_ == mom.energy_ );
  }

  double
  Momentum::threeProduct( const Momentum& mom ) const
  {
    CG_DEBUG_LOOP( "Momentum" )
      << "  (" << px_ << ", " << py_ << ", " << pz_ << ")\n\t"
      << "* (" << mom.px_ << ", " << mom.py_ << ", " << mom.pz_ << ")\n\t"
      << "= " << px_*mom.px_+py_*mom.py_+pz_*mom.pz_;
    return px_*mom.px_+py_*mom.py_+pz_*mom.pz_;
  }

  double
  Momentum::fourProduct( const Momentum& mom ) const
  {
    CG_DEBUG_LOOP( "Momentum" )
      << "  (" << px_ << ", " << py_ << ", " << pz_ << ", " << energy_ << ")\n\t"
      << "* (" << mom.px_ << ", " << mom.py_ << ", " << mom.pz_ << ", " << mom.energy_ << ")\n\t"
      << "= " << energy_*mom.energy_-threeProduct(mom);
    return energy_*mom.energy_-threeProduct(mom);
  }

  double
  Momentum::crossProduct( const Momentum& mom ) const
  {
    return px_*mom.py_-py_*mom.px_;
  }

  double
  Momentum::operator*=( const Momentum& mom )
  {
    return threeProduct( mom );
  }

  Momentum&
  Momentum::operator*=( double c )
  {
    px_ *= c;
    py_ *= c;
    pz_ *= c;
    energy_ *= c;
    computeP();
    return *this;
  }

  Momentum
  operator*( const Momentum& mom, double c )
  {
    Momentum out = mom;
    out *= c;
    return out;
  }

  Momentum
  operator*( double c, const Momentum& mom )
  {
    Momentum out = mom;
    out *= c;
    return out;
  }

  Momentum
  operator+( const Momentum& mom1, const Momentum& mom2 )
  {
    Momentum out = mom1;
    out += mom2;
    return out;
  }

  Momentum
  operator-( const Momentum& mom1, const Momentum& mom2 )
  {
    Momentum out = mom1;
    out -= mom2;
    return out;
  }

  Momentum
  operator-( const Momentum& mom )
  {
    return Momentum()-mom;
  }

  double
  operator*( const Momentum& mom1, const Momentum& mom2 )
  {
    Momentum tmp = mom1;
    return tmp.threeProduct( mom2 );
  }

  //--- various setters

  void
  Momentum::setMass2( double m2 )
  {
    energy_ = sqrt( p2()+m2 );
  }

  void
  Momentum::setP( double px, double py, double pz, double e )
  {
    setP( px, py, pz );
    setEnergy( e );
  }

  void
  Momentum::setP( double px, double py, double pz )
  {
    px_ = px;
    py_ = py;
    pz_ = pz;
    computeP();
  }

  void
  Momentum::computeP()
  {
    p_ = std::hypot( pt(), pz_ );
  }

  void
  Momentum::truncate( double tolerance )
  {
    if ( px_ <= tolerance )
      px_ = 0.;
    if ( py_ <= tolerance )
      py_ = 0.;
    if ( pz_ <= tolerance )
      pz_ = 0.;
    if ( energy_ <= tolerance )
      energy_ = 0.;
    computeP();
  }

  //--- various getters

  double
  Momentum::operator[]( unsigned int i ) const
  {
    switch ( i ) {
      case 0: return px_;
      case 1: return py_;
      case 2: return pz_;
      case 3: return energy_;
      default:
        throw CG_FATAL( "Momentum" ) << "Failed to retrieve the component " << i << "!";
    }
  }

  double&
  Momentum::operator[]( const unsigned int i )
  {
    switch ( i ) {
      case 0: return px_;
      case 1: return py_;
      case 2: return pz_;
      case 3: return energy_;
      default:
        throw CG_FATAL( "Momentum" ) << "Failed to retrieve the component " << i << "!";
    }
  }

  const std::vector<double>
  Momentum::pVector() const
  {
    return std::vector<double>{ px(), py(), pz(), energy(), mass() };
  }

  double
  Momentum::mass() const
  {
    if ( mass2() >= 0. )
      return sqrt( mass2() );
    return -sqrt( -mass2() );
  }

  double
  Momentum::theta() const
  {
    return atan2( pt(), pz() );
  }

  double
  Momentum::phi() const
  {
    return atan2( py(), px() );
  }

  double
  Momentum::pt() const
  {
    return std::hypot( px_, py_ );
  }

  double
  Momentum::pt2() const
  {
    return px_*px_+py_*py_;
  }

  double
  Momentum::eta() const
  {
    const int sign = ( pz()/fabs( pz() ) );
    return ( pt() != 0. )
      ? log( ( p()+fabs( pz() ) )/pt() )*sign
      : 9999.*sign;
  }

  double
  Momentum::rapidity() const
  {
    const int sign = ( pz()/fabs( pz() ) );
    return ( energy() >= 0. )
      ? log( ( energy()+pz() )/( energy()-pz() ) )*0.5
      : 999.*sign;
  }

  //--- boosts/rotations

  Momentum&
  Momentum::betaGammaBoost( double gamma, double betagamma )
  {
    if ( gamma == 1. && betagamma == 0. ) return *this; // trivial case
    const double pz = pz_, e = energy_;
    pz_ = gamma*pz+betagamma*e;
    energy_  = gamma*e +betagamma*pz;
    computeP();
    return *this;
  }

  Momentum&
  Momentum::lorentzBoost( const Momentum& p )
  {
    //--- do not boost on a system at rest
    if ( p.p() == 0. )
      return *this;

    const double m = p.mass();
    const double pf4 = fourProduct( p )/m;
    const double fn = ( pf4+energy_ )/( p.energy_+m );
    *this -= p*fn;
    energy_ = pf4;
    return *this;
  }

  Momentum&
  Momentum::rotatePhi( double phi, double sign )
  {
    const double px = px_*cos( phi )+py_*sin( phi )*sign,
                 py =-px_*sin( phi )+py_*cos( phi )*sign;
    px_ = px;
    py_ = py;
    return *this;
  }

  Momentum&
  Momentum::rotateThetaPhi( double theta, double phi )
  {
    double rotmtx[3][3], mom[3]; //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
    rotmtx[0][0] = -sin( phi ); rotmtx[0][1] = -cos( theta )*cos( phi ); rotmtx[0][2] =  sin( theta )*cos( phi );
    rotmtx[1][0] =  cos( phi ); rotmtx[1][1] = -cos( theta )*sin( phi ); rotmtx[1][2] =  sin( theta )*sin( phi );
    rotmtx[2][0] =  0.;         rotmtx[2][1] =  sin( theta );            rotmtx[2][2] =  cos( theta );

    for ( unsigned short i = 0; i < 3; ++i ) {
      mom[i] = 0.;
      for ( unsigned short j = 0; j < 3; ++j )
        mom[i] += rotmtx[i][j]*operator[]( j );
    }
    setP( mom[0], mom[1], mom[2] );
    return *this;
  }

  //--- printout

  std::ostream&
  operator<<( std::ostream& os, const Momentum& mom )
  {
    return os << "(" << std::fixed
      << std::setw( 9 ) << mom.energy_ << "|"
      << std::setw( 9 ) << mom.px_ << " "
      << std::setw( 9 ) << mom.py_ << " "
      << std::setw( 9 ) << mom.pz_ << ")" << std::defaultfloat;
  }
}
