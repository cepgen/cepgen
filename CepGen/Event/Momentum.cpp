#include "Particle.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <math.h>

namespace CepGen
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

  //--- arithmetic operators

  Momentum&
  Momentum::operator+=( const Momentum& mom )
  {
    px_ += mom.px_;
    py_ += mom.py_;
    pz_ += mom.pz_;
    energy_ += mom.energy_; //FIXME not supposed to be this way!
    computeP();
    return *this;
  }

  Momentum&
  Momentum::operator-=( const Momentum& mom )
  {
    px_ -= mom.px_;
    py_ -= mom.py_;
    pz_ -= mom.pz_;
    energy_ -= mom.energy_; //FIXME not supposed to be this way!
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
    DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
      px_, py_, pz_, energy_,
      mom.px_, mom.py_, mom.pz_, mom.energy_,
      px_*mom.px_+py_*mom.py_+pz_*mom.pz_ ) );
    return px_*mom.px_+py_*mom.py_+pz_*mom.pz_;
  }

  double
  Momentum::fourProduct( const Momentum& mom ) const
  {
    DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
      px_, py_, pz_, energy_,
      mom.px_, mom.py_, mom.pz_, mom.energy_,
      px_*mom.px_+py_*mom.py_+pz_*mom.pz_ ) );
    return energy_*mom.energy_-threeProduct(mom);
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
  Momentum::setP( unsigned int i, double p )
  {
    switch ( i ) {
      case 0: px_ = p; break;
      case 1: py_ = p; break;
      case 2: pz_ = p; break;
      case 3: energy_ = p; break;
      default: return;
    }
    computeP();
  }

  void
  Momentum::computeP()
  {
    p_ = sqrt( px_*px_ + py_*py_ + pz_*pz_ );
  }

  //--- various getters

  double
  Momentum::operator[]( const unsigned int i ) const
  {
    switch ( i ) {
      case 0: return px_;
      case 1: return py_;
      case 2: return pz_;
      case 3: return energy_;
      default:
        throw Exception( __PRETTY_FUNCTION__, Form( "Failed to retrieve the component %d", i ), FatalError );
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
        throw Exception( __PRETTY_FUNCTION__, Form( "Failed to retrieve the component %d", i ), FatalError );
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
    return sqrt( pt2() );
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
    const double m = p.mass();
    if ( m == p.energy() )
      return *this;

    const double pf4 = ( ( *this )*p ) / m,
                 fn = ( pf4+energy() )/( p.energy()+m );
    *this -= p*fn;
    setEnergy( pf4 );
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
    return os << "(E ; p) = (" << mom.energy_ << " ; " << mom.px_ << ", " << mom.py_ << ", " << mom.pz_ << ")";
  }
}
