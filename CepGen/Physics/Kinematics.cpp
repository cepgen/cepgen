#include "Kinematics.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Particle.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    inp( { 6500., 6500. } ), inpdg( { PDG::Proton, PDG::Proton } ),
    mode( Mode::ElasticElastic ), structure_functions( StructureFunctions::SuriYennie )
  {}

  Kinematics::Kinematics( const Kinematics& kin ) :
    inp( kin.inp ), inpdg( kin.inpdg ),
    central_system( kin.central_system ), minimum_final_state( kin.minimum_final_state ),
    mode( kin.mode ), structure_functions( kin.structure_functions ),
    cuts( kin.cuts )
  {}

  Kinematics&
  Kinematics::operator=( const Kinematics& kin )
  {
    inp = kin.inp; inpdg = kin.inpdg;
    central_system = kin.central_system;
    minimum_final_state = kin.minimum_final_state;
    mode = kin.mode;
    structure_functions = kin.structure_functions;
    cuts = kin.cuts;
    return *this;
  }

  Kinematics::~Kinematics()
  {}

  void
  Kinematics::setSqrtS( double sqrts )
  {
    const double pin = 0.5 * sqrts;
    inp = { pin, pin };
  }

  double
  Kinematics::sqrtS() const
  {
    return inp.first + inp.second;
  }

  void
  Kinematics::dump( std::ostream& os ) const
  {
    os << std::setfill(' ');
    os << "===== Central system\n";
    for ( const auto& pdg_lim : cuts.central )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;

    os << "===== Initial state\n";
    for ( const auto& pdg_lim : cuts.initial )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;

    os << "===== Remnants\n";
    for ( const auto& pdg_lim : cuts.remnants )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Mode& pm )
  {
    switch ( pm ) {
      case Kinematics::Mode::ElectronElectron:
        return os << "electron/electron";
      case Kinematics::Mode::ElectronProton:
        return os << "electron/proton";
      case Kinematics::Mode::ProtonElectron:
        return os << "proton/electron";
      case Kinematics::Mode::ElasticElastic:
        return os << "elastic/elastic";
      case Kinematics::Mode::InelasticElastic:
        return os << "inelastic/elastic";
      case Kinematics::Mode::ElasticInelastic:
        return os << "elastic/inelastic";
      case Kinematics::Mode::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }

  //------------------------------------------------------------------------------------------------
  // Kinematics limits
  //------------------------------------------------------------------------------------------------

  Kinematics::Limits::Limits( double min, double max ) :
    std::pair<double,double>( min, max )
  {}

  Kinematics::Limits::Limits( const Limits& rhs ) :
    std::pair<double,double>( rhs.first, rhs.second )
  {}

  void
  Kinematics::Limits::in( double low, double up )
  {
    first = low;
    second = up;
  }

  double
  Kinematics::Limits::range() const
  {
    if ( !hasMin() || !hasMax() )
      return 0.;
    return second-first;
  }

  bool
  Kinematics::Limits::hasMin() const
  {
    return first != kInvalid;
  }

  bool
  Kinematics::Limits::hasMax() const
  {
    return second != kInvalid;
  }

  bool
  Kinematics::Limits::passes( double val ) const
  {
    if ( hasMin() && val < min() )
      return false;
    if ( hasMax() && val > max() )
      return false;
    return true;
  }

  bool
  Kinematics::Limits::valid() const
  {
    return hasMin() || hasMax();
  }

  double
  Kinematics::Limits::x( double v ) const
  {
    if ( v < 0. || v > 1. )
      CG_ERROR( "Kinematics:shoot" )
        << "x must be comprised between 0 and 1; x value = " << v << ".";
    if ( !valid() )
      return kInvalid;

    return first + ( second-first ) * v;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Limits& lim )
  {
    if ( !lim.hasMin() && !lim.hasMax() )
      return os << "no cuts";
    if ( !lim.hasMin() )
      return os << Form( "≤ %g", lim.max() );
    if ( !lim.hasMax() )
      return os << Form( "≥ %g", lim.min() );
    return os << Form( "%g → %g", lim.min(), lim.max() );
  }

  //------------------------------------------------------------------------------------------------
  // List of kinematics limits
  //------------------------------------------------------------------------------------------------

  Kinematics::CutsList::CutsList() :
    initial( { { Cuts::q2, { 0., 1.e5 } } } ),
    central( { { Cuts::pt_single, 0. } } ),
    remnants( { { Cuts::mass_single, { 1.07, 320. } } } )
  {}

  Kinematics::CutsList::CutsList( const CutsList& cuts ) :
    initial( cuts.initial ),
    central( cuts.central ), central_particles( cuts.central_particles ),
    remnants( cuts.remnants )
  {}

  Kinematics::CutsList&
  Kinematics::CutsList::operator=( const CutsList& cuts )
  {
    initial = cuts.initial;
    central = cuts.central;
    central_particles = cuts.central_particles;
    remnants = cuts.remnants;
    return *this;
  }
}

