#include "Kinematics.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    inp( { 6500., 6500. } ), inpdg( { Particle::Proton, Particle::Proton } ),
    central_system( {} ),
    mode( ElasticElastic ), structure_functions( SuriYennie ),
    central_cuts( { { Cuts::pt_single, 3.0 }, { Cuts::pt_diff, { 0., 400.0 } } } ),
    remnant_cuts( { { Cuts::mass, { 1.07, 320.0 } } } ),
    initial_cuts( { { Cuts::q2, { 0.0, 1.0e5 } }, { Cuts::qt, { 0.0, 500.0 } } } )
  {}

  Kinematics::~Kinematics()
  {}

  void
  Kinematics::dump( std::ostream& os ) const
  {
    os << std::setfill(' ');
    os << "===== Central system\n";
    for ( std::map<Cuts::Central,Limits>::const_iterator lim = central_cuts.begin(); lim != central_cuts.end(); ++lim ) {
      os << std::setw(30) << lim->first << ": " << lim->second;
    }
    os << "===== Initial state\n";
    for ( std::map<Cuts::InitialState,Limits>::const_iterator lim = initial_cuts.begin(); lim != initial_cuts.end(); ++lim ) {
      os << std::setw(30) << lim->first << ": " << lim->second;
    }
    os << "===== Remnants\n";
    for ( std::map<Cuts::Remnants,Limits>::const_iterator lim = remnant_cuts.begin(); lim != remnant_cuts.end(); ++lim ) {
      os << std::setw(30) << lim->first << ": " << lim->second;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::ProcessMode& pm )
  {
    switch ( pm ) {
      case Kinematics::ElectronElectron:   return os << "electron/electron";
      case Kinematics::ElectronProton:     return os << "electron/proton";
      case Kinematics::ProtonElectron:     return os << "proton/electron";
      case Kinematics::ElasticElastic:     return os << "elastic/elastic";
      case Kinematics::InelasticElastic:   return os << "inelastic/elastic";
      case Kinematics::ElasticInelastic:   return os << "elastic/inelastic";
      case Kinematics::InelasticInelastic: return os << "inelastic/inelastic";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Limits& lim )
  {
    if ( !lim.hasMin() && !lim.hasMax() ) return os << "no cuts";
    if ( !lim.hasMin() ) return os << Form( "≤ %.3f", lim.max() );
    if ( !lim.hasMax() ) return os << Form( "≥ %.3f", lim.min() );
    return os << Form( "%.3f → %.3f", lim.min(), lim.max() );
  }
}

