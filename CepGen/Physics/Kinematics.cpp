#include "Kinematics.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    in1p( 6500. ), in2p( 6500. ), in1pdg( Particle::Proton ), in2pdg( Particle::Proton ),
    pair( Particle::Muon ),
    mode( ElasticElastic ), remnant_mode( SuriYennie ),
    cuts_mode( BothParticles ),
    pt_single_central( 3. ), mass_remnants( 1.07, 320. ),
    q2( 0., 1.e5 ),
    pt_diff_central( 0., 300. ), qt( 0., 500 )
  {}

  Kinematics::~Kinematics()
  {}

  void
  Kinematics::dump( std::ostream& os )
  {
    os
      << std::setfill(' ')
      << __PRETTY_FUNCTION__ << " Dump" << std::endl
      << std::setw(25) << "Cuts mode :" << std::setw(2) << cuts_mode << "->" << std::setw(4) << cuts_mode << std::endl    
      << "===== Single leptons" << std::endl
      << std::setw(25) << "pT range: " << pt_single_central << std::endl
      << std::setw(25) << "Energy range: " << e_single_central << std::endl
      << std::setw(25) << "Pseudorapidity range:" << eta_single_central << std::endl
      << "===== Central kinematics" << std::endl
      << std::setw(25) << "Q**2 range:" << q2 << std::endl
      << std::setw(25) << "W range:" << w << std::endl;
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
  operator<<( std::ostream& os, const Kinematics::Cuts& cut )
  {
    switch ( cut ) {
      case Kinematics::NoCuts:         return os << "no cuts";
      case Kinematics::BothParticles:  return os << "both outgoing particles";
      case Kinematics::OneParticle:    return os << "single outgoing particle";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Limits& lim )
  {
    return os << ( ( !lim.hasUpper() ) ? Form( ">= %.3f", lim.lower() ) : Form( "%.3f → %.3f", lim.lower(), lim.upper() ) );
  }
}

