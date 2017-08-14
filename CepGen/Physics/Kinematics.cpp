#include "Kinematics.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    in1p( 6500. ), in2p( 6500. ), in1pdg( Particle::Proton ), in2pdg( Particle::Proton ),
    pair( Particle::Muon ),
    kinematics( ElasticElastic ), remnant_mode( SuriYennie ), cuts_mode( BothParticles ),
    pt_min( 3. ), pt_max( -1. ), e_min( 0. ), e_max( -1. ), eta_min( -999. ), eta_max( 999. ),
    mass_min( 0. ), mass_max( -1. ),
    mx_min( 1.07 ), mx_max( 320. ),
    q2_min( 0. ), q2_max( 1.e5 ), w_min( 0. ), w_max( -1. ),
    ptdiff_min( 0. ), ptdiff_max( 300. ),
    qt_min( 0. ), qt_max( 500. )
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
      << std::setw(25) << "Minimal pT :" << std::setw(8) << pt_min << std::endl
      << std::setw(25) << "Maximal pT :" << std::setw(8) << pt_max << std::endl
      << std::setw(25) << "Minimal energy :" << std::setw(8) << e_min << std::endl
      << std::setw(25) << "Maximal energy :" << std::setw(8) << e_max << std::endl
      << std::setw(25) << "Minimal pseudorapidity :" << std::setw(8) << eta_min << std::endl
      << std::setw(25) << "Maximal pseudorapidity :" << std::setw(8) << eta_max << std::endl
      << "===== Central kinematics" << std::endl
      << std::setw(25) << "Minimal Q**2 :" << std::setw(8) << q2_min << std::endl
      << std::setw(25) << "Maximal Q**2 :" << std::setw(8) << q2_max << std::endl
      << std::setw(25) << "Minimal W :" << std::setw(8) << w_min << std::endl
      << std::setw(25) << "Maximal W :" << std::setw(8) << w_max << std::endl;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::ProcessMode& pm )
  {
    switch ( pm ) {
      case Kinematics::ElectronElectron:    os << "electron/electron"; break;
      case Kinematics::ElectronProton:      os << "electron/proton"; break;
      case Kinematics::ProtonElectron:      os << "proton/electron"; break;
      case Kinematics::ElasticElastic:      os << "elastic/elastic"; break;
      case Kinematics::InelasticElastic:    os << "inelastic/elastic"; break;
      case Kinematics::ElasticInelastic:    os << "elastic/inelastic"; break;
      case Kinematics::InelasticInelastic:  os << "inelastic/inelastic"; break;    
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Cuts& cut )
  {
    switch ( cut ) {
      case Kinematics::NoCuts:         os << "no cuts"; break;
      case Kinematics::BothParticles:  os << "both outgoing particles"; break;
      case Kinematics::OneParticle:    os << "single outgoing particle"; break;
    }
    return os;
  }
}

