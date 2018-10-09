#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace particleproperties
  {
    double
    mass( const PDG& pdg_id )
    {
      switch ( pdg_id ) {
        case PDG::electron:     return 0.510998928e-3;
        case PDG::muon:         return 0.1056583715;
        case PDG::tau:          return 1.77682;
        case PDG::down:         return 0.0048;
        case PDG::up:           return 0.0023;
        case PDG::strange:      return 0.095;
        case PDG::charm:        return 1.29;
        case PDG::bottom:       return 4.18;
        case PDG::top:          return 172.44;
        case PDG::electronNeutrino:
        case PDG::muonNeutrino:
        case PDG::tauNeutrino:  return 0.;
        case PDG::gluon:
        case PDG::photon:
        case PDG::pomeron:      return 0.;
        case PDG::Z:            return 91.1876;
        case PDG::W:            return 80.385;
        case PDG::piPlus:       return 0.13957018;
        case PDG::piZero:       return 0.1349766;
        case PDG::KPlus:        return 0.49368;
        case PDG::DPlus:        return 1.86962;
        case PDG::Jpsi:         return 3.0969;
        case PDG::proton:       return 0.938272046;
        case PDG::neutron:      return 0.939565346;
        case PDG::Upsilon1S:    return 9.46030;
        case PDG::Upsilon2S:    return 10.02326;
        case PDG::Upsilon3S:    return 10.3552;
        case PDG::rho770_0:     return 0.77526;
        case PDG::rho1450_0:    return 1.465;
        case PDG::rho1700_0:    return 1.720;
        case PDG::h1380_1:      return 1.38619;
        case PDG::eta:          return 0.547862;
        case PDG::invalid:
        default:
          return -1.;
      }
    }

    double
    charge( const PDG& pdg_id )
    {
      switch ( pdg_id ) {
        case PDG::proton: case PDG::diffractiveProton:
          return +1.;
        case PDG::electron: case PDG::muon: case PDG::tau:
          return -1.;
        case PDG::down: case PDG::strange: case PDG::bottom:
          return -1./3;
        case PDG::up: case PDG::charm: case PDG::top:
          return +2./3;
        case PDG::W:
          return +1.;
        case PDG::piPlus: case PDG::KPlus: case PDG::DPlus:
          return +1.;
        default:
          return 0.;
      }
    }

    double
    charge( int id )
    {
      const short sign = id / abs( id );
      return sign * charge( (PDG)abs( id ) );
    }

    unsigned short
    colours( const PDG& pdg_id )
    {
      switch ( pdg_id ) {
        case PDG::top: return 3;
        default: return 1;
      }
    }

    double
    width( const PDG& pdg_id )
    {
      switch ( pdg_id ) {
        case PDG::Jpsi:      return 92.9e-6; //FIXME
        case PDG::Z:         return 2.4952;
        case PDG::W:         return 2.085;
        case PDG::Upsilon1S: return 54.02e-6;
        case PDG::Upsilon2S: return 31.98e-6;
        case PDG::Upsilon3S: return 20.32e-6;
        case PDG::rho770_0:  return 0.150; // PDG
        case PDG::rho1450_0: return 0.400; // PDG
        case PDG::rho1700_0: return 0.250; // PDG
        default:
          CG_WARNING( "ParticleProperties:width" )
            << "Particle " << pdg_id << " has no registered width.";
          return -1.;
      }
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const PDG& pc )
  {
    switch ( pc ) {
      case PDG::electron:         return os << "e± ";
      case PDG::electronNeutrino: return os << "ν_e ";
      case PDG::muon:             return os << "µ±  ";
      case PDG::muonNeutrino:     return os << "ν_µ  ";
      case PDG::tau:              return os << "τ±  ";
      case PDG::tauNeutrino:      return os << "ν_τ  ";
      case PDG::gluon:            return os << "gluon";
      case PDG::photon:           return os << "ɣ ";
      case PDG::Z:                return os << "Z";
      case PDG::W:                return os << "W± ";
      case PDG::piPlus:           return os << "π±  ";
      case PDG::piZero:           return os << "π⁰\t";
      case PDG::KPlus:            return os << "K± ";
      case PDG::DPlus:            return os << "D± ";
      case PDG::rho770_0:         return os << "ρ(770)₀  ";
      case PDG::rho1450_0:        return os << "ρ(1450)₀  ";
      case PDG::rho1700_0:        return os << "ρ(1700)₀  ";
      case PDG::h1380_1:          return os << "h(1380)₁ ";
      case PDG::eta:              return os << "η meson";
      case PDG::omega782:         return os << "ω(782) ";
      case PDG::Jpsi:             return os << "J/ψ ";
      case PDG::phi1680:          return os << "ɸ(1680) ";
      case PDG::Upsilon1S:        return os << "Υ(1S) ";
      case PDG::Upsilon2S:        return os << "Υ(2S) ";
      case PDG::Upsilon3S:        return os << "Υ(3S) ";;
      case PDG::proton:           return os << "proton";
      case PDG::diffractiveProton:return os << "diffr.proton";
      case PDG::neutron:          return os << "neutron";
      case PDG::pomeron:          return os << "IP";
      case PDG::reggeon:          return os << "IR";
      case PDG::down:             return os << "d";
      case PDG::up:               return os << "u";
      case PDG::strange:          return os << "s";
      case PDG::charm:            return os << "c";
      case PDG::bottom:           return os << "b";
      case PDG::top:              return os << "t";
      case PDG::invalid:          return os << "[...]";
    }
    return os;
  }
}
