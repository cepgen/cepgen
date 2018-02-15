#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <iostream>

namespace CepGen
{
  /** Unique identifier for a particle type. From \cite Beringer:1900zz :
   * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
   * \brief PDG ids of all known particles
   */
  enum ParticleCode {
    invalidParticle = 0,
    //--- fundamental particles
    dQuark = 1, uQuark = 2,
    Electron = 11, ElectronNeutrino = 12,
    Muon = 13, MuonNeutrino = 14,
    Tau = 15, TauNeutrino = 16,
    Gluon = 21, Photon = 22, Z = 23, W = 24,
    //--- composite particles
    PiPlus = 211, PiZero = 111,
    Rho770_0 = 113, Rho1450_0 = 100113, Rho1700_0 = 30113,
    Eta = 221, Omega782 = 223,
    h1380_1 = 10333,
    JPsi= 443,
    Phi1680 = 100333,
    Upsilon1S = 553, Upsilon2S = 100553, Upsilon3S = 200553,
    Proton = 2212, Neutron = 2112,
    Pomeron = 990, Reggeon = 110
  };
  namespace ParticleProperties
  {
    /**
     * \brief Gets the mass of a particle
     * \param pdgId ParticleCode (PDG ID)
     * \return Mass of the particle in \f$\textrm{GeV}/c^2\f$
     */
    inline double mass( const ParticleCode& pdgId ) {
      switch ( pdgId ) {
        case dQuark:       return 0.33;           // mass from PYTHIA6.4
        case uQuark:       return 0.33;           // mass from PYTHIA6.4
        case Electron:     return 0.510998928e-3;
        case Muon:         return 0.1056583715;
        case Tau:          return 1.77682;
        case ElectronNeutrino: case MuonNeutrino: case TauNeutrino: return 0.;
        case Gluon: case Photon: return 0.;
        case Z:            return 91.1876;
        case W:            return 80.385;
        case PiPlus:       return 0.13957018;
        case PiZero:       return 0.1349766;
        case JPsi:         return 20.;            //FIXME FIXME FIXME
        case Proton:       return 0.938272046;
        case Neutron:      return 0.939565346;
        case Upsilon1S:    return 9.46030;
        case Upsilon2S:    return 10.02326;
        case Upsilon3S:    return 10.3552;
        case Rho770_0:     return 0.77526;
        case Rho1450_0:    return 1.465;
        case Rho1700_0:    return 1.720;
        case h1380_1:      return 1.38619;
        case Eta:          return 0.547862;
        case invalidParticle:
        default:           return -1.;
      }
    }
    /**
     * \brief Gets the electric charge of a particle
     * \param pdgId ParticleCode (PDG ID)
     * \return Charge of the particle in \f$e\f$
     */
    inline double charge( const ParticleCode& pdgId ) {
      switch ( pdgId ) {
        case Proton: return +1.;
        case dQuark: return -1./3;
        case uQuark: return +2./3;
        case Electron: case Muon: case Tau: return -1.;
        case ElectronNeutrino: case MuonNeutrino: case TauNeutrino: return 0.;
        case Gluon: case Z: case Photon: return 0.;
        case W: return +1.;
        case PiPlus: return +1.;
        case PiZero: return 0.;
        case Neutron: return 0.;
        case Eta: return 0.;
        default: return 0.;
      }
    }
    /**
     * \brief Total decay width of one unstable particle
     * \param[in] pdgId ParticleCode (PDG ID)
     * \return Decay width in GeV
     */
    inline double width( const ParticleCode& pdgId ) {
      switch ( pdgId ) {
        case JPsi:      return 5.; //FIXME
        case Z:         return 2.4952;
        case W:         return 2.085;
        case Upsilon1S: return 54.02e-6;
        case Upsilon2S: return 31.98e-6;
        case Upsilon3S: return 20.32e-6;
        case Rho770_0:  return 0.150; // PDG
        case Rho1450_0: return 0.400; // PDG
        case Rho1700_0: return 0.250; // PDG
        default:        return -1.;
      }
    }
  }

  inline std::ostream& operator<<( std::ostream& os, const ParticleCode& pc ) {
    switch ( pc ) {
      case dQuark:           return os << "d quark";
      case uQuark:           return os << "u quark";
      case Electron:         return os << "electron";
      case ElectronNeutrino: return os << "nu_e";
      case Muon:             return os << "muon";
      case MuonNeutrino:     return os << "nu_mu";
      case Tau:              return os << "tau";
      case TauNeutrino:      return os << "nu_tau";
      case Gluon:            return os << "gluon";
      case Photon:           return os << "photon";
      case Z:                return os << "Z";
      case W:                return os << "W+-";
      case PiPlus:           return os << "pi+";
      case PiZero:           return os << "pi0";
      case Rho770_0:         return os << "rho(770)0";
      case Rho1450_0:        return os << "rho(1450)0";
      case Rho1700_0:        return os << "rho(1700)0";
      case h1380_1:          return os << "h(1380)1";
      case Eta:              return os << "eta meson";
      case Omega782:         return os << "omega(782)";
      case JPsi:             return os << "J/Psi";
      case Phi1680:          return os << "phi(1680)";
      case Upsilon1S:        return os << "Upsilon(1S)";
      case Upsilon2S:        return os << "Upsilon(2S)";
      case Upsilon3S:        return os << "Upsilon(3S)";;
      case Proton:           return os << "proton";
      case Neutron:          return os << "neutron";
      case Pomeron:          return os << "pomeron";
      case Reggeon:          return os << "reggeon";
      case invalidParticle:  return os << "[...]";
    }
    return os;
  }
}

#endif
