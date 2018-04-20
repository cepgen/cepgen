#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <iostream>

namespace CepGen
{
  /** Unique identifier for a particle type. From \cite Beringer:1900zz :
   * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
   * \brief PDG ids of all known particles
   */
  enum class PDG {
    invalid = 0,
    //--- fundamental particles
    TopQuark = 6,
    Electron = 11, ElectronNeutrino = 12,
    Muon = 13, MuonNeutrino = 14,
    Tau = 15, TauNeutrino = 16,
    Gluon = 21, Photon = 22, Z = 23, W = 24,
    //--- composite particles
    PiPlus = 211, PiZero = 111,
    KPlus = 321, DPlus = 411,
    Rho770_0 = 113, Rho1450_0 = 100113, Rho1700_0 = 30113,
    Eta = 221, Omega782 = 223,
    h1380_1 = 10333,
    JPsi= 443,
    Phi1680 = 100333,
    Upsilon1S = 553, Upsilon2S = 100553, Upsilon3S = 200553,
    Proton = 2212, Neutron = 2112,
    Pomeron = 990, Reggeon = 110,
    DiffractiveProton = 9902210
  };
  struct PDGHash
  {
    template <typename T> std::size_t operator()( T t ) const {
      return static_cast<std::size_t>( t );
    }
  };
  namespace ParticleProperties
  {
    /**
     * \brief Gets the mass of a particle
     * \param pdgId PDG identifier
     * \return Mass of the particle in \f$\textrm{GeV}/c^2\f$
     */
    inline double mass( const PDG& pdgId ) {
      switch ( pdgId ) {
        case PDG::Electron:     return 0.510998928e-3;
        case PDG::Muon:         return 0.1056583715;
        case PDG::Tau:          return 1.77682;
        case PDG::TopQuark:     return 172.44;
        case PDG::ElectronNeutrino:
        case PDG::MuonNeutrino:
        case PDG::TauNeutrino:  return 0.;
        case PDG::Gluon:
        case PDG::Photon:       return 0.;
        case PDG::Z:            return 91.1876;
        case PDG::W:            return 80.385;
        case PDG::PiPlus:       return 0.13957018;
        case PDG::PiZero:       return 0.1349766;
        case PDG::KPlus:        return 0.49368;
        case PDG::DPlus:        return 1.86962;
        case PDG::JPsi:         return 20.;            //FIXME FIXME FIXME
        case PDG::Proton:       return 0.938272046;
        case PDG::Neutron:      return 0.939565346;
        case PDG::Upsilon1S:    return 9.46030;
        case PDG::Upsilon2S:    return 10.02326;
        case PDG::Upsilon3S:    return 10.3552;
        case PDG::Rho770_0:     return 0.77526;
        case PDG::Rho1450_0:    return 1.465;
        case PDG::Rho1700_0:    return 1.720;
        case PDG::h1380_1:      return 1.38619;
        case PDG::Eta:          return 0.547862;
        case PDG::invalid:
        default:           return -1.;
      }
    }
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] pdg_id PDG id
    inline double charge( const PDG& pdg_id ) {
      switch ( pdg_id ) {
        case PDG::Proton:
        case PDG::DiffractiveProton: return +1.;
        case PDG::Electron:
        case PDG::Muon:
        case PDG::Tau:               return -1.;
        case PDG::TopQuark:          return +2./3;
        case PDG::W:                 return +1.;
        case PDG::PiPlus:
        case PDG::KPlus:
        case PDG::DPlus:             return +1.;
        default: return 0.;
      }
    }
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] id integer PDG id
    inline double charge( int id ) {
      const short sign = id / abs( id );
      return sign * charge( (PDG)abs( id ) );
    }
    /// Total decay width of an unstable particle, in GeV
    /// \param[in] pdg_id PDG (PDG ID)
    inline double width( const PDG& pdg_id ) {
      switch ( pdg_id ) {
        case PDG::JPsi:      return 5.; //FIXME
        case PDG::Z:         return 2.4952;
        case PDG::W:         return 2.085;
        case PDG::Upsilon1S: return 54.02e-6;
        case PDG::Upsilon2S: return 31.98e-6;
        case PDG::Upsilon3S: return 20.32e-6;
        case PDG::Rho770_0:  return 0.150; // PDG
        case PDG::Rho1450_0: return 0.400; // PDG
        case PDG::Rho1700_0: return 0.250; // PDG
        default:             return -1.;
      }
    }
  }

  inline std::ostream& operator<<( std::ostream& os, const PDG& pc ) {
    switch ( pc ) {
      case PDG::Electron:         return os << "e± ";
      case PDG::ElectronNeutrino: return os << "ν_e ";
      case PDG::Muon:             return os << "µ±  ";
      case PDG::MuonNeutrino:     return os << "ν_µ  ";
      case PDG::Tau:              return os << "τ±  ";
      case PDG::TauNeutrino:      return os << "ν_τ  ";
      case PDG::Gluon:            return os << "gluon";
      case PDG::Photon:           return os << "ɣ ";
      case PDG::Z:                return os << "Z";
      case PDG::W:                return os << "W± ";
      case PDG::PiPlus:           return os << "π±  ";
      case PDG::PiZero:           return os << "π⁰\t";
      case PDG::KPlus:            return os << "K± ";
      case PDG::DPlus:            return os << "D± ";
      case PDG::Rho770_0:         return os << "ρ(770)₀  ";
      case PDG::Rho1450_0:        return os << "ρ(1450)₀  ";
      case PDG::Rho1700_0:        return os << "ρ(1700)₀  ";
      case PDG::h1380_1:          return os << "h(1380)₁ ";
      case PDG::Eta:              return os << "η meson";
      case PDG::Omega782:         return os << "ω(782) ";
      case PDG::JPsi:             return os << "J/ψ ";
      case PDG::Phi1680:          return os << "ɸ(1680) ";
      case PDG::Upsilon1S:        return os << "Υ(1S) ";
      case PDG::Upsilon2S:        return os << "Υ(2S) ";
      case PDG::Upsilon3S:        return os << "Υ(3S) ";;
      case PDG::Proton:           return os << "proton";
      case PDG::DiffractiveProton:return os << "diffr.proton";
      case PDG::Neutron:          return os << "neutron";
      case PDG::Pomeron:          return os << "pomeron";
      case PDG::Reggeon:          return os << "reggeon";
      case PDG::TopQuark:         return os << "t";
      case PDG::invalid:          return os << "[...]";
    }
    return os;
  }
}

#endif
