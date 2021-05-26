#include "CepGen/Event/Particle.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <iomanip>

namespace cepgen {
  Particle::Particle()
      : id_(-1),
        charge_sign_(1),
        mass_(-1.),
        helicity_(0.),
        role_(UnknownRole),
        status_((int)Status::Undefined),
        pdg_id_(PDG::invalid) {}

  Particle::Particle(Role role, pdgid_t pdgId, Status st)
      : id_(-1), charge_sign_(1), mass_(-1.), helicity_(0.), role_(role), status_((int)st), pdg_id_(pdgId) {
    try {
      phys_prop_ = PDG::get()(pdg_id_);
    } catch (const Exception&) {
    }
    if (pdg_id_ != PDG::invalid)
      computeMass();
  }

  Particle::Particle(const Particle& part)
      : id_(part.id_),
        charge_sign_(part.charge_sign_),
        momentum_(part.momentum_),
        mass_(part.mass_),
        helicity_(part.helicity_),
        role_(part.role_),
        status_(part.status_),
        mothers_(part.mothers_),
        daughters_(part.daughters_),
        pdg_id_(part.pdg_id_) {
    try {
      phys_prop_ = PDG::get()(pdg_id_);
    } catch (const Exception&) {
    }
  }

  bool Particle::operator<(const Particle& rhs) const { return id_ >= 0 && rhs.id_ > 0 && id_ < rhs.id_; }

  double Particle::thetaToEta(double theta) { return -log(tan(0.5 * theta * M_PI / 180.)); }

  double Particle::etaToTheta(double eta) { return 2. * atan(exp(-eta)) * 180. * M_1_PI; }

  bool Particle::valid() {
    if (pdg_id_ == PDG::invalid)
      return false;
    if (momentum_.p() == 0. && mass_ == 0.)
      return false;
    return true;
  }

  float Particle::charge() const { return charge_sign_ * phys_prop_.charge / 3.; }

  Particle& Particle::computeMass(bool off_shell) {
    if (!off_shell && pdg_id_ != PDG::invalid)  // retrieve the mass from the on-shell particle's properties
      mass_ = phys_prop_.mass;
    else if (momentum_.energy() >= 0.)
      mass_ = sqrt(energy2() - momentum_.p2());
    //--- finish by setting the energy accordingly
    if (momentum_.energy() < 0.)  // invalid energy
      momentum_.setEnergy(sqrt(momentum_.p2() + mass2()));

    return *this;
  }

  Particle& Particle::setMass(double m) {
    if (m < 0.)
      return computeMass();
    mass_ = m;
    return *this;
  }

  Particle& Particle::clearMothers() {
    mothers_.clear();
    return *this;
  }

  Particle& Particle::addMother(Particle& part) {
    mothers_.insert(part.id());

    CG_DEBUG_LOOP("Particle") << "Particle " << id() << " (pdgId=" << part.integerPdgId() << ") "
                              << "is the new mother of " << id_ << " (pdgId=" << (int)pdg_id_ << ").";

    part.addDaughter(*this);
    return *this;
  }

  Particle& Particle::clearDaughters() {
    daughters_.clear();
    return *this;
  }

  Particle& Particle::addDaughter(Particle& part) {
    const auto ret = daughters_.insert(part.id());

    CG_DEBUG_LOOP("Particle").log([&](auto& dbg) {
      dbg << "Particle " << role_ << " (pdgId=" << (int)pdg_id_ << ")"
          << " has now " << utils::s("daughter", daughters_.size(), true) << ":";
      for (const auto& daughter : daughters_)
        dbg << utils::format("\n\t * id=%d", daughter);
    });

    if (ret.second) {
      CG_DEBUG_LOOP("Particle") << "Particle " << part.role() << " (pdgId=" << part.integerPdgId() << ") "
                                << "is a new daughter of " << role_ << " (pdgId=" << pdg_id_ << ").";

      if (!utils::contains(part.mothers(), id_))
        part.addMother(*this);
    }
    return *this;
  }

  Particle& Particle::setMomentum(const Momentum& mom, bool offshell) {
    momentum_ = mom;
    if (offshell || mom.mass() <= 0.)
      return computeMass();
    mass_ = momentum_.mass();
    return *this;
  }

  Particle& Particle::setMomentum(double px, double py, double pz, double e) {
    momentum_.setP(px, py, pz);
    setEnergy(e);
    if (fabs(e - momentum_.energy()) > 1.e-6)  // more than 1 eV difference
      CG_WARNING("Particle") << "Energy difference: " << e - momentum_.energy();
    return *this;
  }

  double Particle::energy() const {
    return (momentum_.energy() < 0. ? std::hypot(mass_, momentum_.p()) : momentum_.energy());
  }

  Particle& Particle::setEnergy(double e) {
    if (e < 0. && mass_ >= 0.)
      e = std::hypot(mass_, momentum_.p());
    momentum_.setEnergy(e);
    return *this;
  }

  pdgid_t Particle::pdgId() const { return pdg_id_; }

  Particle& Particle::setPdgId(long pdg) {
    pdg_id_ = labs(pdg);
    try {
      phys_prop_ = PDG::get()(pdg_id_);
    } catch (const Exception&) {
    }
    switch (pdg_id_) {
      case PDG::electron:
      case PDG::muon:
      case PDG::tau:
        charge_sign_ = -pdg / labs(pdg);
        break;
      default:
        charge_sign_ = pdg / labs(pdg);
        break;
    }
    return *this;
  }

  Particle& Particle::setPdgId(pdgid_t pdg, short ch) {
    pdg_id_ = pdg;
    switch (pdg_id_) {
      case PDG::electron:
      case PDG::muon:
      case PDG::tau:
        charge_sign_ = -ch;
        break;
      default:
        charge_sign_ = ch;
        break;
    }
    return *this;
  }

  int Particle::integerPdgId() const {
    const float ch = phys_prop_.charge / 3.;
    if (ch == 0)
      return static_cast<int>(pdg_id_);
    return static_cast<int>(pdg_id_) * charge_sign_ * (ch / fabs(ch));
  }

  double Particle::etaToY(double eta_, double m_, double pt_) {
    const double m2 = m_ * m_, mt = std::hypot(m_, pt_);
    return asinh(sqrt(((mt * mt - m2) * cosh(2. * eta_) + m2) / mt * mt - 1.) * M_SQRT1_2);
  }

  std::ostream& operator<<(std::ostream& os, const Particle& part) {
    os << std::resetiosflags(std::ios::showbase) << "Particle[" << part.id_ << "]{role=" << part.role_
       << ", status=" << (int)part.status_ << ", "
       << "pdg=" << part.integerPdgId() << ", p4=" << part.momentum_ << " GeV, m=" << part.mass_ << " GeV, "
       << "p⟂=" << part.momentum_.pt() << " GeV, eta=" << part.momentum_.eta() << ", phi=" << part.momentum_.phi();
    if (part.primary())
      os << ", primary";
    else {
      os << ", " << utils::s("mother", part.mothers_.size()) << "=";
      std::string delim;
      for (const auto& moth : part.mothers_)
        os << delim << moth, delim = ",";
    }
    const auto& daughters_list = part.daughters();
    if (!daughters_list.empty()) {
      os << ", " << utils::s("daughter", daughters_list.size()) << "=";
      std::string delim;
      for (const auto& daughter : daughters_list)
        os << delim << daughter, delim = ",";
    }
    return os << "}";
  }

  std::ostream& operator<<(std::ostream& os, const Particle::Role& rl) {
    switch (rl) {
      case Particle::UnknownRole:
        return os << "unknown";
      case Particle::IncomingBeam1:
        return os << "i.beam 1";
      case Particle::IncomingBeam2:
        return os << "i.beam 2";
      case Particle::OutgoingBeam1:
        return os << "o.beam 1";
      case Particle::OutgoingBeam2:
        return os << "o.beam 2";
      case Particle::Parton1:
        return os << "parton 1";
      case Particle::Parton2:
        return os << "parton 2";
      case Particle::Intermediate:
        return os << "hard pr.";
      case Particle::CentralSystem:
        return os << "central";
    }
    return os;
  }

  double CMEnergy(const Particle& p1, const Particle& p2) {
    if (p1.mass() * p2.mass() < 0. || p1.energy() * p2.energy() < 0.)
      return 0.;
    return sqrt(p1.mass2() + p2.mass2() + 2. * p1.energy() * p2.energy() - 2. * (p1.momentum() * p2.momentum()));
  }
}  // namespace cepgen
