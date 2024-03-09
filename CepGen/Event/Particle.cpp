/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iomanip>

#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  Particle::Particle(Role role, pdgid_t pdgId, Status st) : role_(role), status_((int)st), pdg_id_(pdgId) {
    if (PDG::get().has(pdg_id_))
      phys_prop_ = PDG::get()(pdg_id_);
  }

  Particle::Particle(const Particle& part)
      : id_(part.id_),
        charge_sign_(part.charge_sign_),
        momentum_(part.momentum_),
        helicity_(part.helicity_),
        role_(part.role_),
        status_(part.status_),
        mothers_(part.mothers_),
        daughters_(part.daughters_),
        pdg_id_(part.pdg_id_) {
    if (PDG::get().has(pdg_id_))
      phys_prop_ = PDG::get()(pdg_id_);
  }

  bool Particle::operator<(const Particle& rhs) const { return id_ >= 0 && rhs.id_ > 0 && id_ < rhs.id_; }

  bool Particle::valid() {
    if (pdg_id_ == PDG::invalid)
      return false;
    if (momentum_.p() == 0. && momentum_.mass() == 0.)
      return false;
    return true;
  }

  float Particle::charge() const { return charge_sign_ * phys_prop_.integerCharge() / 3.; }

  Particle& Particle::clearMothers() {
    mothers_.clear();
    return *this;
  }

  Particle& Particle::addMother(Particle& part) {
    const auto ret = mothers_.insert(part.id_);
    if (part.status_ > 0)
      part.status_ = (int)Status::Propagator;

    if (ret.second) {
      CG_DEBUG_LOOP("Particle") << "Particle " << id() << " (pdgId=" << part.pdg_id_ << ") "
                                << "is a new mother of " << id_ << " (pdgId=" << pdg_id_ << ").";
      if (!utils::contains(part.daughters_, id_))
        part.addDaughter(*this);
    }
    return *this;
  }

  Particle& Particle::clearDaughters() {
    daughters_.clear();
    return *this;
  }

  Particle& Particle::addDaughter(Particle& part) {
    const auto ret = daughters_.insert(part.id_);
    if (status_ > 0)
      status_ = (int)Status::Propagator;

    CG_DEBUG_LOOP("Particle").log([&](auto& dbg) {
      dbg << "Particle " << role_ << " (pdgId=" << (int)pdg_id_ << ")"
          << " has now " << utils::s("daughter", daughters_.size(), true) << ":";
      for (const auto& daughter : daughters_)
        dbg << utils::format("\n\t * id=%d", daughter);
    });

    if (ret.second) {
      CG_DEBUG_LOOP("Particle") << "Particle " << part.role_ << " (pdgId=" << part.pdg_id_ << ") "
                                << "is a new daughter of " << role_ << " (pdgId=" << pdg_id_ << ").";
      if (!utils::contains(part.mothers_, id_))
        part.addMother(*this);
    }
    return *this;
  }

  Particle& Particle::setMomentum(const Momentum& mom, bool offshell) {
    momentum_ = mom;
    if (!offshell)
      momentum_.computeEnergyFromMass(phys_prop_.mass);
    return *this;
  }

  Particle& Particle::setMomentum(double px, double py, double pz, double e) {
    momentum_.setP(px, py, pz).setEnergy(e);
    if (fabs(e - momentum_.energy()) > 1.e-6)  // more than 1 eV difference
      CG_WARNING("Particle") << "Energy difference: " << e - momentum_.energy();
    return *this;
  }

  pdgid_t Particle::pdgId() const { return pdg_id_; }

  Particle& Particle::setPdgId(pdgid_t pdg, short ch) { return setIntegerPdgId(pdg * (ch == 0 ? 1 : ch / abs(ch))); }

  Particle& Particle::setIntegerPdgId(long pdg) {
    pdg_id_ = labs(pdg);
    if (PDG::get().has(pdg_id_)) {
      phys_prop_ = PDG::get()(pdg_id_);
      CG_DEBUG("Particle:setIntegerPdgId") << "Particle PDG id set to " << pdg_id_ << ", "
                                           << "properties set " << phys_prop_ << ".";
    }
    switch (pdg_id_) {
      case 0:
        charge_sign_ = 0.;
        break;
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

  long Particle::integerPdgId() const {
    if (const auto ch = phys_prop_.integerCharge(); ch != 0)
      return static_cast<long>(pdg_id_) * charge_sign_ * (ch / std::abs(ch));
    return static_cast<long>(pdg_id_);
  }

  std::ostream& operator<<(std::ostream& os, const Particle& part) {
    os << std::resetiosflags(std::ios::showbase) << "Particle[" << part.id_ << "]{role=" << part.role_
       << ", status=" << (int)part.status_ << ", "
       << "pdg=" << part.integerPdgId() << ", p4=" << part.momentum_ << " GeV, m=" << part.momentum_.mass() << " GeV, "
       << "p⟂=" << part.momentum_.pt() << " GeV, eta=" << part.momentum_.eta() << ", phi=" << part.momentum_.phi();
    if (part.primary())
      os << ", primary";
    else {
      os << ", " << utils::s("mother", part.mothers_.size());
      if (!part.mothers_.empty()) {
        os << "=";
        std::string delim;
        for (const auto moth : part.mothers_)
          os << delim << moth, delim = ",";
      }
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

  std::ostream& operator<<(std::ostream& os, const Particle::Status& st) {
    switch (st) {
      case Particle::Status::PrimordialIncoming:
        return os << "incoming beam particle";
      case Particle::Status::DebugResonance:
        return os << "intermediate resonance";
      case Particle::Status::Resonance:
        return os << "decayed intermediate resonance";
      case Particle::Status::Fragmented:
        return os << "fragmented outgoing beam";
      case Particle::Status::Propagator:
        return os << "propagator";
      case Particle::Status::Incoming:
        return os << "incoming parton";
      case Particle::Status::Undefined:
        return os << "undefined";
      case Particle::Status::FinalState:
        return os << "final state particle";
      case Particle::Status::Undecayed:
        return os << "particle to be decayed externally";
      case Particle::Status::Unfragmented:
        return os << "particle to be hadronised externally";
    }
    return os;
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

  ParticlesMap::ParticlesMap(const ParticlesMap& oth)
      : std::unordered_map<Particle::Role, Particles, utils::EnumHash<Particle::Role> >() {
    *this = oth;
  }

  ParticlesMap& ParticlesMap::operator=(const ParticlesMap& oth) {
    for (const auto& parts_vs_role : oth)
      for (const auto& part : parts_vs_role.second)
        (*this)[parts_vs_role.first].emplace_back(Particle(part));
    return *this;
  }
}  // namespace cepgen
