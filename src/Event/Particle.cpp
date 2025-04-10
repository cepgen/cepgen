/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include <cmath>
#include <iomanip>

#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

Particle::Particle(Role role, pdgid_t pdgId, Status st) : role_(role), status_(static_cast<int>(st)), pdg_id_(pdgId) {}

bool Particle::operator<(const Particle& rhs) const { return id_ >= 0 && rhs.id_ > 0 && id_ < rhs.id_; }

bool Particle::operator==(const Particle& oth) const {
  return id_ == oth.id_ && pdg_id_ == oth.pdg_id_ && antiparticle_ == oth.antiparticle_ && helicity_ == oth.helicity_ &&
         status_ == oth.status_ && momentum_ == oth.momentum_;
}

bool Particle::valid() const {
  if (pdg_id_ == PDG::invalid)
    return false;
  if (momentum_.p() == 0. && momentum_.mass() == 0.)
    return false;
  return true;
}

float Particle::charge() const { return (antiparticle_ ? -1 : +1) * PDG::get()(pdg_id_).integerCharge() / 3.; }

Particle& Particle::addMother(Particle& mother_particle) {
  if (const auto& [it, inserted] = mothers_.insert(mother_particle.id_); inserted) {
    CG_DEBUG_LOOP("Particle:addMother") << "Particle (id=" << *it << ", role=" << mother_particle.role_
                                        << ", pdgId=" << mother_particle.pdg_id_ << ") is a new mother of (id=" << id_
                                        << ", role=" << role_ << ", pdgId=" << pdg_id_ << ").";
    if (mother_particle.children_.empty() ||
        !utils::contains(mother_particle.children_, id_))  // not yet in particle's children
      mother_particle.addChild(*this);
    else if (mother_particle.status_ > 0)
      mother_particle.status_ = static_cast<int>(Status::Propagator);
  }
  return *this;
}

Particle& Particle::addChild(Particle& child_particle) {
  if (const auto& [it, inserted] = children_.insert(child_particle.id_); inserted) {
    CG_DEBUG_LOOP("Particle:addChild") << "Particle (id=" << *it << ", role=" << child_particle.role_
                                       << ", pdgId=" << child_particle.pdg_id_ << ") is a new child of (id=" << id_
                                       << ", role=" << role_ << ", pdgId=" << pdg_id_ << ").";
    if (child_particle.mothers_.empty() ||
        !utils::contains(child_particle.mothers_, id_))  // not yet in particle's mothers
      child_particle.addMother(*this);
    if (status_ > 0)
      status_ = static_cast<int>(Status::Propagator);
  }
  CG_DEBUG_LOOP("Particle:addChild").log([this](auto& dbg) {
    dbg << "Particle (id=" << id_ << ", role=" << role_ << ", pdgId=" << static_cast<int>(pdg_id_) << ")"
        << " has now " << children_.size() << "child(ren):";
    for (const auto& child : children_)
      dbg << utils::format("\n\t * id=%d", child);
  });
  return *this;
}

Particle& Particle::setMomentum(const Momentum& mom, bool off_shell) {
  momentum_ = mom;
  if (!off_shell)
    momentum_.computeEnergyFromMass(PDG::get().mass(pdg_id_));
  return *this;
}

Particle& Particle::setMomentum(double px, double py, double pz, double energy) {
  momentum_.setP(px, py, pz).setEnergy(energy);
  if (std::fabs(energy - momentum_.energy()) > 1.e-6)  // more than 1 keV difference
    CG_WARNING("Particle") << "Energy difference: " << energy - momentum_.energy();
  return *this;
}

pdgid_t Particle::pdgId() const { return pdg_id_; }

Particle& Particle::setPdgId(pdgid_t pdg, short ch) { return setIntegerPdgId(pdg * (ch == 0 ? 1 : ch / abs(ch))); }

Particle& Particle::setIntegerPdgId(long pdg) {
  if (pdg_id_ = std::labs(pdg); PDG::get().has(pdg_id_))
    CG_DEBUG("Particle:setIntegerPdgId") << "Particle PDG id set to " << pdg_id_ << ".";
  antiparticle_ = pdg < 0;
  return *this;
}

long Particle::integerPdgId() const { return static_cast<long>(pdg_id_) * (antiparticle_ ? -1 : +1); }

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const Particle& part) {
    os << std::resetiosflags(std::ios::showbase) << "Particle[" << part.id_ << "]{role=" << part.role_
       << ", status=" << part.status_ << ", "
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
    if (!part.children_.empty()) {
      os << ", " << utils::s("child", part.children_.size()) << "=";
      std::string delim;
      for (const auto& child : part.children_)
        os << delim << child, delim = ",";
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
      case Particle::Role::UnknownRole:
        return os << "unknown";
      case Particle::Role::IncomingBeam1:
        return os << "i.beam 1";
      case Particle::Role::IncomingBeam2:
        return os << "i.beam 2";
      case Particle::Role::OutgoingBeam1:
        return os << "o.beam 1";
      case Particle::Role::OutgoingBeam2:
        return os << "o.beam 2";
      case Particle::Role::Parton1:
        return os << "parton 1";
      case Particle::Role::Parton2:
        return os << "parton 2";
      case Particle::Role::Intermediate:
        return os << "hard pr.";
      case Particle::Role::CentralSystem:
        return os << "central";
    }
    return os;
  }
}  // namespace cepgen

ParticlesMap::ParticlesMap(const ParticlesMap& oth)
    : std::unordered_map<Particle::Role, Particles, utils::EnumHash<Particle::Role> >() {
  *this = oth;
}

ParticlesMap& ParticlesMap::operator=(const ParticlesMap& oth) {
  for (const auto& [role, particles] : oth)
    for (const auto& part : particles)
      (*this)[role].emplace_back(Particle(part));
  return *this;
}
