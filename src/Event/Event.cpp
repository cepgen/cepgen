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

#include <algorithm>
#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

Event::Event(bool compressed) : compressed_(compressed) {}

Event::Event(const Event& oth) { *this = oth; }

Event& Event::operator=(const Event& oth) {
  clear();
  particles_ = oth.particles_;
  event_content_ = oth.event_content_;
  compressed_ = oth.compressed_;
  metadata = oth.metadata;
  return *this;
}

bool Event::operator==(const Event& oth) const {
  if (compressed_ != oth.compressed_)
    return false;
  if (metadata != oth.metadata)
    CG_WARNING("Event:operator==") << "Comparison of two events with different metadata.";
  if (particles_ != oth.particles_)
    return false;
  return true;
}

Event Event::minimal(size_t num_out_particles) {
  auto evt = Event();
  // add the two incoming beam particles
  auto ib1 = evt.addParticle(Particle::Role::IncomingBeam1);
  ib1.get().setStatus(Particle::Status::PrimordialIncoming);
  auto ib2 = evt.addParticle(Particle::Role::IncomingBeam2);
  ib2.get().setStatus(Particle::Status::PrimordialIncoming);

  // add the two incoming partons
  auto part1 = evt.addParticle(Particle::Role::Parton1);
  part1.get().setStatus(Particle::Status::Incoming);
  part1.get().addMother(ib1);
  auto part2 = evt.addParticle(Particle::Role::Parton2);
  part2.get().setStatus(Particle::Status::Incoming);
  part2.get().addMother(ib2);
  // add the two-parton system
  auto two_parton = evt.addParticle(Particle::Role::Intermediate);
  two_parton.get().setStatus(Particle::Status::Propagator);
  two_parton.get().addMother(part1);
  two_parton.get().addMother(part2);

  // add the two outgoing beam particles
  auto ob1 = evt.addParticle(Particle::Role::OutgoingBeam1);
  ob1.get().setStatus(Particle::Status::FinalState);
  ob1.get().addMother(ib1);
  auto ob2 = evt.addParticle(Particle::Role::OutgoingBeam2);
  ob2.get().setStatus(Particle::Status::FinalState);
  ob2.get().addMother(ib2);

  // finally add the central system
  for (size_t i = 0; i < num_out_particles; ++i) {
    auto cs = evt.addParticle(Particle::Role::CentralSystem);
    cs.get().setStatus(Particle::Status::FinalState);
    cs.get().addMother(two_parton);
  }
  return evt;
}

void Event::clear() {
  particles_.clear();
  metadata.clear();
}

void Event::freeze() {
  if (particles_.count(Particle::Role::CentralSystem) > 0)
    event_content_.cs = particles_[Particle::Role::CentralSystem].size();
  if (particles_.count(Particle::Role::OutgoingBeam1) > 0)
    event_content_.op1 = particles_[Particle::Role::OutgoingBeam1].size();
  if (particles_.count(Particle::Role::OutgoingBeam2) > 0)
    event_content_.op2 = particles_[Particle::Role::OutgoingBeam2].size();
}

void Event::restore() {
  if (particles_.count(Particle::Role::CentralSystem) > 0)
    particles_[Particle::Role::CentralSystem].resize(event_content_.cs);
  if (particles_.count(Particle::Role::OutgoingBeam1) > 0)
    particles_[Particle::Role::OutgoingBeam1].resize(event_content_.op1);
  if (particles_.count(Particle::Role::OutgoingBeam2) > 0)
    particles_[Particle::Role::OutgoingBeam2].resize(event_content_.op2);
}

bool Event::compressed() const { return compressed_; }

Event Event::compress() const {
  if (compressed_)
    return *this;
  Event out(/*compressed=*/true);
  int i = 0;
  //--- add all necessary particles
  for (const auto& role : {Particle::Role::IncomingBeam1,
                           Particle::Role::IncomingBeam2,
                           Particle::Role::OutgoingBeam1,
                           Particle::Role::OutgoingBeam2,
                           Particle::Role::Parton1,
                           Particle::Role::Parton2,
                           Particle::Role::CentralSystem}) {
    if (particles_.count(role) == 0)
      continue;
    for (const auto& old_part : operator()(role)) {
      auto& new_part = out.addParticle(role).get();
      new_part = old_part;  // copy all attributes
      new_part.setId(i++);
      new_part.mothers().clear();
    }
  }
  //--- fix parentage for outgoing beam particles
  if (out[Particle::Role::OutgoingBeam1].size() > 1 || out[Particle::Role::OutgoingBeam2].size() > 1)
    CG_WARNING("Event:compress") << "Event compression not designed for already fragmented beam remnants!\n\t"
                                 << "Particles parentage is not guaranteed to be conserved.";
  if (particles_.count(Particle::Role::OutgoingBeam1) > 0)
    for (auto& part : out[Particle::Role::OutgoingBeam1])
      part.get().addMother(out[Particle::Role::IncomingBeam1][0].get());
  if (particles_.count(Particle::Role::OutgoingBeam2) > 0)
    for (auto& part : out[Particle::Role::OutgoingBeam2])
      part.get().addMother(out[Particle::Role::IncomingBeam2][0].get());
  //--- fix parentage for incoming partons
  for (auto& part : out[Particle::Role::Parton1])
    if (particles_.count(Particle::Role::IncomingBeam1) > 0)
      part.get().addMother(out[Particle::Role::IncomingBeam1][0].get());
  if (particles_.count(Particle::Role::IncomingBeam2) > 0)
    for (auto& part : out[Particle::Role::Parton2])
      part.get().addMother(out[Particle::Role::IncomingBeam2][0].get());
  //--- fix parentage for central system
  if (particles_.count(Particle::Role::Parton1) > 0 && particles_.count(Particle::Role::Parton2) > 0)
    for (auto& part : out[Particle::Role::CentralSystem]) {
      part.get().addMother(out[Particle::Role::Parton1][0]);
      part.get().addMother(out[Particle::Role::Parton2][0]);
    }
  return out;
}

ParticlesRefs Event::operator[](Particle::Role role) {
  ParticlesRefs out;
  //--- retrieve all particles with a given role
  for (auto& part : particles_[role])
    out.emplace_back(std::ref(part));
  return out;
}

const Particles& Event::operator()(Particle::Role role) const {
  if (particles_.count(role) == 0)
    throw CG_FATAL("Event") << "Failed to retrieve a particle with " << role << " role.";
  //--- retrieve all particles with a given role
  return particles_.at(role);
}

ParticlesIds Event::ids(Particle::Role role) const {
  ParticlesIds out;
  //--- retrieve all particles ids with a given role
  if (particles_.count(role) == 0)
    return out;

  for (const auto& part : particles_.at(role))
    out.insert(part.id());

  return out;
}

Particle& Event::oneWithRole(Particle::Role role) {
  //--- retrieve the first particle of a given role
  auto parts_by_role = operator[](role);
  if (parts_by_role.empty())
    throw CG_FATAL("Event") << "No particle retrieved with " << role << " role.";
  if (parts_by_role.size() > 1)
    throw CG_FATAL("Event") << "More than one particle with " << role << " role: " << parts_by_role.size()
                            << " particles.";
  return parts_by_role.front().get();
}

const Particle& Event::oneWithRole(Particle::Role role) const {
  //--- retrieve the first particle of a given role
  const Particles& parts_by_role = operator()(role);
  if (parts_by_role.empty())
    throw CG_FATAL("Event") << "No particle retrieved with " << role << " role.";
  if (parts_by_role.size() > 1)
    throw CG_FATAL("Event") << "More than one particle with " << role << " role: " << parts_by_role.size()
                            << " particles";
  return *parts_by_role.begin();
}

Particle& Event::operator[](int id) {
  for (auto& role_part : particles_)
    for (auto& part : role_part.second)
      if (part.id() == id)
        return part;

  throw CG_FATAL("Event") << "Failed to retrieve the particle with id=" << id << ".";
}

const Particle& Event::operator()(int id) const {
  for (const auto& role_part : particles_) {
    if (auto it = std::find_if(
            role_part.second.begin(), role_part.second.end(), [&id](const auto& part) { return part.id() == id; });
        it != role_part.second.end())
      return *it;
  }
  throw CG_FATAL("Event") << "Failed to retrieve the particle with id=" << id << ".";
}

ParticlesRefs Event::operator[](const ParticlesIds& ids) {
  ParticlesRefs out;
  std::transform(ids.begin(), ids.end(), std::back_inserter(out), [this](const auto& id) {
    return std::ref(this->operator[](id));
  });
  return out;
}

Particles Event::operator()(const ParticlesIds& ids) const {
  Particles out;
  std::transform(
      ids.begin(), ids.end(), std::back_inserter(out), [this](const auto& id) { return this->operator()(id); });
  return out;
}

Particles Event::mothers(const Particle& part) const { return operator()(part.mothers()); }

ParticlesRefs Event::mothers(const Particle& part) { return operator[](part.mothers()); }

void Event::clearMothers(Particle& part) {
  for (const auto& mid : part.mothers())
    operator[](mid).daughters().erase(part.id());
  part.mothers().clear();
}

Particles Event::daughters(const Particle& part) const { return operator()(part.daughters()); }

ParticlesRefs Event::daughters(const Particle& part) { return operator[](part.daughters()); }

Particles Event::stableDaughters(const Particle& part, bool recursive) const {
  Particles parts;
  for (const auto& daughter : operator()(part.daughters())) {
    if (daughter.status() == Particle::Status::FinalState)
      parts.emplace_back(daughter);
    else if (recursive) {
      const auto stable_daughters = stableDaughters(daughter, recursive);
      parts.insert(parts.end(),
                   std::make_move_iterator(stable_daughters.begin()),
                   std::make_move_iterator(stable_daughters.end()));
    }
  }
  std::sort(parts.begin(), parts.end());
  parts.erase(std::unique(parts.begin(), parts.end()), parts.end());
  return parts;
}

ParticlesRefs Event::stableDaughters(const Particle& part, bool recursive) {
  ParticlesRefs parts;
  for (const auto& daughter : operator[](part.daughters())) {
    if (daughter.get().status() == Particle::Status::FinalState)
      parts.emplace_back(daughter);
    else if (recursive) {
      const auto stable_daughters = stableDaughters(daughter, recursive);
      parts.insert(parts.end(),
                   std::make_move_iterator(stable_daughters.begin()),
                   std::make_move_iterator(stable_daughters.end()));
    }
  }
  std::sort(
      parts.begin(), parts.end(), [](const auto& lhs, const auto& rhs) { return lhs.get().id() < rhs.get().id(); });
  parts.erase(
      std::unique(parts.begin(), parts.end(), [](const auto& lhs, const auto& rhs) { return lhs.get() == rhs.get(); }),
      parts.end());
  return parts;
}

void Event::clearDaughters(Particle& part) {
  for (const auto& did : part.daughters())
    operator[](did).mothers().erase(part.id());
  part.daughters().clear();
}

ParticleRoles Event::roles() const {
  ParticleRoles out;
  std::transform(
      particles_.begin(), particles_.end(), std::back_inserter(out), [](const auto& pr) { return pr.first; });
  return out;
}

void Event::updateRoles() {
  for (auto& parts_vs_role : particles_)  // 1st loop to copy particles to correct role container
    for (const auto& part : parts_vs_role.second)
      if (part.role() != parts_vs_role.first)
        particles_[part.role()].emplace_back(part);
  for (auto& parts_vs_role : particles_)  // 2nd loop to remove wrongly-assigned particles
    parts_vs_role.second.erase(
        std::remove_if(parts_vs_role.second.begin(),
                       parts_vs_role.second.end(),
                       [&parts_vs_role](const auto& part) { return part.role() != parts_vs_role.first; }),
        parts_vs_role.second.end());
}

ParticleRef Event::addParticle(Particle& part, bool replace) {
  CG_DEBUG_LOOP("Event") << "Particle with PDGid = " << part.integerPdgId() << " has role " << part.role();
  if (static_cast<int>(part.role()) <= 0)
    throw CG_FATAL("Event") << "Trying to add a particle with role=" << static_cast<int>(part.role()) << ".";

  auto& part_with_same_role = particles_[part.role()];  // list of particles with the same role
  if (part.id() < 0)
    part.setId(part_with_same_role.empty() || !replace
                   ? size()                         // set the id if previously invalid/non-existent
                   : part_with_same_role[0].id());  // set the previous id if replacing a particle
  if (replace)
    part_with_same_role = {part};
  else
    part_with_same_role.emplace_back(part);
  return std::ref(part_with_same_role.back());
}

ParticleRef Event::addParticle(Particle::Role role, bool replace) {
  Particle np(role, PDG::invalid);
  return addParticle(np, replace);
}

size_t Event::size() const {
  return std::accumulate(particles_.begin(), particles_.end(), 0, [](size_t size, const auto& role_part) {
    return size + role_part.second.size();
  });
}

bool Event::empty() const { return particles_.empty(); }

Particles Event::particles() const {
  Particles out;
  for (const auto& role_part : particles_)
    out.insert(out.end(), role_part.second.begin(), role_part.second.end());

  std::sort(out.begin(), out.end());
  return out;
}

Particles Event::stableParticles() const {
  Particles out;
  for (const auto& role_part : particles_)
    std::copy_if(role_part.second.begin(), role_part.second.end(), std::back_inserter(out), [](const auto& part) {
      return static_cast<short>(part.status()) > 0;
    });

  std::sort(out.begin(), out.end());
  return out;
}

Particles Event::stableParticlesWithRole(Particle::Role role) const {
  Particles out;
  const auto& parts_role = particles_.at(role);
  std::copy_if(parts_role.begin(), parts_role.end(), std::back_inserter(out), [](const auto& part) {
    return static_cast<short>(part.status()) > 0;
  });
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

Momentum Event::missingMomentum() const {
  Momentum missing_momentum;
  for (const auto& cp : operator()(Particle::Role::CentralSystem))
    if (cp.status() == Particle::Status::FinalState) {
      const auto pdg = cp.integerPdgId();
      if (pdg == 12 || pdg == 14 || pdg == 16)  // neutrinos
        missing_momentum += cp.momentum();
      if (pdg == 1000022 || pdg == 1000023 || pdg == 1000025 || pdg == 1000035)  // neutralinos
        missing_momentum += cp.momentum();
    }
  return missing_momentum;
}

void Event::checkKinematics() const {
  // check the kinematics through parentage
  for (const auto& part : particles()) {
    ParticlesIds daughters = part.daughters();
    if (daughters.empty())
      continue;
    Momentum total_momentum;
    for (const auto& daughter : daughters) {
      const auto& d = operator()(daughter);
      const auto mothers = d.mothers();
      total_momentum += d.momentum();
      if (mothers.size() < 2)
        continue;
      for (const auto& moth : mothers)
        if (moth != part.id())
          total_momentum -= operator()(moth).momentum();
    }
    const double mass_diff = (total_momentum - part.momentum()).mass();
    if (fabs(mass_diff) > MIN_PRECISION) {
      dump();
      throw CG_FATAL("Event") << "Error in momentum balance for particle " << part.id() << ": mdiff = " << mass_diff
                              << ".";
    }
  }
}

void Event::dump() const { CG_INFO("Event") << *this; }

double Event::cmEnergy() const {
  return (oneWithRole(Particle::Role::IncomingBeam1).momentum() + oneWithRole(Particle::Role::IncomingBeam2).momentum())
      .mass();
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& out, const Event& ev) {
    const auto parts = ev.particles();
    std::ostringstream os;

    Momentum p_total;
    for (const auto& part : parts) {
      const ParticlesIds mothers = part.mothers();
      {
        std::ostringstream oss_pdg;
        if (part.pdgId() == PDG::invalid && !mothers.empty()) {
          //--- if particles compound
          std::string delim;
          for (size_t i = 0; i < mothers.size(); ++i)
            try {
              oss_pdg << delim << static_cast<PDG::Id>(ev(*std::next(mothers.begin(), i)).pdgId()), delim = "/";
            } catch (const Exception&) {
              oss_pdg << delim << ev(*std::next(mothers.begin(), i)).pdgId(), delim = "/";
            }
          os << utils::format("\n %2d\t\t   %-7s", part.id(), oss_pdg.str().c_str());
        } else {
          //--- if single particle/HI
          if (HeavyIon::isHI(part.pdgId()))
            oss_pdg << HeavyIon::fromPdgId(part.pdgId());
          else
            try {
              oss_pdg << static_cast<PDG::Id>(part.pdgId());
            } catch (const Exception&) {
              oss_pdg << "?";
            }
          os << utils::format("\n %2d\t%-+10d %-7s", part.id(), part.integerPdgId(), oss_pdg.str().c_str());
        }
      }
      os << "\t";
      if (part.charge() != static_cast<int>(part.charge())) {
        if (part.charge() * 2 == static_cast<int>(part.charge() * 2))
          os << utils::format("%-d/2", static_cast<int>(part.charge() * 2));
        else if (part.charge() * 3 == static_cast<int>(part.charge() * 3))
          os << utils::format("%-d/3", static_cast<int>(part.charge() * 3));
        else
          os << utils::format("%-.2f", part.charge());
      } else
        os << utils::format("%-g", part.charge());
      os << "\t";
      {
        std::ostringstream oss;
        oss << part.role();
        os << utils::format("%-8s %6d\t", oss.str().c_str(), part.status());
      }
      if (!mothers.empty()) {
        std::ostringstream oss;
        unsigned short i = 0;
        for (const auto& moth : mothers) {
          oss << (i > 0 ? "+" : "") << moth;
          ++i;
        }
        os << utils::format("%6s ", oss.str().c_str());
      } else
        os << "       ";
      const auto& mom = part.momentum();
      os << utils::format(
          "% 9.6e % 9.6e % 9.6e % 9.6e % 12.5f", mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass());

      // discard non-primary, decayed particles
      if (part.status() >= Particle::Status::Undefined) {
        const int sign = (part.status() == Particle::Status::Undefined) ? -1 : +1;
        p_total += sign * mom;
      }
    }
    //--- set a threshold to the computation precision
    p_total.truncate();
    return out << utils::format(
               "Event content:\n"
               " Id\tPDG id\t   Name\t\tCharge\tRole\t Status\tMother\tpx            py            pz            E   "
               "  "
               " \t M         \n"
               " --\t------\t   ----\t\t------\t----\t ------\t------\t----GeV/c---  ----GeV/c---  ----GeV/c---  "
               "----GeV/c---\t --GeV/c²--"
               "%s\n"
               " ----------------------------------------------------------------------------------------------------"
               "--"
               "----------------------------\n"
               "\t\t\t\t\t\t\tBalance% 9.6e % 9.6e % 9.6e % 9.6e",
               os.str().c_str(),
               p_total.px(),
               p_total.py(),
               p_total.pz(),
               p_total.energy());
  }
}  // namespace cepgen

Event::EventMetadata::EventMetadata()
    : std::unordered_map<std::string, float>{{"time:generation", -1.f},
                                             {"time:total", -1.f},
                                             {"weight", 1.f},
                                             {"alphaEM", static_cast<float>(constants::ALPHA_EM)},
                                             {"alphaS", static_cast<float>(constants::ALPHA_QCD)}} {}
