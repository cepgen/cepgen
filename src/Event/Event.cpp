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
using namespace std::string_literals;

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
  auto& ib1 = evt.addParticle(Particle::Role::IncomingBeam1).get();
  ib1.setStatus(Particle::Status::PrimordialIncoming);
  auto& ib2 = evt.addParticle(Particle::Role::IncomingBeam2).get();
  ib2.setStatus(Particle::Status::PrimordialIncoming);

  // add the two incoming partons
  auto& part1 = evt.addParticle(Particle::Role::Parton1).get();
  part1.setStatus(Particle::Status::Incoming);
  part1.addMother(ib1);
  auto& part2 = evt.addParticle(Particle::Role::Parton2).get();
  part2.setStatus(Particle::Status::Incoming);
  part2.addMother(ib2);
  // add the two-parton system
  auto& two_parton = evt.addParticle(Particle::Role::Intermediate).get();
  two_parton.setStatus(Particle::Status::Propagator);
  two_parton.addMother(part1);
  two_parton.addMother(part2);

  // add the two outgoing beam particles
  evt.addParticle(Particle::Role::OutgoingBeam1).get().setStatus(Particle::Status::FinalState).addMother(ib1);
  evt.addParticle(Particle::Role::OutgoingBeam2).get().setStatus(Particle::Status::FinalState).addMother(ib2);

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
  for (const auto& role : {Particle::Role::IncomingBeam1,
                           Particle::Role::IncomingBeam2,
                           Particle::Role::OutgoingBeam1,
                           Particle::Role::OutgoingBeam2,
                           Particle::Role::Parton1,
                           Particle::Role::Parton2,
                           Particle::Role::CentralSystem}) {  // add all necessary particles
    if (particles_.count(role) == 0)
      continue;
    for (const auto& old_particle : operator()(role)) {
      auto& new_particle = out.addParticle(role).get();
      new_particle = old_particle;  // copy all attributes
      new_particle.setId(i++);
      new_particle.mothers().clear();
    }
  }
  // fix parentage for outgoing beam particles
  if (out[Particle::Role::OutgoingBeam1].size() > 1 || out[Particle::Role::OutgoingBeam2].size() > 1)
    CG_WARNING("Event:compress") << "Event compression not designed for already fragmented beam remnants!\n\t"
                                 << "Particles parentage is not guaranteed to be conserved.";
  if (particles_.count(Particle::Role::OutgoingBeam1) > 0)
    for (const auto& part : out[Particle::Role::OutgoingBeam1])
      part.get().addMother(out[Particle::Role::IncomingBeam1][0].get());
  if (particles_.count(Particle::Role::OutgoingBeam2) > 0)
    for (const auto& part : out[Particle::Role::OutgoingBeam2])
      part.get().addMother(out[Particle::Role::IncomingBeam2][0].get());

  // fix parentage for incoming partons
  for (const auto& part : out[Particle::Role::Parton1])
    if (particles_.count(Particle::Role::IncomingBeam1) > 0)
      part.get().addMother(out[Particle::Role::IncomingBeam1][0].get());
  if (particles_.count(Particle::Role::IncomingBeam2) > 0)
    for (const auto& part : out[Particle::Role::Parton2])
      part.get().addMother(out[Particle::Role::IncomingBeam2][0].get());

  // fix parentage for central system
  if (particles_.count(Particle::Role::Parton1) > 0 && particles_.count(Particle::Role::Parton2) > 0)
    for (const auto& part : out[Particle::Role::CentralSystem]) {
      part.get().addMother(out[Particle::Role::Parton1][0]);
      part.get().addMother(out[Particle::Role::Parton2][0]);
    }
  return out;
}

ParticlesRefs Event::operator[](Particle::Role role) {
  ParticlesRefs out;
  for (auto& part : particles_[role])
    out.emplace_back(std::ref(part));
  return out;
}

const Particles& Event::operator()(Particle::Role role) const {
  if (particles_.count(role) == 0)
    throw CG_FATAL("Event") << "Failed to retrieve a particle with " << role << " role.";
  return particles_.at(role);
}

ParticlesIds Event::ids(Particle::Role role) const {
  ParticlesIds out;
  if (particles_.count(role) == 0)
    return out;
  for (const auto& part : particles_.at(role))
    out.insert(part.id());
  return out;
}

Particle& Event::oneWithRole(Particle::Role role) {
  const auto parts_by_role = operator[](role);
  if (parts_by_role.empty())
    throw CG_FATAL("Event") << "No particle retrieved with " << role << " role.";
  if (parts_by_role.size() > 1)
    throw CG_FATAL("Event") << "More than one particle with " << role << " role: " << parts_by_role.size()
                            << " particles.";
  return parts_by_role.front().get();
}

const Particle& Event::oneWithRole(Particle::Role role) const {
  const Particles& parts_by_role = operator()(role);
  if (parts_by_role.empty())
    throw CG_FATAL("Event") << "No particle retrieved with " << role << " role.";
  if (parts_by_role.size() > 1)
    throw CG_FATAL("Event") << "More than one particle with " << role << " role: " << parts_by_role.size()
                            << " particles";
  return *parts_by_role.begin();
}

Particle& Event::operator[](int id) {
  for (auto& [role, particles] : particles_)
    for (auto& particle : particles)
      if (particle.id() == id)
        return particle;
  throw CG_FATAL("Event") << "Failed to retrieve the particle with id=" << id << ".";
}

const Particle& Event::operator()(int id) const {
  for (const auto& [role, particles] : particles_) {
    if (auto it = std::find_if(
            particles.begin(), particles.end(), [&id](const auto& particle) { return particle.id() == id; });
        it != particles.end())
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

Particles Event::mothers(const Particle& particle) const { return operator()(particle.mothers()); }

ParticlesRefs Event::mothers(const Particle& particle) { return operator[](particle.mothers()); }

void Event::clearMothers(Particle& particle) {
  for (const auto& mother_id : particle.mothers())
    operator[](mother_id).children().erase(particle.id());
  particle.mothers().clear();
}

Particles Event::children(const Particle& particle) const { return operator()(particle.children()); }

ParticlesRefs Event::children(const Particle& particle) { return operator[](particle.children()); }

Particles Event::stableChildren(const Particle& particle, bool recursive) const {
  Particles children;
  for (const auto& child : operator()(particle.children())) {
    if (child.status() == Particle::Status::FinalState)
      children.emplace_back(child);
    else if (recursive) {
      const auto stable_children = stableChildren(child, recursive);
      children.insert(children.end(),
                      std::make_move_iterator(stable_children.begin()),
                      std::make_move_iterator(stable_children.end()));
    }
  }
  std::sort(children.begin(), children.end());
  children.erase(std::unique(children.begin(), children.end()), children.end());
  return children;
}

ParticlesRefs Event::stableChildren(const Particle& particle, bool recursive) {
  ParticlesRefs children;
  for (const auto& child : operator[](particle.children())) {
    if (child.get().status() == Particle::Status::FinalState)
      children.emplace_back(child);
    else if (recursive) {
      const auto stable_children = stableChildren(child, recursive);
      children.insert(children.end(),
                      std::make_move_iterator(stable_children.begin()),
                      std::make_move_iterator(stable_children.end()));
    }
  }
  std::sort(children.begin(), children.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.get().id() < rhs.get().id();
  });
  children.erase(
      std::unique(
          children.begin(), children.end(), [](const auto& lhs, const auto& rhs) { return lhs.get() == rhs.get(); }),
      children.end());
  return children;
}

void Event::clearChildren(Particle& particle) {
  for (const auto& child_id : particle.children())
    operator[](child_id).mothers().erase(particle.id());
  particle.children().clear();
}

ParticleRoles Event::roles() const {
  ParticleRoles particle_roles;
  std::transform(particles_.begin(), particles_.end(), std::back_inserter(particle_roles), [](const auto& pr) {
    return pr.first;
  });
  return particle_roles;
}

void Event::updateRoles() {
  for (auto& [role, particles] : particles_)  // 1st loop to copy particles to correct role container
    for (const auto& particle : particles)
      if (particle.role() != role)
        particles_[particle.role()].emplace_back(particle);
  for (auto& parts_vs_role : particles_)  // 2nd loop to remove wrongly-assigned particles
    parts_vs_role.second.erase(
        std::remove_if(parts_vs_role.second.begin(),
                       parts_vs_role.second.end(),
                       [&parts_vs_role](const auto& part) { return part.role() != parts_vs_role.first; }),
        parts_vs_role.second.end());
}

ParticleRef Event::addParticle(Particle& particle, bool replace) {
  CG_DEBUG_LOOP("Event") << "Particle with PDGid = " << particle.integerPdgId() << " has role " << particle.role();
  if (static_cast<int>(particle.role()) <= 0)
    throw CG_FATAL("Event") << "Trying to add a particle with role=" << static_cast<int>(particle.role()) << ".";

  auto& part_with_same_role = particles_[particle.role()];  // list of particles with the same role
  if (particle.id() < 0)
    particle.setId(part_with_same_role.empty() || !replace
                       ? size()                         // set the id if previously invalid/non-existent
                       : part_with_same_role[0].id());  // set the previous id if replacing a particle
  if (replace)
    part_with_same_role = {particle};
  else
    part_with_same_role.emplace_back(particle);
  return std::ref(part_with_same_role.back());
}

ParticleRef Event::addParticle(Particle::Role role, bool replace) {
  Particle new_particle(role, PDG::invalid);
  return addParticle(new_particle, replace);
}

size_t Event::size() const {
  return std::accumulate(particles_.begin(), particles_.end(), 0, [](size_t size, const auto& role_part) {
    return size + role_part.second.size();
  });
}

bool Event::empty() const { return particles_.empty(); }

Particles Event::particles() const {
  Particles all_particles;
  for (const auto& [role, particles] : particles_)
    all_particles.insert(all_particles.end(), particles.begin(), particles.end());
  std::sort(all_particles.begin(), all_particles.end());
  return all_particles;
}

Particles Event::stableParticles() const {
  Particles stable_particles;
  for (const auto& [role, particles] : particles_)
    std::copy_if(particles.begin(), particles.end(), std::back_inserter(stable_particles), [](const auto& particle) {
      return static_cast<short>(particle.status()) > 0;
    });
  std::sort(stable_particles.begin(), stable_particles.end());
  return stable_particles;
}

Particles Event::stableParticlesWithRole(Particle::Role role) const {
  Particles particles_with_role;
  if (particles_.count(role) == 0)
    throw CG_FATAL("Event:stableParticlesWithRole") << "Failed to retrieve role '" << role << "' from event.";
  const auto& parts_role = particles_.at(role);
  std::copy_if(parts_role.begin(), parts_role.end(), std::back_inserter(particles_with_role), [](const auto& part) {
    return static_cast<short>(part.status()) > 0;
  });
  std::sort(particles_with_role.begin(), particles_with_role.end());
  particles_with_role.erase(std::unique(particles_with_role.begin(), particles_with_role.end()),
                            particles_with_role.end());
  return particles_with_role;
}

Momentum Event::missingMomentum() const {
  Momentum missing_momentum;
  for (const auto& central_particle : operator()(Particle::Role::CentralSystem))
    if (central_particle.status() == Particle::Status::FinalState) {
      const auto pdg = central_particle.integerPdgId();
      if (pdg == 12 || pdg == 14 || pdg == 16)  // neutrinos
        missing_momentum += central_particle.momentum();
      if (pdg == 1000022 || pdg == 1000023 || pdg == 1000025 || pdg == 1000035)  // neutralinos
        missing_momentum += central_particle.momentum();
    }
  return missing_momentum;
}

void Event::checkKinematics() const {
  for (const auto& particle : particles()) {  // check the kinematics through parentage
    const auto children = particle.children();
    if (children.empty())
      continue;
    Momentum total_momentum;
    for (const auto& child : children) {
      const auto& d = operator()(child);
      total_momentum += d.momentum();
      if (const auto mothers = d.mothers(); mothers.size() >= 2)
        for (const auto& moth : mothers)
          if (moth != particle.id())
            total_momentum -= operator()(moth).momentum();
    }
    if (const auto mass_diff = (total_momentum - particle.momentum()).mass(); std::fabs(mass_diff) > MIN_PRECISION)
      throw CG_FATAL("Event") << "Error in momentum balance for particle " << particle.id()
                              << ": mass difference=" << mass_diff << " GeV. Event dump:\n"
                              << *this;
  }
}

void Event::dump() const { CG_INFO("Event") << *this; }

double Event::cmEnergy() const {
  return (oneWithRole(Particle::Role::IncomingBeam1).momentum() + oneWithRole(Particle::Role::IncomingBeam2).momentum())
      .mass();
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& out, const Event& event) {
    std::ostringstream os;
    Momentum p_total;
    for (const auto& particle : event.particles()) {
      const ParticlesIds mothers = particle.mothers();
      {
        std::ostringstream oss_pdg;
        if (particle.pdgId() == PDG::invalid && !mothers.empty()) {  // particles compound
          std::string delim;
          for (size_t i = 0; i < mothers.size(); ++i)
            try {
              oss_pdg << delim << static_cast<PDG::Id>(event(*std::next(mothers.begin(), i)).pdgId()), delim = "/";
            } catch (const Exception&) {
              oss_pdg << delim << event(*std::next(mothers.begin(), i)).pdgId(), delim = "/";
            }
          os << utils::format("\n %2d\t\t   %-7s", particle.id(), oss_pdg.str().c_str());
        } else {
          if (HeavyIon::isHI(particle.pdgId()))
            oss_pdg << HeavyIon::fromPdgId(particle.pdgId());
          else  // single particle
            try {
              oss_pdg << static_cast<PDG::Id>(particle.pdgId());
            } catch (const Exception&) {
              oss_pdg << "?";
            }
          os << utils::format("\n %2d\t%-+10d %-7s", particle.id(), particle.integerPdgId(), oss_pdg.str().c_str());
        }
      }
      os << "\t";
      if (particle.charge() == static_cast<int>(particle.charge()))  // integer charge
        os << utils::format("%-g", particle.charge());
      else if (particle.charge() * 2 == static_cast<int>(particle.charge() * 2))  // n/2 charge
        os << utils::format("%-d/2", static_cast<int>(particle.charge() * 2));
      else if (particle.charge() * 3 == static_cast<int>(particle.charge() * 3))  // n/3 charge
        os << utils::format("%-d/3", static_cast<int>(particle.charge() * 3));
      else
        os << utils::format("%-.2f", particle.charge());
      os << "\t";
      {
        std::ostringstream oss;
        oss << particle.role();
        os << utils::format("%-8s %6d\t", oss.str().c_str(), particle.status());
      }
      const auto& momentum = particle.momentum();
      os << utils::format("%6s % 9.6e % 9.6e % 9.6e % 9.6e % 12.5f",
                          (!mothers.empty() ? utils::repr(mothers, "+"s) : ""s).data(),
                          momentum.px(),
                          momentum.py(),
                          momentum.pz(),
                          momentum.energy(),
                          momentum.mass());
      if (particle.status() >= Particle::Status::Undefined)  // discard non-primary, decayed particles
        p_total += (particle.status() == Particle::Status::Undefined ? -1 : +1) * momentum;
    }
    p_total.truncate();  // set a threshold to the computation precision
    return out << utils::format(
               "Event content:\n"
               " Id\tPDG id\t   Name\t\tCharge\tRole\t Status\tMother\tpx            py            pz            E     "
               " \t M         \n"
               " --\t------\t   ----\t\t------\t----\t ------\t------\t----GeV/c---  ----GeV/c---  ----GeV/c---  "
               "----GeV/c---\t --GeV/c²--%s\n"
               " %s\n"
               "\t\t\t\t\t\t\tBalance% 9.6e % 9.6e % 9.6e % 9.6e",
               os.str().c_str(),
               std::string(130, '-').data(),
               p_total.px(),
               p_total.py(),
               p_total.pz(),
               p_total.energy());
  }
}  // namespace cepgen

Event::EventMetadata::EventMetadata()
    : std::unordered_map<std::string, float>{{"time:generation"s, -1.f},
                                             {"time:total"s, -1.f},
                                             {"weight"s, 1.f},
                                             {"alphaEM"s, constants::ALPHA_EM},
                                             {"alphaS"s, constants::ALPHA_QCD}} {}
