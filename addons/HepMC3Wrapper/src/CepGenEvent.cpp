/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2026  Laurent Forthomme
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

#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>

#include <functional>
#include <list>
#include <numeric>
#include <unordered_set>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Collections.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace HepMC3;

CepGenEvent::CepGenEvent(const cepgen::Event& event) : GenEvent(Units::GEV, Units::MM) {
  add_attribute("AlphaQCD", make_shared<DoubleAttribute>(event.metadata("alphaS")));
  add_attribute("AlphaEM", make_shared<DoubleAttribute>(event.metadata("alphaEM")));

  weights().push_back(event.metadata("weight"));  // unit weight for unweighted events

  const FourVector origin(0., 0., 0., 0.);
  int central_system_id = 0;
  auto vertex_beam1 = make_shared<GenVertex>(origin), vertex_beam2 = make_shared<GenVertex>(origin),
       vertex_central_system = make_shared<GenVertex>(origin);
  std::vector<GenVertexPtr> decay_vertices;  // added after the three primordial vertices, expected first by merge()
  for (const auto& cepgen_particle : event.particles()) {  // filling the particles content
    const auto& cepgen_momentum = cepgen_particle.momentum();
    const auto momentum =
        FourVector(cepgen_momentum.px(), cepgen_momentum.py(), cepgen_momentum.pz(), cepgen_momentum.energy());
    const auto hepmc_particle =
        make_shared<GenParticle>(momentum, cepgen_particle.integerPdgId(), static_cast<int>(cepgen_particle.status()));
    hepmc_particle->set_generated_mass(cepgen::PDG::get().mass(cepgen_particle.pdgId()));
    cepgen_id_vs_hepmc_particle_[cepgen_particle.id()] = hepmc_particle;

    switch (cepgen_particle.role()) {
      case cepgen::Particle::Role::IncomingBeam1:
        vertex_beam1->add_particle_in(hepmc_particle);
        break;
      case cepgen::Particle::Role::IncomingBeam2:
        vertex_beam2->add_particle_in(hepmc_particle);
        break;
      case cepgen::Particle::Role::OutgoingBeam1:
        vertex_beam1->add_particle_out(hepmc_particle);
        break;
      case cepgen::Particle::Role::OutgoingBeam2:
        vertex_beam2->add_particle_out(hepmc_particle);
        break;
      case cepgen::Particle::Role::Parton1:
        vertex_beam1->add_particle_out(hepmc_particle);
        vertex_central_system->add_particle_in(hepmc_particle);
        break;
      case cepgen::Particle::Role::Parton2:
        vertex_beam2->add_particle_out(hepmc_particle);
        vertex_central_system->add_particle_in(hepmc_particle);
        break;
      case cepgen::Particle::Role::Intermediate:  // skip the two-parton system and propagate the parentage
        central_system_id = cepgen_particle.id();
        continue;
      case cepgen::Particle::Role::CentralSystem:
      default: {
        const auto& mothers = cepgen_particle.mothers();
        if (mothers.empty())
          continue;  // skip disconnected lines
        // check if particle is connected to the two-parton system
        if (const auto m1 = *mothers.begin(), m2 = mothers.size() > 1 ? *mothers.rbegin() : -1;  // get mother(s) id(s)
            m1 == central_system_id ||
            (m2 >= 0 && (m1 < central_system_id && central_system_id <= m2)))  // also supports range
          vertex_central_system->add_particle_out(hepmc_particle);
        else if (cepgen_id_vs_hepmc_particle_.count(m1) !=
                 0) {  // if part of the decay chain of central system, find parents
          auto production_vertex = cepgen_id_vs_hepmc_particle_.at(m1)->end_vertex();
          std::list ids{m1};  // list of mother particles
          if (cepgen_id_vs_hepmc_particle_.count(m2) != 0 && m2 > m1) {
            ids.resize(m2 - m1 + 1);
            std::iota(ids.begin(), ids.end(), m1);
          }
          if (!production_vertex) {
            production_vertex = make_shared<GenVertex>();
            for (const auto& id : ids)
              production_vertex->add_particle_in(cepgen_id_vs_hepmc_particle_.at(id));
            decay_vertices.emplace_back(production_vertex);
          }
          production_vertex->add_particle_out(hepmc_particle);
        } else
          throw CG_FATAL("HepMC3:fillEvent") << "Other particle requested! Not yet implemented!";
      } break;
    }
  }
  add_vertex(vertex_beam1);
  add_vertex(vertex_beam2);
  add_vertex(vertex_central_system);
  for (const auto& decay_vertex : decay_vertices)
    add_vertex(decay_vertex);
}

CepGenEvent::operator cepgen::Event() const {
  cepgen::Event event;
  const auto convert_particle = [](const GenParticle& hepmc_particle, const cepgen::Particle::Role cepgen_role) {
    const auto convert_momentum = [](const FourVector& momentum) -> cepgen::Momentum {
      return cepgen::Momentum::fromPxPyPzE(momentum.px(), momentum.py(), momentum.pz(), momentum.e());
    };
    auto cepgen_particle =
        cepgen::Particle(cepgen_role, 0, static_cast<cepgen::Particle::Status>(hepmc_particle.status()));
    cepgen_particle.setPdgId(hepmc_particle.pdg_id());
    cepgen_particle.setMomentum(convert_momentum(hepmc_particle.momentum()));
    return cepgen_particle;
  };

  auto incoming_beam1 = beams().at(0), incoming_beam2 = beams().at(1);
  std::vector<int> beam_vtx_ids;
  for (const auto& vertex : vertices()) {
    if (vertex->particles_in().size() == 1) {  // single particle decay
      auto incoming_role = cepgen::Particle::Role::UnknownRole, intermediate_role = cepgen::Particle::Role::UnknownRole,
           outgoing_role = cepgen::Particle::Role::UnknownRole;
      size_t id_beam_in = 0;
      if (auto& hepmc_particle = vertex->particles_in().at(0); hepmc_particle) {  // retrieve the origin particle
        if (hepmc_particle->id() == incoming_beam1->id()) {                       // positive-z beam
          incoming_role = cepgen::Particle::Role::IncomingBeam1;
          intermediate_role = cepgen::Particle::Role::Parton1;
          outgoing_role = cepgen::Particle::Role::OutgoingBeam1;
        } else if (hepmc_particle->id() == incoming_beam2->id()) {  // negative-z system
          incoming_role = cepgen::Particle::Role::IncomingBeam2;
          intermediate_role = cepgen::Particle::Role::Parton2;
          outgoing_role = cepgen::Particle::Role::OutgoingBeam2;
        }
        auto cepgen_particle =
            convert_particle(*hepmc_particle, incoming_role).setStatus(cepgen::Particle::Status::PrimordialIncoming);
        event.addParticle(cepgen_particle);
        id_beam_in = cepgen_particle.id();
      }
      if (vertex->particles_out_size() >= 2) {  //FIXME handle cases with multiple partons?
        size_t num_outgoing_particles = 0;
        for (auto outgoing_particle : vertex->particles_out()) {
          auto cepgen_particle =
              convert_particle(*outgoing_particle, num_outgoing_particles == 0 ? intermediate_role : outgoing_role);
          cepgen_particle.setStatus(num_outgoing_particles == 0 ? cepgen::Particle::Status::Incoming
                                                                : cepgen::Particle::Status::Unfragmented);
          auto added_particle = event.addParticle(cepgen_particle);
          event[id_beam_in].addChild(added_particle);
          ++num_outgoing_particles;
        }
      }
      beam_vtx_ids.emplace_back(vertex->id());
    }
  }

  auto cepgen_intermediate =
      cepgen::Particle(cepgen::Particle::Role::Intermediate, 0, cepgen::Particle::Status::Propagator);
  auto &parton1 = event.oneWithRole(cepgen::Particle::Role::Parton1),
       &parton2 = event.oneWithRole(cepgen::Particle::Role::Parton2);
  cepgen_intermediate.setMomentum(parton1.momentum() + parton2.momentum(), true);
  parton1.addChild(cepgen_intermediate);
  parton2.addChild(cepgen_intermediate);
  event.addParticle(cepgen_intermediate);

  for (const auto& vtx : vertices()) {
    if (cepgen::utils::contains(beam_vtx_ids, vtx->id()))
      continue;
    for (auto outgoing_particles : vtx->particles_out()) {
      auto cepgen_particle = convert_particle(*outgoing_particles, cepgen::Particle::Role::CentralSystem);
      auto added_particle = event.addParticle(cepgen_particle);
      event.oneWithRole(cepgen::Particle::Role::Intermediate).addChild(added_particle);
    }
  }
  return event;
}

bool CepGenEvent::merge(cepgen::Event& event) const {
  // set of sanity checks to ensure the compatibility between the HepMC and CepGen event records
  for (const auto role : {cepgen::Particle::Role::IncomingBeam1,
                          cepgen::Particle::Role::IncomingBeam2,
                          cepgen::Particle::Role::Intermediate})
    if (!event.hasRole(role) || event[role].size() != 1) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Event record has no unique particle with " << role << " role.";
      return false;
    }
  if (vertices().size() < 3) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Failed to retrieve the three primordial vertices in event.";
    return false;
  }
  if (const auto vertex_incoming_beam1 = vertices().at(0);
      !vertex_incoming_beam1 || vertex_incoming_beam1->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particles multiplicity: found "
                                         << vertex_incoming_beam1->particles_in().size()
                                         << ", expecting one single beam particle.";
    return false;
  } else {
    const auto incoming_beam1 = vertex_incoming_beam1->particles_in().at(0);
    if (const auto& cepgen_incoming_beam1 = event.oneWithRole(cepgen::Particle::Role::IncomingBeam1);
        std::fabs(incoming_beam1->momentum().x() - cepgen_incoming_beam1.momentum().px()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().y() - cepgen_incoming_beam1.momentum().py()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().z() - cepgen_incoming_beam1.momentum().pz()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().t() - cepgen_incoming_beam1.momentum().energy()) > kTolerance) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particle kinematics.";
      return false;
    }
  }
  if (const auto vertex_incoming_beam2 = vertices().at(1);
      !vertex_incoming_beam2 || vertex_incoming_beam2->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particles multiplicity: found "
                                         << vertex_incoming_beam2->particles_in().size()
                                         << ", expecting one single beam particle.";
    return false;
  } else {
    const auto incoming_beam2 = vertex_incoming_beam2->particles_in().at(0);
    if (const auto& cepgen_incoming_beam2 = event.oneWithRole(cepgen::Particle::Role::IncomingBeam2);
        std::fabs(incoming_beam2->momentum().x() - cepgen_incoming_beam2.momentum().px()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().y() - cepgen_incoming_beam2.momentum().py()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().z() - cepgen_incoming_beam2.momentum().pz()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().t() - cepgen_incoming_beam2.momentum().energy()) > kTolerance) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particle kinematics.";
      return false;
    }
  }
  const auto vertex_central_system = vertices().at(2);

  const auto convert_momentum = [](const ConstGenParticlePtr& hepmc_particle) {
    return cepgen::Momentum::fromPxPyPzE(hepmc_particle->momentum().x(),
                                         hepmc_particle->momentum().y(),
                                         hepmc_particle->momentum().z(),
                                         hepmc_particle->momentum().t());
  };
  // particles already present in the CepGen event record; anything else was produced by the external algorithm
  std::unordered_set<const GenParticle*> known_particles;
  for (const auto& [_, hepmc_particle] : cepgen_id_vs_hepmc_particle_)
    known_particles.insert(hepmc_particle.get());

  // helper function to store a new particle into the CepGen event content
  const auto add_particle = [&](const ConstGenParticlePtr& hepmc_particle, int mother_id, cepgen::Particle::Role role) {
    auto cepgen_particle = cepgen::Particle(role, 0, cepgen::Particle::Status::FinalState);
    cepgen_particle.setIntegerPdgId(hepmc_particle->pdg_id());
    cepgen_particle.setMomentum(convert_momentum(hepmc_particle), true);
    const auto particle_id = event.addParticle(cepgen_particle).get().id();
    event[mother_id].addChild(event[particle_id]);
    return particle_id;
  };
  // helper function to browse particles decay products and store the new ones into the CepGen event content
  std::function<void(const ConstGenParticlePtr&, int)> browse_children = [&](const ConstGenParticlePtr& hepmc_particle,
                                                                             int mother_id) {
    const auto children = hepmc_particle->children();
    for (const auto& hepmc_child : children)
      if (known_particles.insert(hepmc_child.get()).second)  // this child is not yet in the CepGen event
        browse_children(hepmc_child,
                        add_particle(hepmc_child, mother_id, event[mother_id].role()));  // launch recursion
    if (!children.empty())
      event[mother_id].setStatus(cepgen::Particle::Status::Propagator);
  };

  std::vector<int> central_system_ids;  // freeze the list of "primordial" central system particles
  for (const auto& central_particle : event[cepgen::Particle::Role::CentralSystem])
    central_system_ids.emplace_back(central_particle.get().id());
  for (const auto id : central_system_ids) {
    const auto hepmc_particle = cepgen_id_vs_hepmc_particle_.find(id);
    if (hepmc_particle == cepgen_id_vs_hepmc_particle_.end() || !hepmc_particle->second->in_event())
      continue;  // particle was not propagated to the HepMC3 record (e.g. a disconnected line)
    event[id].setMomentum(convert_momentum(hepmc_particle->second), true);  // e.g. recoil from a radiated photon
    browse_children(hepmc_particle->second, id);
  }
  // particles directly attached to the central system vertex (e.g. a photon radiated in the two-parton system decay)
  for (const auto& central_particle : vertex_central_system->particles_out())
    if (known_particles.insert(central_particle.get()).second)
      browse_children(central_particle,
                      add_particle(central_particle,
                                   event.oneWithRole(cepgen::Particle::Role::Intermediate).id(),
                                   cepgen::Particle::Role::CentralSystem));
  return true;
}

void CepGenEvent::dump() const {
  CG_LOG.log([&](auto& log) {
    log << "HepMC3::CepGenEvent\n"
        << " Attributes:\n";
    for (const auto& [name, _] : attributes())
      log << " * " << name << " = " << attribute_as_string(name) << "\n";
    log << " Vertices:";
    for (const auto& vertex : vertices()) {
      FourVector incoming_momentum, outgoing_momentum;
      log << "\n  * vertex#" << -vertex->id() << " (status: " << vertex->status() << ")" << "\n     in: ";
      for (const auto& incoming_particle : vertex->particles_in())
        log << "\n      * " << incoming_particle->pdg_id() << " (status: " << incoming_particle->status()
            << "): " << incoming_particle->momentum(),
            incoming_momentum += incoming_particle->momentum();
      log << "\n     total: " << incoming_momentum << "\n     out:";
      for (const auto& outgoing_particle : vertex->particles_out())
        log << "\n      * " << outgoing_particle->pdg_id() << " (status: " << outgoing_particle->status()
            << "): " << outgoing_particle->momentum(),
            outgoing_momentum += outgoing_particle->momentum();
      const auto momentum_imbalance = incoming_momentum - outgoing_momentum;
      log << "\n     total: " << outgoing_momentum << "\n    (im)balance: " << momentum_imbalance
          << " (norm: " << momentum_imbalance.length() << ").";
    }
    log << "\n" << std::string(70, '-');
  });
}

namespace HepMC3 {
  std::ostream& operator<<(std::ostream& os, const FourVector& vec) {
    return os << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << "; " << vec.t() << ")";
  }
}  // namespace HepMC3
