/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/SimpleVector.h>

#include <list>
#include <numeric>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/String.h"
#include "CepGenHepMC2/CepGenEvent.h"

using namespace HepMC;

CepGenEvent::CepGenEvent(const cepgen::Event& event) : GenEvent(Units::GEV, Units::MM) {
  set_alphaQCD(event.metadata("alphaS"));
  set_alphaQED(event.metadata("alphaEM"));

  weights().push_back(1.);  // unweighted events

  const FourVector origin(0., 0., 0., 0.);
  int central_system_id = 0;

  auto convert_particle = [](const cepgen::Particle& cg_part) -> GenParticle* {
    const auto cg_mom = cg_part.momentum();
    auto* part = new GenParticle(FourVector(cg_mom.px(), cg_mom.py(), cg_mom.pz(), cg_mom.energy()),
                                 cg_part.integerPdgId(),
                                 static_cast<int>(cg_part.status()));
    part->set_generated_mass(cepgen::PDG::get().mass(cg_part.pdgId()));
    return part;
  };

  auto vertex_beam1 = std::make_unique<GenVertex>(origin), vertex_beam2 = std::make_unique<GenVertex>(origin),
       vertex_central_system = std::make_unique<GenVertex>(origin);
  unsigned short idx = 1;
  for (const auto& cepgen_particle : event.particles()) {  // filling the particles content
    auto* hepmc_particle = convert_particle(cepgen_particle);
    hepmc_particle->suggest_barcode(idx);
    cepgen_id_vs_hepmc_particle_[idx] = hepmc_particle;

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
        central_system_id = idx;
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
            production_vertex = new GenVertex();
            for (const auto& id : ids)
              production_vertex->add_particle_in(cepgen_id_vs_hepmc_particle_.at(id));
            add_vertex(production_vertex);
          }
          production_vertex->add_particle_out(hepmc_particle);
        } else
          throw CG_FATAL("HepMC2:fillEvent") << "Other particle requested! Not yet implemented!";
      } break;
    }
    idx++;
  }
  auto *v1 = vertex_beam1.release(), *v2 = vertex_beam2.release(), *vcm = vertex_central_system.release();
  add_vertex(v1);
  add_vertex(v2);
  add_vertex(vcm);
  if (v1->particles_in_size() > 0 && v2->particles_in_size() > 0)
    set_beam_particles(*v1->particles_in_const_begin(), *v2->particles_in_const_begin());
  if (event.hasRole(cepgen::Particle::Role::Intermediate))
    set_event_scale(event.oneWithRole(cepgen::Particle::Role::Intermediate).momentum().mass());
  set_signal_process_vertex(vcm);
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

  auto [incoming_beam1, incoming_beam2] = beam_particles();
  std::unordered_map<size_t, size_t> hepmc_to_cepgen;
  std::vector<int> beam_vtx_barcodes;
  for (auto it_vertices = vertices_begin(); it_vertices != vertices_end(); ++it_vertices) {
    if (const auto* vertex = *it_vertices; vertex->particles_in_size() == 1) {
      auto incoming_role = cepgen::Particle::Role::UnknownRole, intermediate_role = cepgen::Particle::Role::UnknownRole,
           outgoing_role = cepgen::Particle::Role::UnknownRole;
      size_t id_beam_in = 0;
      if (const auto* hepmc_particle = *vertex->particles_in_const_begin(); hepmc_particle) {
        if (hepmc_particle->barcode() == incoming_beam1->barcode()) {  // first branch: positive-z system
          incoming_role = cepgen::Particle::Role::IncomingBeam1;
          intermediate_role = cepgen::Particle::Role::Parton1;
          outgoing_role = cepgen::Particle::Role::OutgoingBeam1;
        } else if (hepmc_particle->barcode() == incoming_beam2->barcode()) {  // second branch: negative-z system
          incoming_role = cepgen::Particle::Role::IncomingBeam2;
          intermediate_role = cepgen::Particle::Role::Parton2;
          outgoing_role = cepgen::Particle::Role::OutgoingBeam2;
        }
        auto cepgen_particle = convert_particle(*hepmc_particle, incoming_role);
        cepgen_particle.setStatus(cepgen::Particle::Status::PrimordialIncoming);
        event.addParticle(cepgen_particle);
        hepmc_to_cepgen[hepmc_particle->barcode()] = cepgen_particle.id();
        id_beam_in = cepgen_particle.id();
      }
      if (vertex->particles_out_size() >= 2) {  //FIXME handle cases with multiple partons?
        size_t num_outgoing_particles = 0;
        for (auto it_outgoing_particles = vertex->particles_out_const_begin();
             it_outgoing_particles != vertex->particles_out_const_end();
             ++it_outgoing_particles, ++num_outgoing_particles) {
          const auto* outgoing_particle = *it_outgoing_particles;
          auto cepgen_particle =
              convert_particle(*outgoing_particle, num_outgoing_particles == 0 ? intermediate_role : outgoing_role);
          cepgen_particle.setStatus(num_outgoing_particles == 0 ? cepgen::Particle::Status::Incoming
                                                                : cepgen::Particle::Status::Unfragmented);
          cepgen_particle.addMother(event[id_beam_in]);
          event.addParticle(cepgen_particle);
          hepmc_to_cepgen[outgoing_particle->barcode()] = cepgen_particle.id();  // bookkeeping of the two IDs
        }
      }
      beam_vtx_barcodes.emplace_back(vertex->barcode());
    }
  }

  auto cepgen_intermediate =
      cepgen::Particle(cepgen::Particle::Role::Intermediate, 0, cepgen::Particle::Status::Propagator);
  auto &parton1 = event.oneWithRole(cepgen::Particle::Role::Parton1),
       &parton2 = event.oneWithRole(cepgen::Particle::Role::Parton2);
  cepgen_intermediate.setMomentum(parton1.momentum() + parton2.momentum(), true);
  cepgen_intermediate.addMother(parton1);
  cepgen_intermediate.addMother(parton2);
  event.addParticle(cepgen_intermediate);

  for (auto it_vertices = vertices_begin(); it_vertices != vertices_end(); ++it_vertices) {
    const auto* vertex = *it_vertices;
    if (cepgen::utils::contains(beam_vtx_barcodes, vertex->barcode()))
      continue;
    if (vertex->barcode() == signal_process_vertex()->barcode()) {
      for (auto it_outgoing_particles = vertex->particles_out_const_begin();
           it_outgoing_particles != vertex->particles_out_const_end();
           ++it_outgoing_particles) {
        auto* outgoing_particle = *it_outgoing_particles;
        auto cepgen_particle = convert_particle(*outgoing_particle, cepgen::Particle::Role::CentralSystem);
        cepgen_particle.addMother(event.oneWithRole(cepgen::Particle::Role::Intermediate));
        event.addParticle(cepgen_particle);
        hepmc_to_cepgen[outgoing_particle->barcode()] = cepgen_particle.id();
      }
    } else
      throw CG_FATAL("CepGenEvent").log([&vertex](auto& log) {
        log << "Not yet supporting secondary decay of central system. Problematic vertex:\n";
        vertex->print(log.stream());
      });
  }
  return event;
}
