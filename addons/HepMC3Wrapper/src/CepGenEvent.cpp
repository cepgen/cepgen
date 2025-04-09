/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <list>
#include <numeric>

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

  weights().push_back(1.);  // unweighted events

  // filling the particles content
  const FourVector origin(0., 0., 0., 0.);
  int central_system_id = 0;

  auto vertex_beam1 = make_shared<GenVertex>(origin), vertex_beam2 = make_shared<GenVertex>(origin),
       vertex_central_system = make_shared<GenVertex>(origin);
  size_t idx = 0;
  for (const auto& cepgen_particle : event.particles()) {
    const auto& mom_orig = cepgen_particle.momentum();
    FourVector momentum(mom_orig.px(), mom_orig.py(), mom_orig.pz(), mom_orig.energy());
    auto hepmc_particle =
        make_shared<GenParticle>(momentum, cepgen_particle.integerPdgId(), static_cast<int>(cepgen_particle.status()));
    hepmc_particle->set_generated_mass(cepgen::PDG::get().mass(cepgen_particle.pdgId()));
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
        const auto& mother = cepgen_particle.mothers();
        if (mother.empty())
          continue;  // skip disconnected lines
        // check if particle is connected to the two-parton system
        if (const short m1 = *mother.begin(), m2 = mother.size() > 1 ? *mother.rbegin() : -1;  // get mother(s) id(s)
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
            add_vertex(production_vertex);
          }
          production_vertex->add_particle_out(hepmc_particle);
        } else
          throw CG_FATAL("HepMC3:fillEvent") << "Other particle requested! Not yet implemented!";
      } break;
    }
    idx++;
  }
  add_vertex(vertex_beam1);
  add_vertex(vertex_beam2);
  add_vertex(vertex_central_system);
}

CepGenEvent::operator cepgen::Event() const {
  cepgen::Event event;
  auto convert_particle = [](const GenParticle& hepmc_particle,
                             const cepgen::Particle::Role& cepgen_role =
                                 cepgen::Particle::Role::UnknownRole) -> cepgen::Particle {
    auto convert_momentum = [](const FourVector& mom) -> cepgen::Momentum {
      return cepgen::Momentum::fromPxPyPzE(mom.px(), mom.py(), mom.pz(), mom.e());
    };
    auto cepgen_particle =
        cepgen::Particle(cepgen_role, 0, static_cast<cepgen::Particle::Status>(hepmc_particle.status()));
    cepgen_particle.setPdgId(hepmc_particle.pdg_id());
    cepgen_particle.setMomentum(convert_momentum(hepmc_particle.momentum()));
    return cepgen_particle;
  };

  auto ip1 = beams().at(0), ip2 = beams().at(1);
  std::unordered_map<size_t, size_t> hepmc_to_cepgen;
  std::vector<int> beam_vtx_ids;
  for (const auto& vertex : vertices()) {
    if (vertex->particles_in().size() == 1) {
      auto role1 = cepgen::Particle::Role::UnknownRole, role2 = role1, role3 = role1;
      size_t id_beam_in = 0;
      if (auto& part = vertex->particles_in().at(0); part) {
        if (part->id() == ip1->id()) {  // first branch: positive-z system
          role1 = cepgen::Particle::Role::IncomingBeam1;
          role2 = cepgen::Particle::Role::Parton1;
          role3 = cepgen::Particle::Role::OutgoingBeam1;
        } else if (part->id() == ip2->id()) {  // second branch: negative-z system
          role1 = cepgen::Particle::Role::IncomingBeam2;
          role2 = cepgen::Particle::Role::Parton2;
          role3 = cepgen::Particle::Role::OutgoingBeam2;
        }
        auto cepgen_particle = convert_particle(*part, role1);
        cepgen_particle.setStatus(cepgen::Particle::Status::PrimordialIncoming);
        event.addParticle(cepgen_particle);
        hepmc_to_cepgen[part->id()] = cepgen_particle.id();
        id_beam_in = cepgen_particle.id();
      }
      if (vertex->particles_out_size() == 2) {  //FIXME handle cases with multiple partons?
        size_t num_outgoing_particles = 0;
        for (auto outgoing_particle : vertex->particles_out()) {
          auto cepgen_particle = convert_particle(*outgoing_particle, num_outgoing_particles == 0 ? role2 : role3);
          cepgen_particle.setStatus(num_outgoing_particles == 0 ? cepgen::Particle::Status::Incoming
                                                                : cepgen::Particle::Status::Unfragmented);
          cepgen_particle.addMother(event[id_beam_in]);
          event.addParticle(cepgen_particle);
          hepmc_to_cepgen[outgoing_particle->id()] = cepgen_particle.id();  // bookkeeping of the two IDs
          ++num_outgoing_particles;
        }
      }
      beam_vtx_ids.emplace_back(vertex->id());
    }
  }

  auto cg_intermediate =
      cepgen::Particle(cepgen::Particle::Role::Intermediate, 0, cepgen::Particle::Status::Propagator);
  auto &part1 = event.oneWithRole(cepgen::Particle::Role::Parton1),
       &part2 = event.oneWithRole(cepgen::Particle::Role::Parton2);
  cg_intermediate.setMomentum(part1.momentum() + part2.momentum(), true);
  cg_intermediate.addMother(part1);
  cg_intermediate.addMother(part2);
  event.addParticle(cg_intermediate);

  for (const auto& vtx : vertices()) {
    if (cepgen::utils::contains(beam_vtx_ids, vtx->id()))
      continue;
    for (auto op : vtx->particles_out()) {
      auto cg_part = convert_particle(*op, cepgen::Particle::Role::CentralSystem);
      cg_part.addMother(event.oneWithRole(cepgen::Particle::Role::Intermediate));
      event.addParticle(cg_part);
      hepmc_to_cepgen[op->id()] = cg_part.id();
    }
  }
  return event;
}

void CepGenEvent::merge(cepgen::Event& event) const {
  // set of sanity checks to perform on the HepMC event content
  if (vertices().size() < 3) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Failed to retrieve the three primordial vertices in event.";
    return;
  }
  if (const auto vertex_incoming_beam1 = vertices().at(0); vertex_incoming_beam1->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particles multiplicity: found "
                                         << vertex_incoming_beam1->particles_in().size() << ", expecting one.";
    return;
  } else {  // set of sanity checks to ensure the compatibility between the HepMC and CepGen event records
    const auto incoming_beam1 = vertex_incoming_beam1->particles_in().at(0);
    if (const auto& cepgen_incoming_beam1 = event.oneWithRole(cepgen::Particle::Role::IncomingBeam1);
        std::fabs(incoming_beam1->momentum().x() - cepgen_incoming_beam1.momentum().px()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().y() - cepgen_incoming_beam1.momentum().py()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().z() - cepgen_incoming_beam1.momentum().pz()) > kTolerance ||
        std::fabs(incoming_beam1->momentum().t() - cepgen_incoming_beam1.momentum().energy()) > kTolerance) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particle kinematics.";
      return;
    }
  }
  if (const auto vertex_incoming_beam2 = vertices().at(1); vertex_incoming_beam2->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particles multiplicity: found "
                                         << vertex_incoming_beam2->particles_in().size() << ", expecting one.";
    return;
  } else {
    const auto incoming_beam2 = vertex_incoming_beam2->particles_in().at(0);
    if (const auto& cepgen_incoming_beam2 = event.oneWithRole(cepgen::Particle::Role::IncomingBeam2);
        std::fabs(incoming_beam2->momentum().x() - cepgen_incoming_beam2.momentum().px()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().y() - cepgen_incoming_beam2.momentum().py()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().z() - cepgen_incoming_beam2.momentum().pz()) > kTolerance ||
        std::fabs(incoming_beam2->momentum().t() - cepgen_incoming_beam2.momentum().energy()) > kTolerance) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particle kinematics.";
      return;
    }
  }
  const auto vertex_central_system = vertices().at(2);
  const auto central_system = event[cepgen::Particle::Role::CentralSystem];
  if (central_system.size() != vertex_central_system->particles_out().size()) {
    CG_ERROR("HepMC3:CepGenEvent:merge")
        << "Central system particles multiplicities differ between CepGen and HepMC3 event records.";
    return;
  }
  // freeze the "primordial" central system size
  const auto cs_size = central_system.size();

  // helper function to browse particles decay products and store them into the CepGen event content
  std::function<void(const ConstGenParticlePtr&, cepgen::ParticleRef)> browse_children =
      [&](const ConstGenParticlePtr& hepmc_particle, cepgen::ParticleRef cepgen_particle) {
        if (hepmc_particle->children().empty())
          return;
        cepgen_particle.get().setStatus(cepgen::Particle::Status::Propagator);
        for (const auto& hepmc_child : hepmc_particle->children()) {
          cepgen::Particle cepgen_child(cepgen_particle.get().role(), 0);
          cepgen_child.setPdgId(hepmc_child->pdg_id());
          const auto& hepmc_child_momentum = hepmc_child->momentum();
          cepgen_child.setStatus(cepgen::Particle::Status::FinalState);
          cepgen_child.setMomentum(cepgen::Momentum::fromPxPyPzE(
              hepmc_child_momentum.x(), hepmc_child_momentum.y(), hepmc_child_momentum.z(), hepmc_child_momentum.t()));
          cepgen_child.addMother(cepgen_particle);
          browse_children(hepmc_child, event.addParticle(cepgen_child));  // launch recursion
        }
      };

  for (size_t icg = 0; icg < cs_size; ++icg) {  // try to find the associated CepGen event particle
    const auto central_particle_momentum = central_system.at(icg).get().momentum().p();
    // loop over the central system particles
    for (const auto& central_particle : vertex_central_system->particles_out()) {
      if (std::fabs(central_particle_momentum - central_particle->momentum().length()) > kTolerance)
        continue;
      browse_children(central_particle, central_system[icg]);
      break;  // found the association between the HepMC and CepGen particles kinematics
    }
  }
}

void CepGenEvent::dump() const {
  CG_LOG.log([&](auto& log) {
    log << "HepMC3::CepGenEvent\n"
        << " Attributes:\n";
    for (const auto& attr : {"AlphaEM", "AlphaQCD"})
      log << " * " << attr << " = " << attribute_as_string(attr) << "\n";
    log << " Vertices:";
    for (const auto& vtxPtr : vertices()) {
      FourVector in_sys, out_sys;
      log << "\n  * vertex#" << -vtxPtr->id() << " (status: " << vtxPtr->status() << ")"
          << "\n     in: ";
      for (const auto& ipPtr : vtxPtr->particles_in())
        log << "\n      * " << ipPtr->pdg_id() << " (status: " << ipPtr->status() << "): " << ipPtr->momentum(),
            in_sys += ipPtr->momentum();
      log << "\n     total: " << in_sys << "\n     out:";
      for (const auto& opPtr : vtxPtr->particles_out())
        log << "\n      * " << opPtr->pdg_id() << " (status: " << opPtr->status() << "): " << opPtr->momentum(),
            out_sys += opPtr->momentum();
      const auto momentum_imbalance = in_sys - out_sys;
      log << "\n     total: " << out_sys << "\n    (im)balance: " << momentum_imbalance
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
