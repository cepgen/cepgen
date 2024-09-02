/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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
#include <HepMC/Version.h>

#include <list>
#include <numeric>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/HepMC2Wrapper/HepMC2EventInterface.h"

namespace HepMC {
  CepGenEvent::CepGenEvent(const cepgen::Event& evt) : GenEvent(Units::GEV, Units::MM) {
    set_alphaQCD(evt.metadata("alphaS"));
    set_alphaQED(evt.metadata("alphaEM"));

    weights().push_back(1.);  // unweighted events

    // filling the particles content
    const FourVector origin(0., 0., 0., 0.);
    int cm_id = 0;

    auto convert_particle = [](const cepgen::Particle& cg_part) -> GenParticle* {
      const auto cg_mom = cg_part.momentum();
      auto* part = new GenParticle(FourVector(cg_mom.px(), cg_mom.py(), cg_mom.pz(), cg_mom.energy()),
                                   cg_part.integerPdgId(),
                                   (int)cg_part.status());
      part->set_generated_mass(cepgen::PDG::get().mass(cg_part.pdgId()));
      return part;
    };

    auto v1 = new GenVertex(origin), v2 = new GenVertex(origin), vcm = new GenVertex(origin);
    unsigned short idx = 1;
    for (const auto& part_orig : evt.particles()) {
      auto* part = convert_particle(part_orig);
      part->suggest_barcode(idx);
      assoc_map_[idx] = part;

      switch (part_orig.role()) {
        case cepgen::Particle::Role::IncomingBeam1:
          v1->add_particle_in(part);
          break;
        case cepgen::Particle::Role::IncomingBeam2:
          v2->add_particle_in(part);
          break;
        case cepgen::Particle::Role::OutgoingBeam1:
          v1->add_particle_out(part);
          break;
        case cepgen::Particle::Role::OutgoingBeam2:
          v2->add_particle_out(part);
          break;
        case cepgen::Particle::Role::Parton1:
          v1->add_particle_out(part);
          vcm->add_particle_in(part);
          break;
        case cepgen::Particle::Role::Parton2:
          v2->add_particle_out(part);
          vcm->add_particle_in(part);
          break;
        case cepgen::Particle::Role::Intermediate:
          // skip the two-parton system and propagate the parentage
          cm_id = idx;
          continue;
        case cepgen::Particle::Role::CentralSystem:
        default: {
          const auto& moth = part_orig.mothers();
          if (moth.empty())
            // skip disconnected lines
            continue;
          // get mother(s) id(s)
          const short m1 = *moth.begin();
          const short m2 = moth.size() > 1 ? *moth.rbegin() : -1;
          // check if particle is connected to the two-parton system
          if (m1 == cm_id || (m2 >= 0 && (m1 < cm_id && cm_id <= m2)))  // also supports range
            vcm->add_particle_out(part);
          // if part of the decay chain of central system, find parents
          else if (assoc_map_.count(m1) != 0) {
            auto vprod = assoc_map_.at(m1)->end_vertex();
            std::list<short> ids{m1};  // list of mother particles
            if (assoc_map_.count(m2) != 0 && m2 > m1) {
              ids.resize(m2 - m1 + 1);
              std::iota(ids.begin(), ids.end(), m1);
            }
            if (!vprod) {
              vprod = new GenVertex();
              for (const auto& id : ids)
                vprod->add_particle_in(assoc_map_.at(id));
              add_vertex(vprod);
            }
            vprod->add_particle_out(part);
          } else {
            if (v1)
              delete v1;
            if (v2)
              delete v2;
            if (vcm)
              delete vcm;
            throw CG_FATAL("HepMC2:fillEvent") << "Other particle requested! Not yet implemented!";
          }
        } break;
      }
      idx++;
    }
    add_vertex(v1);
    add_vertex(v2);
    add_vertex(vcm);

    if (v1->particles_in_size() > 0 && v2->particles_in_size() > 0)
      set_beam_particles(*v1->particles_in_const_begin(), *v2->particles_in_const_begin());
    if (evt.hasRole(cepgen::Particle::Role::Intermediate))
      set_event_scale(evt.oneWithRole(cepgen::Particle::Role::Intermediate).momentum().mass());
    set_signal_process_vertex(vcm);
  }

  CepGenEvent::operator cepgen::Event() const {
    cepgen::Event evt;
    auto convert_particle = [](const GenParticle& part,
                               const cepgen::Particle::Role& role =
                                   cepgen::Particle::Role::UnknownRole) -> cepgen::Particle {
      auto convert_momentum = [](const FourVector& mom) -> cepgen::Momentum {
        return cepgen::Momentum::fromPxPyPzE(mom.px(), mom.py(), mom.pz(), mom.e());
      };
      auto cg_part = cepgen::Particle(role, 0, (cepgen::Particle::Status)part.status());
      cg_part.setPdgId((long)part.pdg_id());
      cg_part.setMomentum(convert_momentum(part));
      return cg_part;
    };

    auto [ip1, ip2] = beam_particles();
    std::unordered_map<size_t, size_t> h_to_cg;
    std::vector<int> beam_vtx_barcodes;
    for (auto it_vtx = vertices_begin(); it_vtx != vertices_end(); ++it_vtx) {
      if ((*it_vtx)->particles_in_size() == 1) {
        auto role1 = cepgen::Particle::Role::UnknownRole, role2 = role1, role3 = role1;
        size_t id_beam_in = 0;
        if (auto* part = *(*it_vtx)->particles_in_const_begin(); part) {
          if (part->barcode() == ip1->barcode()) {
            role1 = cepgen::Particle::Role::IncomingBeam1;
            role2 = cepgen::Particle::Role::Parton1;
            role3 = cepgen::Particle::Role::OutgoingBeam1;
          } else if (part->barcode() == ip2->barcode()) {
            role1 = cepgen::Particle::Role::IncomingBeam2;
            role2 = cepgen::Particle::Role::Parton2;
            role3 = cepgen::Particle::Role::OutgoingBeam2;
          }
          auto cg_part = convert_particle(*part, role1);
          cg_part.setStatus(cepgen::Particle::Status::PrimordialIncoming);
          evt.addParticle(cg_part);
          h_to_cg[part->barcode()] = cg_part.id();
          id_beam_in = cg_part.id();
        }
        if ((*it_vtx)->particles_out_size() == 2) {  //FIXME handle cases with multiple partons?
          size_t num_op = 0;
          for (auto it_op = (*it_vtx)->particles_out_const_begin(); it_op != (*it_vtx)->particles_out_const_end();
               ++it_op, ++num_op) {
            auto cg_part = convert_particle(*(*it_op), num_op == 0 ? role2 : role3);
            cg_part.setStatus(num_op == 0 ? cepgen::Particle::Status::Incoming
                                          : cepgen::Particle::Status::Unfragmented);
            cg_part.addMother(evt[id_beam_in]);
            evt.addParticle(cg_part);
            h_to_cg[(*it_op)->barcode()] = cg_part.id();
          }
        }
        beam_vtx_barcodes.emplace_back((*it_vtx)->barcode());
      }
    }

    auto cg_interm = cepgen::Particle(cepgen::Particle::Role::Intermediate, 0, cepgen::Particle::Status::Propagator);
    auto &part1 = evt.oneWithRole(cepgen::Particle::Role::Parton1),
         &part2 = evt.oneWithRole(cepgen::Particle::Role::Parton2);
    cg_interm.setMomentum(part1.momentum() + part2.momentum(), true);
    cg_interm.addMother(part1);
    cg_interm.addMother(part2);
    evt.addParticle(cg_interm);

    for (auto it_vtx = vertices_begin(); it_vtx != vertices_end(); ++it_vtx) {
      if (cepgen::utils::contains(beam_vtx_barcodes, (*it_vtx)->barcode()))
        continue;
      if ((*it_vtx)->barcode() == signal_process_vertex()->barcode()) {
        for (auto it_op = (*it_vtx)->particles_out_const_begin(); it_op != (*it_vtx)->particles_out_const_end();
             ++it_op) {
          auto cg_part = convert_particle(*(*it_op), cepgen::Particle::Role::CentralSystem);
          cg_part.addMother(evt.oneWithRole(cepgen::Particle::Role::Intermediate));
          evt.addParticle(cg_part);
          h_to_cg[(*it_op)->barcode()] = cg_part.id();
        }
      } else {
        (*it_vtx)->print(std::cout);
        throw CG_FATAL("CepGenEvent") << "Not yet supporting secondary decay of central system.";
      }
    }
    return evt;
  }
}  // namespace HepMC
