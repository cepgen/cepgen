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
#include "CepGenAddOns/HepMC2Wrapper/HepMC2EventInterface.h"

namespace HepMC {
  CepGenEvent::CepGenEvent(const cepgen::Event& evt) : GenEvent(Units::GEV, Units::MM) {
    set_alphaQCD(evt.alpha_s);
    set_alphaQED(evt.alpha_em);

    weights().push_back(1.);  // unweighted events

    // filling the particles content
    const FourVector origin(0., 0., 0., 0.);
    int cm_id = 0;

    auto v1 = new GenVertex(origin), v2 = new GenVertex(origin), vcm = new GenVertex(origin);
    unsigned short idx = 1;
    for (const auto& part_orig : evt.particles()) {
      const auto& mom_orig = part_orig.momentum();
      FourVector pmom(mom_orig.px(), mom_orig.py(), mom_orig.pz(), mom_orig.energy());
      auto part = new GenParticle(pmom, part_orig.integerPdgId(), (int)part_orig.status());
      part->set_generated_mass(cepgen::PDG::get().mass(part_orig.pdgId()));
      part->suggest_barcode(idx);
      assoc_map_[idx].reset(part);

      switch (part_orig.role()) {
        case cepgen::Particle::IncomingBeam1:
          v1->add_particle_in(part);
          break;
        case cepgen::Particle::IncomingBeam2:
          v2->add_particle_in(part);
          break;
        case cepgen::Particle::OutgoingBeam1:
          v1->add_particle_out(part);
          break;
        case cepgen::Particle::OutgoingBeam2:
          v2->add_particle_out(part);
          break;
        case cepgen::Particle::Parton1:
          v1->add_particle_out(part);
          vcm->add_particle_in(part);
          break;
        case cepgen::Particle::Parton2:
          v2->add_particle_out(part);
          vcm->add_particle_in(part);
          break;
        case cepgen::Particle::Intermediate:
          // skip the two-parton system and propagate the parentage
          cm_id = idx;
          continue;
        case cepgen::Particle::CentralSystem:
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
                vprod->add_particle_in(assoc_map_.at(id).get());
              add_vertex(vprod);
            }
            vprod->add_particle_out(part);
          } else
            throw CG_FATAL("HepMC2:fillEvent") << "Other particle requested! Not yet implemented!";
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
}  // namespace HepMC
