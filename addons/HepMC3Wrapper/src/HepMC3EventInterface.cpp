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
#include "CepGenHepMC3/HepMC3EventInterface.h"

using namespace HepMC3;

CepGenEvent::CepGenEvent(const cepgen::Event& evt) : GenEvent(Units::GEV, Units::MM) {
  add_attribute("AlphaQCD", make_shared<DoubleAttribute>(evt.metadata("alphaS")));
  add_attribute("AlphaEM", make_shared<DoubleAttribute>(evt.metadata("alphaEM")));

  weights().push_back(1.);  // unweighted events

  // filling the particles content
  const FourVector origin(0., 0., 0., 0.);
  int cm_id = 0;

  auto v1 = make_shared<GenVertex>(origin), v2 = make_shared<GenVertex>(origin), vcm = make_shared<GenVertex>(origin);
  unsigned short idx = 0;
  for (const auto& part_orig : evt.particles()) {
    const auto& mom_orig = part_orig.momentum();
    FourVector momentum(mom_orig.px(), mom_orig.py(), mom_orig.pz(), mom_orig.energy());
    auto part = make_shared<GenParticle>(momentum, part_orig.integerPdgId(), static_cast<int>(part_orig.status()));
    part->set_generated_mass(cepgen::PDG::get().mass(part_orig.pdgId()));
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
          auto production_vertex = assoc_map_.at(m1)->end_vertex();
          std::list<short> ids{m1};  // list of mother particles
          if (assoc_map_.count(m2) != 0 && m2 > m1) {
            ids.resize(m2 - m1 + 1);
            std::iota(ids.begin(), ids.end(), m1);
          }
          if (!production_vertex) {
            production_vertex = make_shared<GenVertex>();
            for (const auto& id : ids)
              production_vertex->add_particle_in(assoc_map_.at(id));
            add_vertex(production_vertex);
          }
          production_vertex->add_particle_out(part);
        } else
          throw CG_FATAL("HepMC3:fillEvent") << "Other particle requested! Not yet implemented!";
      } break;
    }
    idx++;
  }
  add_vertex(v1);
  add_vertex(v2);
  add_vertex(vcm);
}

CepGenEvent::operator cepgen::Event() const {
  cepgen::Event evt;
  auto convert_particle = [](const GenParticle& part,
                             const cepgen::Particle::Role& role =
                                 cepgen::Particle::Role::UnknownRole) -> cepgen::Particle {
    auto convert_momentum = [](const FourVector& mom) -> cepgen::Momentum {
      return cepgen::Momentum::fromPxPyPzE(mom.px(), mom.py(), mom.pz(), mom.e());
    };
    auto cg_part = cepgen::Particle(role, 0, static_cast<cepgen::Particle::Status>(part.status()));
    cg_part.setPdgId((long)part.pdg_id());
    cg_part.setMomentum(convert_momentum(part.momentum()));
    return cg_part;
  };

  auto ip1 = beams().at(0), ip2 = beams().at(1);
  std::unordered_map<size_t, size_t> h_to_cg;
  std::vector<int> beam_vtx_ids;
  for (const auto& vtx : vertices()) {
    if (vtx->particles_in().size() == 1) {
      auto role1 = cepgen::Particle::Role::UnknownRole, role2 = role1, role3 = role1;
      auto status1 = cepgen::Particle::Status::PrimordialIncoming, status2 = cepgen::Particle::Status::Incoming,
           status3 = cepgen::Particle::Status::Unfragmented;
      size_t id_beam_in = 0;
      if (auto& part = vtx->particles_in().at(0); part) {
        if (part->id() == ip1->id()) {
          role1 = cepgen::Particle::Role::IncomingBeam1;
          role2 = cepgen::Particle::Role::Parton1;
          role3 = cepgen::Particle::Role::OutgoingBeam1;
        } else if (part->id() == ip2->id()) {
          role1 = cepgen::Particle::Role::IncomingBeam2;
          role2 = cepgen::Particle::Role::Parton2;
          role3 = cepgen::Particle::Role::OutgoingBeam2;
        }
        auto cg_part = convert_particle(*part, role1);
        cg_part.setStatus(status1);
        evt.addParticle(cg_part);
        h_to_cg[part->id()] = cg_part.id();
        id_beam_in = cg_part.id();
      }
      if (vtx->particles_out_size() == 2) {  //FIXME handle cases with multiple partons?
        size_t num_op = 0;
        for (auto op : vtx->particles_out()) {
          auto cg_part = convert_particle(*op, num_op == 0 ? role2 : role3);
          cg_part.setStatus(num_op == 0 ? status2 : status3);
          cg_part.addMother(evt[id_beam_in]);
          evt.addParticle(cg_part);
          h_to_cg[op->id()] = cg_part.id();
          ++num_op;
        }
      }
      beam_vtx_ids.emplace_back(vtx->id());
    }
  }

  auto cg_intermediate =
      cepgen::Particle(cepgen::Particle::Role::Intermediate, 0, cepgen::Particle::Status::Propagator);
  auto &part1 = evt.oneWithRole(cepgen::Particle::Role::Parton1),
       &part2 = evt.oneWithRole(cepgen::Particle::Role::Parton2);
  cg_intermediate.setMomentum(part1.momentum() + part2.momentum(), true);
  cg_intermediate.addMother(part1);
  cg_intermediate.addMother(part2);
  evt.addParticle(cg_intermediate);

  for (const auto& vtx : vertices()) {
    if (cepgen::utils::contains(beam_vtx_ids, vtx->id()))
      continue;
    for (auto op : vtx->particles_out()) {
      auto cg_part = convert_particle(*op, cepgen::Particle::Role::CentralSystem);
      cg_part.addMother(evt.oneWithRole(cepgen::Particle::Role::Intermediate));
      evt.addParticle(cg_part);
      h_to_cg[op->id()] = cg_part.id();
    }
  }
  return evt;
}

void CepGenEvent::merge(cepgen::Event& evt) const {
  // set of sanity checks to perform on the HepMC event content
  if (vertices().size() < 3) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Failed to retrieve the three primordial vertices in event.";
    return;
  }
  const auto v1 = vertices().at(0), v2 = vertices().at(1), vcm = vertices().at(2);
  if (v1->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particles multiplicity: found "
                                         << v1->particles_in().size() << ", expecting one.";
    return;
  }
  if (v2->particles_in().size() != 1) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particles multiplicity: found "
                                         << v2->particles_in().size() << ", expecting one.";
    return;
  }
  // set of sanity checks to ensure the compatibility between the HepMC and CepGen event records
  const auto ip1 = v1->particles_in().at(0), ip2 = v2->particles_in().at(0);
  const auto &cg_ip1 = evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1),
             &cg_ip2 = evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2);
  if (ip1->momentum().x() != cg_ip1.momentum().px() || ip1->momentum().y() != cg_ip1.momentum().py() ||
      ip1->momentum().z() != cg_ip1.momentum().pz() || ip1->momentum().t() != cg_ip1.momentum().energy()) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particle kinematics.";
    return;
  }
  if (ip2->momentum().x() != cg_ip2.momentum().px() || ip2->momentum().y() != cg_ip2.momentum().py() ||
      ip2->momentum().z() != cg_ip2.momentum().pz() || ip2->momentum().t() != cg_ip2.momentum().energy()) {
    CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particle kinematics.";
    return;
  }
  auto cs = evt[cepgen::Particle::Role::CentralSystem];
  if (cs.size() != (size_t)vcm->particles_out().size()) {
    CG_ERROR("HepMC3:CepGenEvent:merge")
        << "Central system particles multiplicities differ between CepGen and HepMC3 event records.";
    return;
  }
  // freeze the "primordial" central system size
  const auto cs_size = cs.size();

  // helper function to browse particles decay products and store them into the CepGen event content
  std::function<void(const ConstGenParticlePtr& hp, cepgen::ParticleRef cp)> browse_children =
      [&](const ConstGenParticlePtr& hp, cepgen::ParticleRef cp) {
        if (hp->children().empty())
          return;
        cp.get().setStatus(cepgen::Particle::Status::Propagator);
        for (const auto& h_child : hp->children()) {
          cepgen::Particle cg_child(cp.get().role(), 0);
          cg_child.setPdgId((long)h_child->pdg_id());
          const auto& c_mom = h_child->momentum();
          cg_child.setStatus(cepgen::Particle::Status::FinalState);
          cg_child.setMomentum(cepgen::Momentum::fromPxPyPzE(c_mom.x(), c_mom.y(), c_mom.z(), c_mom.t()));
          cg_child.addMother(cp);
          browse_children(h_child, evt.addParticle(cg_child));
        }
      };

  for (size_t icg = 0; icg < cs_size; ++icg) {  // try to find the associated CepGen event particle
    const auto cg_cp_mom3 = cs[icg].get().momentum().p();
    for (const auto& h_cp : vcm->particles_out()) {  // loop over the central system particles
      if (fabs(cg_cp_mom3 - h_cp->momentum().length()) > 1.e-10)
        continue;
      // found the association between the HepMC and CepGen particles kinematics
      browse_children(h_cp, cs[icg]);
      break;
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
      log << "\n  * vertex#" << (-vtxPtr->id()) << " (status: " << vtxPtr->status() << ")"
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
