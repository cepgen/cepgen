/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include <HepMC3/Version.h>

#include <list>
#include <numeric>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenAddOns/HepMC3Wrapper/HepMC3EventInterface.h"

namespace HepMC3 {
  CepGenEvent::CepGenEvent(const cepgen::Event& evt) : GenEvent(Units::GEV, Units::MM) {
    add_attribute("AlphaQCD", make_shared<DoubleAttribute>(cepgen::constants::ALPHA_QCD));
    add_attribute("AlphaEM", make_shared<DoubleAttribute>(cepgen::constants::ALPHA_EM));

    weights().push_back(1.);  // unweighted events

    // filling the particles content
    const FourVector origin(0., 0., 0., 0.);
    int cm_id = 0;

    auto v1 = make_shared<GenVertex>(origin), v2 = make_shared<GenVertex>(origin), vcm = make_shared<GenVertex>(origin);
    unsigned short idx = 0;
    for (const auto& part_orig : evt.particles()) {
      const auto& mom_orig = part_orig.momentum();
      FourVector pmom(mom_orig.px(), mom_orig.py(), mom_orig.pz(), part_orig.energy());
      auto part = make_shared<GenParticle>(pmom, part_orig.pdgId(), (int)part_orig.status());
      part->set_generated_mass(cepgen::PDG::get().mass(part_orig.pdgId()));
      assoc_map_[idx] = part;

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
              vprod = make_shared<GenVertex>();
              for (const auto& id : ids)
                vprod->add_particle_in(assoc_map_.at(id));
              add_vertex(vprod);
            }
            vprod->add_particle_out(part);
          } else
            throw CG_FATAL("HepMCHandler:fillEvent") << "Other particle requested! Not yet implemented!";
        } break;
      }
      idx++;
    }
    add_vertex(v1);
    add_vertex(v2);
    add_vertex(vcm);
  }

  std::ostream& operator<<(std::ostream& os, const FourVector& vec) {
    return os << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << "; " << vec.t() << ")";
  }

  void CepGenEvent::merge(cepgen::Event& evt) const {
    // set of sanity checks to perform on the HepMC event content
    if (vertices_size() < 3) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Failed to retrieve the three primordial vertices in event.";
      return;
    }
    const auto v1 = vertices().at(0), v2 = vertices().at(1), vcm = vertices().at(2);
    if (v1->particles_in_size() != 1) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid first incoming beam particles multiplicity: found "
                                           << v1->particles_in_size() << ", expecting one.";
      return;
    }
    if (v2->particles_in_size() != 1) {
      CG_ERROR("HepMC3:CepGenEvent:merge") << "Invalid second incoming beam particles multiplicity: found "
                                           << v2->particles_in_size() << ", expecting one.";
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
    if (cs.size() != (size_t)vcm->particles_out_size()) {
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
          log << "\n      * " << ipPtr->pdg_id() << ": " << ipPtr->momentum(), in_sys += ipPtr->momentum();
        log << "\n     total: " << in_sys << "\n     out:";
        for (const auto& opPtr : vtxPtr->particles_out())
          log << "\n      * " << opPtr->pdg_id() << ": " << opPtr->momentum(), out_sys += opPtr->momentum();
        const auto imbal = in_sys - out_sys;
        log << "\n     total: " << out_sys << "\n    (im)balance: " << imbal << " (norm: " << imbal.length() << ").";
      }
      log << "\n" << std::string(70, '-');
    });
  }
}  // namespace HepMC3
