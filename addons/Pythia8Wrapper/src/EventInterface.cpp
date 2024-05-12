/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenPythia8/EventInterface.h"

namespace cepgen::pythia8 {
  /// Convert a CepGen particle momentum into its Pythia 8 counterpart
  Pythia8::Vec4 momToVec4(const Momentum& mom) { return Pythia8::Vec4(mom.px(), mom.py(), mom.pz(), mom.energy()); }
  /// Convert a Pythia 8 particle momentum into its CepGen counterpart
  Momentum vec4ToMom(const Pythia8::Vec4& v4) { return Momentum::fromPxPyPzM(v4.px(), v4.py(), v4.pz(), v4.mCalc()); }

  EventInterface::EventInterface() : Pythia8::LHAup(3), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {}

  void EventInterface::initialise(const RunParameters& params) {
    params_ = &params;
    inel1_ = !params_->kinematics().incomingBeams().positive().elastic();
    inel2_ = !params_->kinematics().incomingBeams().negative().elastic();

    setBeamA((short)params_->kinematics().incomingBeams().positive().integerPdgId(),
             params_->kinematics().incomingBeams().positive().momentum().energy());
    setBeamB((short)params_->kinematics().incomingBeams().negative().integerPdgId(),
             params_->kinematics().incomingBeams().negative().momentum().energy());
    //addProcess( 0, params_->integration().result, params_->integration().err_result, 100. );
  }

  void EventInterface::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
    osLHEF << comments;
#else
    CG_WARNING("pythia8:EventInterface:addComments") << "Pythia 8 is too outdated... Unused comments: " << comments;
#endif
  }

  void EventInterface::setCrossSection(int id, double cross_section, double cross_section_err) {
    addProcess(0, cross_section, cross_section_err, 100.);
    setXSec(id, cross_section);
    setXErr(id, cross_section_err);
    //listInit();
  }

  void EventInterface::feedEvent(const Event& ev) {
    const double scale = ev(Particle::Role::Intermediate)[0].momentum().mass();
    setProcess(0, 1., ev.cmEnergy(), ev.metadata("alphaEM"), ev.metadata("alphaS"));

    const auto &part1 = ev(Particle::Role::Parton1)[0], &part2 = ev(Particle::Role::Parton2)[0];
    const auto &op1 = ev(Particle::Role::OutgoingBeam1)[0], &op2 = ev(Particle::Role::OutgoingBeam2)[0];
    const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
    const double x1 = q2_1 / (q2_1 + op1.momentum().mass2() - mp2_), x2 = q2_2 / (q2_2 + op2.momentum().mass2() - mp2_);

    unsigned short colour_index = MIN_COLOUR_INDEX;

    // incoming partons
    setPdf(part1.integerPdgId(), part2.integerPdgId(), x1, x2, scale, 0., 0., false);
    cm_mom_ = (part1.momentum() + part2.momentum()).transverse();
    const auto parton1_id = addCepGenParticle(part1, -2), parton2_id = addCepGenParticle(part2, -2);

    if (store_remnants_)  // full beam remnants content
      for (const auto& syst : {Particle::Role::OutgoingBeam1, Particle::Role::OutgoingBeam2})
        for (const auto& p : ev(syst))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));

    // central system
    const unsigned short central_colour = colour_index++;
    for (const auto& p : ev(Particle::Role::CentralSystem)) {
      range_t colours = {0, 0}, mothers = {1, 2};
      if (mothers == std::make_pair((int)INVALID_ID, (int)INVALID_ID))
        mothers = std::make_pair(parton1_id, parton2_id);
      if (PDG::get().has(p.pdgId()) && PDG::get().colours(p.pdgId()) > 1) {
        if (p.integerPdgId() > 0)  // particle
          colours.first = central_colour;
        else  // anti-particle
          colours.second = central_colour;
      }
      int status = 1;
      if (p.status() == Particle::Status::Resonance)
        status = 2;
      addCepGenParticle(p, status, mothers, colours);
    }
  }

  void EventInterface::updateEvent(const Pythia8::Event& pyevt, Event& ev, double& weight) const {
    pyevt.list();
    std::map<unsigned short, unsigned short> pyid_vs_cgid{{1, 5}, {2, 6}};  // keep it ordered...
    for (auto& py_cg : lha_cg_corresp_)
      pyid_vs_cgid[py_cg.first + 2 /* Pythia adds the two incoming beam particles to event content */] = py_cg.second;
    if (pyevt.size() <= (int)lha_cg_corresp_.size() + 3) {
      CG_WARNING("pythia6:EventInterface:updateEvent")
          << "Failed to update the event with (possibly invalid) Pythia output.";
      return;
    }
    // 0 = two-beam system
    // 1 = incoming beam 1
    // 2 = incoming beam 2
    for (int i = lha_cg_corresp_.size() + 3; i < pyevt.size(); ++i) {  // 1st loop to add particles contents
      const auto& pypart = pyevt.at(i);
      checkPDGid(pypart);
      auto& cgpart = ev.addParticle(Particle::Role::Intermediate).get();
      cgpart.setStatus(Particle::Status::DebugResonance);
      cgpart.setIntegerPdgId(pypart.id());
      cgpart.setMomentum(vec4ToMom(pypart.p()).lorentzBoost(cm_mom_));
      pyid_vs_cgid[i] = cgpart.id();
    }
    for (const auto& py_cg : pyid_vs_cgid) {  // 2nd loop to establish parentage
      const auto& pypart = pyevt.at(py_cg.first);
      auto& cgpart = ev[py_cg.second];
      if (pyid_vs_cgid.count(pypart.mother1()) > 0)
        cgpart.addMother(ev[pyid_vs_cgid.at(pypart.mother1())]);
      if (pyid_vs_cgid.count(pypart.mother2()) > 0)
        cgpart.addMother(ev[pyid_vs_cgid.at(pypart.mother2())]);
      if (cgpart.role() == Particle::Role::Intermediate) {  // invalid role ; need to update from parentage
        if (const auto& moths = cgpart.mothers(); !moths.empty()) {
          const auto moth_role = ev[*moths.begin()].role();  // we only account for the first mother
          if (pypart.status() == -61) {                      // intermediate partons
            if (moth_role == Particle::Role::OutgoingBeam1) {
              cgpart.setRole(Particle::Role::Parton1);
              ev.clearMothers(cgpart);  // patch to set incoming beam as only mother
              cgpart.addMother(ev[Particle::Role::IncomingBeam1][0]);
            } else if (moth_role == Particle::Role::OutgoingBeam2) {
              cgpart.setRole(Particle::Role::Parton2);
              ev.clearMothers(cgpart);  // patch to set incoming beam as only mother
              cgpart.addMother(ev[Particle::Role::IncomingBeam2][0]);
            }
            cgpart.setStatus(Particle::Status::Incoming);
          } else
            cgpart.setRole(moth_role);  // child inherits its parent's role
        }
      }
      if (cgpart.status() == Particle::Status::DebugResonance) {  // fix whatever status we can fix
        if (pypart.isResonance()) {
          if (cgpart.role() == Particle::Role::CentralSystem && pypart.status() < 0)
            weight *= pypart.particleDataEntry().pickChannel().bRatio();
          cgpart.setStatus(Particle::Status::Resonance);
        } else if (pypart.status() > 0)
          cgpart.setStatus(Particle::Status::FinalState);
      }
    }
    ev.updateRoles();  // update all newly-reassigned roles after 2nd loop
    // post-fix to set outgoing diffractive systems as fragmented
    if (inel1_)
      if (auto& diffx = ev[Particle::Role::OutgoingBeam1][0].get(); !diffx.daughters().empty())
        diffx.setStatus(Particle::Status::Fragmented);
    if (inel2_)
      if (auto& diffy = ev[Particle::Role::OutgoingBeam2][0].get(); !diffy.daughters().empty())
        diffy.setStatus(Particle::Status::Fragmented);
  }

  void EventInterface::setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
    LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
    lha_cg_corresp_.clear();
  }

  unsigned short EventInterface::lhaId(unsigned short cg_id) const {
    if (auto it = std::find_if(lha_cg_corresp_.begin(),
                               lha_cg_corresp_.end(),
                               [&cg_id](const auto& py_cg) { return py_cg.second == cg_id; });
        it != lha_cg_corresp_.end())
      return it->first;
    return INVALID_ID;
  }

  unsigned short EventInterface::addCepGenParticle(const Particle& part,
                                                   int status,
                                                   const EventInterface::range_t& mothers,
                                                   const EventInterface::range_t& colours) {
    if (status == INVALID_ID) {
      if (const auto& ps = part.status(); ps == Particle::Status::Resonance || ps == Particle::Status::Fragmented)
        status = 2;
      else
        status = 1;
    }
    const auto py_id = sizePart();
    addCorresp(py_id, part.id());
    const auto mom = Momentum(part.momentum()).lorentzBoost(-cm_mom_);
    addParticle(part.integerPdgId(),
                status,
                mothers.first,
                mothers.second,
                colours.first,
                colours.second,
                mom.px(),
                mom.py(),
                mom.pz(),
                mom.energy(),
                mom.mass(),
                0.,
                0.,
                0.);
    return py_id;
  }

  void EventInterface::addCorresp(unsigned short py_id, unsigned short cg_id) { lha_cg_corresp_[py_id] = cg_id; }

  void EventInterface::dumpCorresp() const {
    CG_INFO("pythia8:EventInterface:dump").log([&](auto& msg) {
      msg << "List of Pythia ←|→ CepGen particle ids correspondence";
      for (const auto& py_cg : lha_cg_corresp_)
        msg << "\n\t" << py_cg.first << " <-> " << py_cg.second;
    });
  }

  EventInterface::range_t EventInterface::findMothers(const Event& ev, const Particle& p) const {
    range_t out{0, 0};
    const auto& mothers = p.mothers();
    if (mothers.empty())
      return out;
    const unsigned short moth1_cg_id = *mothers.begin();
    if (out.first = lhaId(moth1_cg_id); out.first == INVALID_ID) {  // did not find the Pythia equivalent to mother
      const auto& moth = ev(moth1_cg_id);
      out = {moth.mothers().size() > 0 ? lhaId(*moth.mothers().begin()) : 0,
             moth.mothers().size() > 1 ? lhaId(*moth.mothers().rbegin()) : 0};
    }
    if (mothers.size() > 1)
      if (out.second = lhaId(*mothers.rbegin()); out.second == INVALID_ID)
        out.second = 0;
    return out;
  }

  void EventInterface::checkPDGid(const Pythia8::Particle& part) {
    if (cepgen::PDG::get().has(part.idAbs()))
      return;
    cepgen::ParticleProperties prop;
    prop.pdgid = part.idAbs();
    prop.name = prop.human_name = part.name();
    prop.colours = part.col();  // colour factor
    prop.mass = part.m0();
    prop.width = part.mWidth();
    if (const auto ch = int(part.charge() * 3.); std::abs(ch) > 0)
      prop.charges = {ch, -ch};
    prop.fermion = part.isLepton();
    PDG::get().define(prop);
  }
}  // namespace cepgen::pythia8
