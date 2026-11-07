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

#include "CepGenPythia8/EventInterface.h"

#include <algorithm>
#include <map>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Value.h"

namespace cepgen::pythia8 {
  /// Convert a CepGen particle momentum into its Pythia 8 counterpart
  Pythia8::Vec4 momToVec4(const Momentum& mom) { return Pythia8::Vec4(mom.px(), mom.py(), mom.pz(), mom.energy()); }
  /// Convert a Pythia 8 particle momentum into its CepGen counterpart
  Momentum vec4ToMom(const Pythia8::Vec4& v4) { return Momentum::fromPxPyPzM(v4.px(), v4.py(), v4.pz(), v4.mCalc()); }

  EventInterface::EventInterface() : Pythia8::LHAup(3), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {}

  void EventInterface::initialise(const RunParameters& params) {
    params_ = &params;
    const auto &pos_beam = params_->kinematics().incomingBeams().positive(),
               &neg_beam = params_->kinematics().incomingBeams().negative();
    inel1_ = !pos_beam.elastic();
    inel2_ = !neg_beam.elastic();
    setBeamA(pos_beam.integerPdgId(), pos_beam.momentum().energy());
    setBeamB(neg_beam.integerPdgId(), neg_beam.momentum().energy());
  }

  void EventInterface::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
    osLHEF << comments;
#else
    CG_WARNING("pythia8:EventInterface:addComments")
        << "Pythia >= 8.2 is required to accept comments on run. Unused comments: " << comments;
#endif
  }

  void EventInterface::setCrossSection(int id, const Value& cross_section) {
    addProcess(0, (double)cross_section, cross_section.uncertainty(), 100.);
    setXSec(id, (double)cross_section);
    setXErr(id, cross_section.uncertainty());
    CG_DEBUG("pythia8:EventInterface:setCrossSection").log([=](auto& log) {
      log << "Cross section set to " << cross_section << " pb. Event kinematics:\n";
      listInit();
    });
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

  bool EventInterface::feedEventForHadronisation(const Event& ev) {
    if (!rndm_)
      throw CG_FATAL("pythia8:EventInterface:feedEventForHadronisation")
          << "No random number engine set for the struck quark flavour selection.";
    const double scale = ev(Particle::Role::Intermediate)[0].momentum().mass();
    setProcess(0, 1., scale, ev.metadata("alphaEM"), ev.metadata("alphaS"));
    cm_mom_ = Momentum();  // everything is fed in the laboratory frame: no (un)boosting

    struct Side {
      Particle::Role in_beam, out_beam, parton;
      bool inelastic;
    };
    const std::array<Side, 2> sides{
        {{Particle::Role::IncomingBeam1, Particle::Role::OutgoingBeam1, Particle::Role::Parton1, inel1_},
         {Particle::Role::IncomingBeam2, Particle::Role::OutgoingBeam2, Particle::Role::Parton2, inel2_}}};
    std::array<Pythia8::Vec4, 2> mom_in, mom_out;
    std::array<int, 2> pdg_in{0, 0}, pdg_out{0, 0}, colour{0, 0};
    std::array<double, 2> x{1., 1.}, mass_in{0., 0.}, mass_out{0., 0.};
    unsigned short colour_index = MIN_COLOUR_INDEX;

    for (size_t i = 0; i < 2; ++i) {
      const auto& side = sides.at(i);
      const auto &beam = ev.oneWithRole(side.in_beam), &parton = ev.oneWithRole(side.parton),
                 &outgoing = ev.oneWithRole(side.out_beam);
      incoming_cg_ids_[i] = beam.id();
      const auto mom_beam = momToVec4(beam.momentum()), mom_parton = momToVec4(parton.momentum());
      if (!side.inelastic) {  // unresolved hadron: the whole beam enters, the intact hadron (p - q) leaves
        dissociated_cg_ids_[i] = -1;
        pdg_in[i] = beam.integerPdgId();
        pdg_out[i] = outgoing.integerPdgId();
        mom_in[i] = mom_beam;
        mom_out[i] = momToVec4(outgoing.momentum());
        mass_in[i] = beam.momentum().mass();
        mass_out[i] = outgoing.momentum().mass();
        continue;
      }
      // dissociated beam: massless quark k = x P collinear to the beam, struck quark k - q left on-shell (massless)
      //   (k - q)^2 = 0  =>  x = q^2 / (2 E_b (q^0 - s q_z)), with s the beam direction
      dissociated_cg_ids_[i] = outgoing.id();
      const double sign = mom_beam.pz() > 0. ? +1. : -1., e_beam = mom_beam.e();
      const double denominator = 2. * e_beam * (mom_parton.e() - sign * mom_parton.pz());
      if (denominator >= 0.)
        return false;
      x[i] = mom_parton.m2Calc() / denominator;
      if (x[i] <= 0. || x[i] >= 1.) {
        CG_DEBUG("pythia8:EventInterface:feedEventForHadronisation")
            << "Invalid struck quark momentum fraction for beam " << i + 1 << ": x=" << x[i] << ".";
        return false;
      }
      pdg_in[i] = pdg_out[i] = pickStruckQuark(beam.integerPdgId());
      colour[i] = colour_index++;
      mom_in[i] = Pythia8::Vec4(0., 0., sign * x[i] * e_beam, x[i] * e_beam);
      mom_out[i] = mom_in[i] - mom_parton;
    }

    setIdX(pdg_in[0], pdg_in[1], x[0], x[1]);
    setPdf(pdg_in[0], pdg_in[1], x[0], x[1], scale, 0., 0., false);

    const auto add_vec4 = [this](
                              int pdg, int status, int moth1, int moth2, int col, const Pythia8::Vec4& mom, double m) {
      const auto lha_id = sizePart();
      addParticle(pdg,
                  status,
                  moth1,
                  moth2,
                  pdg > 0 ? col : 0,
                  pdg > 0 ? 0 : col,
                  mom.px(),
                  mom.py(),
                  mom.pz(),
                  mom.e(),
                  m,
                  0.,
                  9.);
      return lha_id;
    };

    // incoming legs (all of them collinear to their respective beam)
    std::array<unsigned short, 2> lha_in{0, 0};
    for (size_t i = 0; i < 2; ++i) {
      lha_in[i] = add_vec4(pdg_in[i], -1, 0, 0, colour[i], mom_in[i], mass_in[i]);
      if (sides.at(i).inelastic)
        lha_forced_[lha_in[i]] = {sides.at(i).parton, incoming_cg_ids_[i], Particle::Status::Incoming};
      else
        addCorresp(lha_in[i], incoming_cg_ids_[i]);
    }
    // outgoing beam-side legs (carrying the photons' kT and virtuality)
    for (size_t i = 0; i < 2; ++i) {
      const auto lha_out = add_vec4(pdg_out[i], 1, lha_in[0], lha_in[1], colour[i], mom_out[i], mass_out[i]);
      if (sides.at(i).inelastic)
        lha_forced_[lha_out] = {sides.at(i).out_beam, dissociated_cg_ids_[i], Particle::Status::Fragmented};
      else
        addCorresp(lha_out, ev.oneWithRole(sides.at(i).out_beam).id());
    }

    // central system
    const unsigned short central_colour = colour_index++;
    for (const auto& p : ev(Particle::Role::CentralSystem)) {
      range_t colours = {0, 0};
      if (PDG::get().has(p.pdgId()) && PDG::get().colours(p.pdgId()) > 1) {
        if (p.integerPdgId() > 0)  // particle
          colours.first = central_colour;
        else  // anti-particle
          colours.second = central_colour;
      }
      addCepGenParticle(p, p.status() == Particle::Status::Resonance ? 2 : 1, {lha_in[0], lha_in[1]}, colours);
    }
    return true;
  }

  int EventInterface::pickStruckQuark(long beam_pdg_id) const {
    // photon couples to the quark charge: P(q) ~ n_q(valence) * e_q^2
    //   proton (uud): u 8/9, d 1/9 ; neutron (udd): u 2/3, d 1/3
    const auto abs_pdg = std::labs(beam_pdg_id);
    double prob_up = 8. / 9.;
    if (abs_pdg == 2112)
      prob_up = 2. / 3.;
    else if (abs_pdg != 2212)
      CG_WARNING("pythia8:EventInterface:pickStruckQuark")
          << "Beam with PDG id " << beam_pdg_id << " is not a nucleon. Assuming a proton valence content.";
    const int quark = rndm_->flat() < prob_up ? 2 : 1;
    return beam_pdg_id > 0 ? quark : -quark;
  }

  void EventInterface::updateEvent(const Pythia8::Event& pyevt, Event& ev, double& weight) const {
    // Pythia event record: 0 = two-beam system, 1/2 = incoming beams, then the LHA entries (k -> k + 2),
    // then everything produced by Pythia (resonance decays, beam remnants, showers, hadrons)
    std::map<int, int> py_cg;  // Pythia index -> CepGen id
    for (const auto& [lha_id, cg_id] : lha_cg_corresp_)
      py_cg[lha_id + 2] = cg_id;
    if (pyevt.size() < sizePart() + 2) {
      CG_WARNING("pythia8:EventInterface:updateEvent")
          << "Failed to update the event with (possibly invalid) Pythia output.";
      return;
    }
    const auto is_beam_side = [](const Particle::Role& role) {
      return role == Particle::Role::OutgoingBeam1 || role == Particle::Role::OutgoingBeam2;
    };
    for (int i = 3; i < pyevt.size(); ++i) {
      const auto& pypart = pyevt.at(i);
      if (auto it = py_cg.find(i); it != py_cg.end()) {  // particle already in the CepGen event
        auto& cgpart = ev[it->second];
        if (cgpart.role() == Particle::Role::CentralSystem && pypart.status() < 0) {
          if (pypart.isResonance()) {  // resonance decayed by Pythia; apply branching ratio for this decay
            weight *= pypart.particleDataEntry().pickChannel().bRatio();
            cgpart.setStatus(Particle::Status::Resonance);
          } else
            cgpart.setStatus(Particle::Status::Propagator);  // e.g. recoiler copied by the shower
        }
        continue;
      }
      checkPDGid(pypart);
      Particle::Role role{Particle::Role::UnknownRole};
      std::vector<int> cg_mothers;
      auto non_final_status = Particle::Status::Propagator;
      if (const auto it_forced = lha_forced_.find(i - 2); i - 2 < sizePart() && it_forced != lha_forced_.end()) {
        // incoming or struck quark of a dissociated beam
        role = it_forced->second.role;
        cg_mothers.emplace_back(it_forced->second.cg_mother);
        non_final_status = it_forced->second.status;
      } else if (const auto moth1 = pypart.mother1(); moth1 == 1 || moth1 == 2) {
        // directly attached to a Pythia beam: initiator copy or beam remnant
        const size_t side = moth1 - 1;
        if (pypart.status() == -61 || pypart.status() == -21) {
          role = side == 0 ? Particle::Role::Parton1 : Particle::Role::Parton2;
          cg_mothers.emplace_back(incoming_cg_ids_.at(side));
          non_final_status = Particle::Status::Incoming;
        } else {
          role = side == 0 ? Particle::Role::OutgoingBeam1 : Particle::Role::OutgoingBeam2;
          cg_mothers.emplace_back(dissociated_cg_ids_.at(side) >= 0 ? dissociated_cg_ids_.at(side)
                                                                    : incoming_cg_ids_.at(side));
          non_final_status = Particle::Status::Fragmented;
        }
      } else {  // child inherits its (first) mother's role
        for (const auto& moth : {pypart.mother1(), pypart.mother2()})
          if (moth > 0 && py_cg.count(moth) > 0 &&
              std::find(cg_mothers.begin(), cg_mothers.end(), py_cg.at(moth)) == cg_mothers.end())
            cg_mothers.emplace_back(py_cg.at(moth));
        if (!cg_mothers.empty())
          role = ev[cg_mothers.front()].role();
        if (is_beam_side(role))
          non_final_status = Particle::Status::Fragmented;
      }
      auto& cgpart = ev.addParticle(role).get();
      cgpart.setIntegerPdgId(pypart.id());
      cgpart.setMomentum(vec4ToMom(pypart.p()), true);
      cgpart.setStatus(pypart.isFinal()       ? Particle::Status::FinalState
                       : pypart.isResonance() ? Particle::Status::Resonance
                                              : non_final_status);
      const auto cg_id = cgpart.id();
      py_cg[i] = cg_id;
      for (const auto& moth_id : cg_mothers)
        ev[cg_id].addMother(ev[moth_id]);
    }
    for (const auto& diff_id : dissociated_cg_ids_)  // outgoing diffractive systems are now fragmented
      if (diff_id >= 0 && !ev[diff_id].daughters().empty())
        ev[diff_id].setStatus(Particle::Status::Fragmented);
  }

  void EventInterface::setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
    LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
    lha_cg_corresp_.clear();
    lha_forced_.clear();
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
      status = 1;
      if (const auto& ps = part.status(); ps == Particle::Status::Resonance || ps == Particle::Status::Fragmented)
        status = 2;
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
    checkPDGid(part.particleDataEntry());
    auto& prop = PDG::get()[part.idAbs()];
    prop.colours = part.col();  // colour factor
  }

  void EventInterface::checkPDGid(const Pythia8::ParticleDataEntry& data) {
    if (cepgen::PDG::get().has(std::abs(data.id())))
      return;
    cepgen::ParticleProperties prop;
    prop.pdgid = std::abs(data.id());
    prop.name = prop.human_name = data.name();
    prop.colours = data.colType();  // colour factor
    prop.mass = data.m0();
    prop.width = data.mWidth();
    if (const short ch = data.charge() * 3; std::abs(ch) > 0)
      prop.charges = {ch, -ch};
    prop.fermion = data.isLepton();
    PDG::get().define(prop);
  }
}  // namespace cepgen::pythia8
