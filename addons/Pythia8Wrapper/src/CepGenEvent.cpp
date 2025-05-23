/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2025  Laurent Forthomme
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
#include "CepGenPythia8/CepGenEvent.h"

using namespace Pythia8;

/// Convert a CepGen particle momentum into its Pythia8 counterpart
Vec4 momToVec4(const cepgen::Momentum& mom) { return {mom.px(), mom.py(), mom.pz(), mom.energy()}; }

CepGenEvent::CepGenEvent() : LHAup(3), mp_(cepgen::PDG::get().mass(cepgen::PDG::proton)), mp2_(mp_ * mp_) {}

void CepGenEvent::initialise(const cepgen::RunParameters& params) {
  params_ = &params;
  inel1_ = !params_->kinematics().incomingBeams().positive().elastic();
  inel2_ = !params_->kinematics().incomingBeams().negative().elastic();

  setBeamA(params_->kinematics().incomingBeams().positive().integerPdgId(),
           params_->kinematics().incomingBeams().positive().momentum().pz());
  setBeamB(params_->kinematics().incomingBeams().negative().integerPdgId(),
           params_->kinematics().incomingBeams().negative().momentum().pz());
  //addProcess( 0, params_->integration().result, params_->integration().err_result, 100. );
}

void CepGenEvent::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
  osLHEF << comments;
#else
  CG_WARNING("CepGenEvent:addComments") << "Pythia 8 is too outdated... Unused comments: " << comments;
#endif
}

void CepGenEvent::setCrossSection(int id, double cross_section, double cross_section_err) {
  addProcess(0, cross_section, cross_section_err, 100.);
  setXSec(id, cross_section);
  setXErr(id, cross_section_err);
  //listInit();
}

void CepGenEvent::feedEvent(const cepgen::Event& ev, const Type& type) {
  const double scale = ev(cepgen::Particle::Role::Intermediate)[0].momentum().mass();
  setProcess(0, 1., scale, ev.metadata("alphaEM"), ev.metadata("alphaS"));

  const auto &part1 = ev(cepgen::Particle::Role::Parton1)[0], &part2 = ev(cepgen::Particle::Role::Parton2)[0];
  const auto &op1 = ev(cepgen::Particle::Role::OutgoingBeam1)[0], &op2 = ev(cepgen::Particle::Role::OutgoingBeam2)[0];
  const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
  const double x1 = q2_1 / (q2_1 + op1.momentum().mass2() - mp2_), x2 = q2_2 / (q2_2 + op2.momentum().mass2() - mp2_);

  unsigned short colour_index = MIN_COLOUR_INDEX;

  const Vec4 mom_part1(momToVec4(part1.momentum())), mom_part2(momToVec4(part2.momentum()));

  if (type == Type::centralAndBeamRemnants) {  // full event content (with collinear partons)
    Vec4 mom_iq1 = mom_part1, mom_iq2 = mom_part2;
    unsigned short parton1_id, parton2_id;
    unsigned short parton1_pdgid = part1.integerPdgId(), parton2_pdgid = part2.integerPdgId();
    unsigned short parton1_colour = 0, parton2_colour = 0;
    //FIXME select quark flavours accordingly
    if (inel1_) {
      parton1_pdgid = 2;
      parton1_colour = colour_index++;
      mom_iq1 = momToVec4(x1 * ev(cepgen::Particle::Role::IncomingBeam1)[0].momentum());
    }
    if (inel2_) {
      parton2_pdgid = 2;
      parton2_colour = colour_index++;
      mom_iq2 = momToVec4(x2 * ev(cepgen::Particle::Role::IncomingBeam2)[0].momentum());
    }

    //--- flavour / x value of hard-process initiators
    setIdX(part1.integerPdgId(), part2.integerPdgId(), x1, x2);
    setPdf(parton1_pdgid, parton2_pdgid, x1, x2, scale, 0., 0., false);

    //===========================================================================================
    // incoming valence quarks
    //===========================================================================================

    parton1_id = sizePart();
    addCorresp(parton1_id, op1.id());
    addParticle(parton1_pdgid,
                -1,
                0,
                0,
                parton1_colour,
                0,
                mom_iq1.px(),
                mom_iq1.py(),
                mom_iq1.pz(),
                mom_iq1.e(),
                mom_iq1.mCalc(),
                0.,
                1.);

    parton2_id = sizePart();
    addCorresp(parton2_id, op2.id());
    addParticle(parton2_pdgid,
                -1,
                0,
                0,
                parton2_colour,
                0,
                mom_iq2.px(),
                mom_iq2.py(),
                mom_iq2.pz(),
                mom_iq2.e(),
                mom_iq2.mCalc(),
                0.,
                1.);

    //===========================================================================================
    // outgoing valence quarks
    //===========================================================================================

    if (inel1_) {
      const Vec4 mom_oq1 = mom_iq1 - mom_part1;
      addParticle(parton1_pdgid,
                  1,
                  parton1_id,
                  parton2_id,
                  parton1_colour,
                  0,
                  mom_oq1.px(),
                  mom_oq1.py(),
                  mom_oq1.pz(),
                  mom_oq1.e(),
                  mom_oq1.mCalc(),
                  0.,
                  1.);
    }
    if (inel2_) {
      const Vec4 mom_oq2 = mom_iq2 - mom_part2;
      addParticle(parton2_pdgid,
                  1,
                  parton1_id,
                  parton2_id,
                  parton2_colour,
                  0,
                  mom_oq2.px(),
                  mom_oq2.py(),
                  mom_oq2.pz(),
                  mom_oq2.e(),
                  mom_oq2.mCalc(),
                  0.,
                  1.);
    }
  } else {
    //===========================================================================================
    // incoming partons
    //===========================================================================================

    addCepGenParticle(part1, -2);
    addCepGenParticle(part2, -2);

    if (type == Type::centralAndFullBeamRemnants)  // full beam remnants content
      for (const auto& forward_system : {cepgen::Particle::Role::OutgoingBeam1, cepgen::Particle::Role::OutgoingBeam2})
        for (const auto& p : ev(forward_system))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));
  }

  //=============================================================================================
  // central system
  //=============================================================================================

  const unsigned short central_colour = colour_index++;
  for (const auto& p : ev(cepgen::Particle::Role::CentralSystem)) {
    std::pair<int, int> colours = {0, 0}, mothers = {1, 2};
    if (type != Type::centralAndBeamRemnants)
      mothers = findMothers(ev, p);
    try {
      if (cepgen::PDG::get().colours(p.pdgId()) > 1) {
        if (p.integerPdgId() > 0)  // particle
          colours.first = central_colour;
        else  // anti-particle
          colours.second = central_colour;
      }
    } catch (const cepgen::Exception&) {
    }
    int status = 1;
    if (type == Type::centralAndFullBeamRemnants && p.status() == cepgen::Particle::Status::Resonance)
      status = 2;
    addCepGenParticle(p, status, mothers, colours);
  }
}

void CepGenEvent::setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
  LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
  py_cg_corresp_.clear();
}

unsigned short CepGenEvent::cepgenId(unsigned short pythia_id) const {
  if (py_cg_corresp_.count(pythia_id) == 0)
    return INVALID_ID;
  return py_cg_corresp_.at(pythia_id);
}

unsigned short CepGenEvent::pythiaId(unsigned short cepgen_id) const {
  const auto it = std::find_if(py_cg_corresp_.begin(), py_cg_corresp_.end(), [&cepgen_id](const auto& py_cg) {
    return py_cg.second == cepgen_id;
  });
  if (it != py_cg_corresp_.end())
    return it->first;
  return INVALID_ID;
}

void CepGenEvent::addCepGenParticle(const cepgen::Particle& part,
                                    int status,
                                    const std::pair<int, int>& mothers,
                                    const std::pair<int, int>& colours) {
  const Vec4 mom_part(momToVec4(part.momentum()));
  int pdg_id = part.integerPdgId();
  if (status == INVALID_ID)
    switch (part.status()) {
      case cepgen::Particle::Status::Resonance:
      case cepgen::Particle::Status::Fragmented:
        status = 2;
        break;
      default: {
        if (part.pdgId() == 21 && static_cast<int>(part.status()) == 12)
          pdg_id = -21;  // workaround for HepMC2 interface
        else
          status = 1;
      } break;
    }
  addCorresp(sizePart(), part.id());
  addParticle(pdg_id,
              status,
              mothers.first,
              mothers.second,
              colours.first,
              colours.second,
              mom_part.px(),
              mom_part.py(),
              mom_part.pz(),
              mom_part.e(),
              mom_part.mCalc(),
              0.,
              0.,
              0.);
}

void CepGenEvent::addCorresp(unsigned short pythia_id, unsigned short cepgen_id) {
  py_cg_corresp_[pythia_id] = cepgen_id;
}

void CepGenEvent::dumpCorresp() const {
  CG_INFO("CepGenEvent:dump").log([&](auto& msg) {
    msg << "List of Pythia ←|→ CepGen particle ids correspondence";
    for (const auto& py_cg : py_cg_corresp_)
      msg << "\n\t" << py_cg.first << " <-> " << py_cg.second;
  });
}

std::pair<int, int> CepGenEvent::findMothers(const cepgen::Event& cepgen_event,
                                             const cepgen::Particle& cepgen_particle) const {
  std::pair out = {0, 0};

  const auto& mothers = cepgen_particle.mothers();
  if (mothers.empty())
    return out;
  const unsigned short moth1_cg_id = *mothers.begin();
  if (out.first = pythiaId(moth1_cg_id); out.first == INVALID_ID) {
    const auto& moth = cepgen_event(moth1_cg_id);
    out = {(moth.mothers().size() > 0) ? pythiaId(*moth.mothers().begin()) : 0,
           (moth.mothers().size() > 1) ? pythiaId(*moth.mothers().rbegin()) : 0};
  }
  if (mothers.size() > 1) {
    const unsigned short moth2_cg_id = *mothers.rbegin();
    out.second = pythiaId(moth2_cg_id);
    if (out.second == INVALID_ID)
      out.second = 0;
  }
  return out;
}
