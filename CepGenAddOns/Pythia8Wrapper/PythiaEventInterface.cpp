#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

namespace Pythia8 {
  /// Convert a CepGen particle momentum into its Pythia8 counterpart
  Vec4 momToVec4(const cepgen::Momentum& mom) { return Vec4(mom.px(), mom.py(), mom.pz(), mom.energy()); }

  CepGenEvent::CepGenEvent()
      : LHAup(3),
        mp_(cepgen::PDG::get().mass(cepgen::PDG::proton)),
        mp2_(mp_ * mp_),
        inel1_(false),
        inel2_(false),
        params_(nullptr) {}

  void CepGenEvent::initialise(const cepgen::Parameters& params) {
    params_ = &params;
    inel1_ = params_->kinematics.incomingBeams().positive().mode == cepgen::mode::Beam::ProtonInelastic;
    inel2_ = params_->kinematics.incomingBeams().negative().mode == cepgen::mode::Beam::ProtonInelastic;

    setBeamA((short)params_->kinematics.incomingBeams().positive().pdg,
             params_->kinematics.incomingBeams().positive().momentum.pz());
    setBeamB((short)params_->kinematics.incomingBeams().negative().pdg,
             params_->kinematics.incomingBeams().negative().momentum.pz());
    //addProcess( 0, params_->integration().result, params_->integration().err_result, 100. );
  }

  void CepGenEvent::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
    osLHEF << comments;
#endif
  }

  void CepGenEvent::setCrossSection(int id, double cross_section, double cross_section_err) {
    addProcess(0, cross_section, cross_section_err, 100.);
    setXSec(id, cross_section);
    setXErr(id, cross_section_err);
    //listInit();
  }

  void CepGenEvent::feedEvent(const cepgen::Event& ev, const Type& type) {
    const double scale = ev[cepgen::Particle::Intermediate][0].mass();
    setProcess(0, 1., scale, cepgen::constants::ALPHA_EM, cepgen::constants::ALPHA_QCD);

    const auto &part1 = ev[cepgen::Particle::Parton1][0], &part2 = ev[cepgen::Particle::Parton2][0];
    const auto &op1 = ev[cepgen::Particle::OutgoingBeam1][0], &op2 = ev[cepgen::Particle::OutgoingBeam2][0];
    const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
    const double x1 = q2_1 / (q2_1 + op1.mass2() - mp2_), x2 = q2_2 / (q2_2 + op2.mass2() - mp2_);

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
        mom_iq1 = momToVec4(x1 * ev[cepgen::Particle::IncomingBeam1][0].momentum());
      }
      if (inel2_) {
        parton2_pdgid = 2;
        parton2_colour = colour_index++;
        mom_iq2 = momToVec4(x2 * ev[cepgen::Particle::IncomingBeam2][0].momentum());
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

      if (type == Type::centralAndFullBeamRemnants) {
        //=========================================================================================
        // full beam remnants content
        //=========================================================================================

        for (const auto& syst : {cepgen::Particle::OutgoingBeam1, cepgen::Particle::OutgoingBeam2}) {
          for (const auto& p : ev[syst])
            addCepGenParticle(p, INVALID_ID, findMothers(ev, p));
        }
      }
    }

    //=============================================================================================
    // central system
    //=============================================================================================

    const unsigned short central_colour = colour_index++;
    for (const auto& p : ev[cepgen::Particle::CentralSystem]) {
      std::pair<int, int> colours = {0, 0}, mothers = {1, 2};
      if (type != Type::centralAndBeamRemnants)
        mothers = findMothers(ev, p);
      try {
        if (cepgen::PDG::get().colours(p.pdgId()) > 1) {
          if (p.integerPdgId() > 0)  //--- particle
            colours.first = central_colour;
          else  //--- anti-particle
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

  unsigned short CepGenEvent::cepgenId(unsigned short py_id) const {
    if (py_cg_corresp_.count(py_id) == 0)
      return INVALID_ID;
    return py_cg_corresp_.at(py_id);
  }

  unsigned short CepGenEvent::pythiaId(unsigned short cg_id) const {
    for (const auto& py_cg : py_cg_corresp_)
      if (py_cg.second == cg_id)
        return py_cg.first;
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
          if (part.pdgId() == 21 && (int)part.status() == 12)
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

  void CepGenEvent::addCorresp(unsigned short py_id, unsigned short cg_id) { py_cg_corresp_[py_id] = cg_id; }

  void CepGenEvent::dumpCorresp() const {
    CG_INFO("CepGenEvent:dump").log([&](auto& msg) {
      msg << "List of Pythia ←|→ CepGen particle ids correspondence";
      for (const auto& py_cg : py_cg_corresp_)
        msg << "\n\t" << py_cg.first << " <-> " << py_cg.second;
    });
  }

  std::pair<int, int> CepGenEvent::findMothers(const cepgen::Event& ev, const cepgen::Particle& p) const {
    std::pair<int, int> out = {0, 0};

    const auto& mothers = p.mothers();
    if (mothers.empty())
      return out;
    const unsigned short moth1_cg_id = *mothers.begin();
    out.first = pythiaId(moth1_cg_id);
    if (out.first == INVALID_ID) {
      const auto& moth = ev[moth1_cg_id];
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
}  // namespace Pythia8
