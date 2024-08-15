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

#include <unordered_map>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Value.h"
#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

namespace cepgen::hadr {
  /// Interface to the Pythia8 hadronisation algorithm
  /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
  class Pythia8Hadroniser : public Hadroniser {
  public:
    explicit Pythia8Hadroniser(const ParametersList& plist)
        : Hadroniser(plist),
          pythia_(new Pythia8::Pythia),
          cg_evt_(new Pythia8::CepGenEvent),
          correct_central_(steer<bool>("correctCentralSystem")),
          debug_lhef_(steer<bool>("debugLHEF")),
          output_config_(steer<std::string>("outputConfig")) {}

    virtual ~Pythia8Hadroniser() {
      if (!output_config_.empty())
        pythia_->settings.writeFile(output_config_, false);
      if (debug_lhef_)
        cg_evt_->closeLHEF(true);
    }

    static ParametersDescription description();

    void readString(const std::string& param) override {
      if (!pythia_->readString(param))
        throw CG_FATAL("Pythia8Hadroniser") << "The Pythia8 core failed to parse the following setting:\n\t" << param;
    }
    void initialise() override;
    bool run(Event& ev, double& weight, bool fast) override;

    void setCrossSection(const Value& cross_section) override {
      cg_evt_->setCrossSection(0, cross_section, cross_section.uncertainty());
    }

  private:
    void* enginePtr() override { return (void*)pythia_.get(); }

    static constexpr unsigned short PYTHIA_STATUS_IN_BEAM = 12;
    static constexpr unsigned short PYTHIA_STATUS_IN_PARTON_KT = 61;

    pdgids_t min_ids_;
    std::unordered_map<short, short> py_cg_corresp_;
    unsigned short findRole(const Event& ev, const Pythia8::Particle& p) const;
    void updateEvent(Event& ev, double& weight) const;
    Particle& addParticle(Event& ev, const Pythia8::Particle&, const Pythia8::Vec4& mom, unsigned short) const;

    const std::unique_ptr<Pythia8::Pythia> pythia_;       ///< Pythia 8 core to be wrapped
    const std::shared_ptr<Pythia8::CepGenEvent> cg_evt_;  ///< Event interface between CepGen and Pythia

    const bool correct_central_;
    const bool debug_lhef_;
    const std::string output_config_;
    bool res_decay_{true};
    bool enable_hadr_{false};
    unsigned short offset_{0};
    bool first_evt_{true};
  };

  void Pythia8Hadroniser::initialise() {
    cg_evt_->initialise(runParameters());
#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER < 8300
    pythia_->setLHAupPtr(cg_evt_.get());
#else
    pythia_->setLHAupPtr(cg_evt_);
#endif
    const auto& kin = runParameters().kinematics();

    pythia_->settings.parm("Beams:idA", (long)kin.incomingBeams().positive().integerPdgId());
    pythia_->settings.parm("Beams:idB", (long)kin.incomingBeams().negative().integerPdgId());
    // specify we will be using a LHA input
    pythia_->settings.mode("Beams:frameType", 5);
    pythia_->settings.parm("Beams:eCM", kin.incomingBeams().sqrtS());
    min_ids_ = kin.minimumFinalState();
    if (debug_lhef_)
      cg_evt_->openLHEF("debug.lhe");
    pythia_->settings.flag("ProcessLevel:resonanceDecays", res_decay_);
    if (pythia_->settings.flag("ProcessLevel:all") != enable_hadr_)
      pythia_->settings.flag("ProcessLevel:all", enable_hadr_);

    if (seed_ == -1ll)
      pythia_->settings.flag("Random:setSeed", false);
    else {
      pythia_->settings.flag("Random:setSeed", true);
      pythia_->settings.mode("Random:seed", seed_);
    }

#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER >= 8226
    switch (kin.incomingBeams().mode()) {
      case mode::Kinematics::ElasticElastic: {
        pythia_->settings.mode("BeamRemnants:unresolvedHadron", 3);
        pythia_->settings.flag("PartonLevel:MPI", false);
      } break;
      case mode::Kinematics::InelasticElastic: {
        pythia_->settings.mode("BeamRemnants:unresolvedHadron", 2);
        pythia_->settings.flag("PartonLevel:MPI", false);
      } break;
      case mode::Kinematics::ElasticInelastic: {
        pythia_->settings.mode("BeamRemnants:unresolvedHadron", 1);
        pythia_->settings.flag("PartonLevel:MPI", false);
      } break;
      case mode::Kinematics::InelasticInelastic:
      default: {
        pythia_->settings.mode("BeamRemnants:unresolvedHadron", 0);
      } break;
    }
#else
    CG_WARNING("Pythia8Hadroniser") << "Beam remnants framework for this version of Pythia "
                                    << "(" << utils::format("%.3f", pythia_->settings.parm("Pythia:versionNumber"))
                                    << ")\n\t"
                                    << "does not support mixing of unresolved hadron states.\n\t"
                                    << "The proton remnants output might hence be wrong.\n\t"
                                    << "Please update the Pythia version or disable this part.";
#endif
    if (correct_central_ && res_decay_)
      CG_WARNING("Pythia8Hadroniser") << "Central system's kinematics correction enabled while resonances are\n\t"
                                      << "expected to be decayed. Please check that this is fully intended.";

    if (!pythia_->init())
      throw CG_FATAL("Pythia8Hadroniser") << "Failed to initialise the Pythia8 core!\n\t"
                                          << "See the message above for more details.";

    if (debug_lhef_)
      cg_evt_->initLHEF();
  }

  bool Pythia8Hadroniser::run(Event& ev, double& weight, bool fast) {
    //--- initialise the event weight before running any decay algorithm
    weight = 1.;

    //--- only launch Pythia if:
    // 1) the full event kinematics (i.e. with remnants) is to be specified,
    // 2) the remnants are to be fragmented, or
    // 3) the resonances are to be decayed.
    if (!fast && !fragment_remnants_ && !res_decay_)
      return true;
    if (fast && !res_decay_)
      return true;

    //--- switch full <-> partial event
    if (!fast != enable_hadr_) {
      enable_hadr_ = !fast;
      initialise();
    }

    //===========================================================================================
    // convert our event into a custom LHA format
    //===========================================================================================

    cg_evt_->feedEvent(
        ev, fast ? Pythia8::CepGenEvent::Type::centralAndPartons : Pythia8::CepGenEvent::Type::centralAndBeamRemnants);
    if (debug_lhef_ && !fast)
      cg_evt_->eventLHEF();

    //===========================================================================================
    // launch the hadronisation / resonances decays, and update the event accordingly
    //===========================================================================================

    auto& num_hadr_trials = ev.metadata["pythia8:num_hadronisation_trials"];
    num_hadr_trials = 0;
    while (true) {
      if (num_hadr_trials++ > max_trials_)
        return false;
      //--- run the hadronisation/fragmentation algorithm
      if (pythia_->next()) {
        //--- hadronisation successful
        if (first_evt_ && !fast) {
          offset_ = 0;
          for (unsigned short i = 1; i < pythia_->event.size(); ++i)
            if (pythia_->event[i].status() == -PYTHIA_STATUS_IN_BEAM)
              //--- no incoming particles in further stages
              offset_++;
          first_evt_ = false;
        }
        break;
      }
    }
    CG_DEBUG("Pythia8Hadroniser") << "Pythia8 hadronisation performed successfully.\n\t"
                                  << "Number of trials: " << num_hadr_trials << "/" << max_trials_ << ".\n\t"
                                  << "Particles multiplicity: " << ev.particles().size() << " → "
                                  << pythia_->event.size() << ".\n\t"
                                  << "  indices offset: " << offset_ << ".";

    //===========================================================================================
    // update the event content with Pythia's output
    //===========================================================================================

    updateEvent(ev, weight);
    return true;
  }

  Particle& Pythia8Hadroniser::addParticle(Event& ev,
                                           const Pythia8::Particle& py_part,
                                           const Pythia8::Vec4& mom,
                                           unsigned short role) const {
    ParticleProperties prop;
    const pdgid_t pdg_id = py_part.idAbs();
    //--- define the particle if not already in the list of handled PDGs
    try {
      prop = PDG::get()(pdg_id);
    } catch (const Exception&) {
      prop.pdgid = pdg_id;
      prop.name = prop.human_name = py_part.name();
      prop.colours = py_part.col();  // colour factor
      prop.mass = py_part.m0();
      prop.width = py_part.mWidth();
      if (const auto ch = int(py_part.charge() * 3.); std::abs(ch) > 0)
        prop.charges = {ch, -ch};
      prop.fermion = py_part.isLepton();
      PDG::get().define(prop);
    }
    //--- add the particle to the event content
    Particle& op = ev.addParticle((Particle::Role)role);
    op.setPdgId((long)py_part.id());
    op.setStatus(py_part.isFinal()                                       ? Particle::Status::FinalState
                 : (Particle::Role)role == Particle::Role::CentralSystem ? Particle::Status::Propagator
                                                                         : Particle::Status::Fragmented);
    op.setMomentum(Momentum(mom.px(), mom.py(), mom.pz(), mom.e()).setMass(mom.mCalc()));
    cg_evt_->addCorresp(py_part.index() - offset_, op.id());
    return op;
  }

  void Pythia8Hadroniser::updateEvent(Event& ev, double& weight) const {
    std::vector<unsigned short> central_parts;

    for (unsigned short i = 1 + offset_; i < pythia_->event.size(); ++i) {
      const Pythia8::Particle& p = pythia_->event[i];
      const unsigned short cg_id = cg_evt_->cepgenId(i - offset_);
      if (cg_id != Pythia8::CepGenEvent::INVALID_ID) {
        //----- particle already in the event
        Particle& cg_part = ev[cg_id];
        //--- fragmentation result
        if (cg_part.role() == Particle::OutgoingBeam1 || cg_part.role() == Particle::OutgoingBeam2) {
          cg_part.setStatus(Particle::Status::Fragmented);
          continue;
        }
        //--- resonance decayed; apply branching ratio for this decay
        if (cg_part.role() == Particle::CentralSystem && p.status() < 0) {
          if (res_decay_)
            weight *= p.particleDataEntry().pickChannel().bRatio();
          cg_part.setStatus(Particle::Status::Resonance);
          central_parts.emplace_back(i);
        }
        //--- particle is not what we expect
        if (p.idAbs() != std::labs(cg_part.integerPdgId())) {
          CG_INFO("Pythia8Hadroniser:update") << "LHAEVT event content:";
          cg_evt_->listEvent();
          CG_INFO("Pythia8Hadroniser:update") << "Pythia event content:";
          pythia_->event.list();
          CG_INFO("Pythia8Hadroniser:update") << "CepGen event content:";
          ev.dump();
          CG_INFO("Pythia8Hadroniser:update") << "Correspondence:";
          cg_evt_->dumpCorresp();

          throw CG_FATAL("Pythia8Hadroniser:update")
              << "Event list corruption detected for (Pythia/CepGen) particle " << i << "/" << cg_id << ":\n\t"
              << "should be " << abs(p.id()) << ", "
              << "got " << cg_part.integerPdgId() << "!";
        }
      }
      //--- check for messed up particles parentage and discard incoming beam particles
      /*else if ( p.mother1() > i || p.mother1() <= offset_ )
          continue;
        else if ( p.mother2() > i || p.mother2() <= offset_ )
          continue;*/
      else {
        //----- new particle to be added
        const unsigned short role = findRole(ev, p);
        switch ((Particle::Role)role) {
          case Particle::OutgoingBeam1:
            ev[Particle::OutgoingBeam1][0].get().setStatus(Particle::Status::Fragmented);
            break;
          case Particle::OutgoingBeam2:
            ev[Particle::OutgoingBeam2][0].get().setStatus(Particle::Status::Fragmented);
            break;
          default:
            break;
        }
        // found the role ; now we can add the particle
        Particle& cg_part = addParticle(ev, p, p.p(), role);
        if (correct_central_ && (Particle::Role)role == Particle::CentralSystem) {
          if (const auto ip = std::find(central_parts.begin(), central_parts.end(), p.mother1());
              ip != central_parts.end())
            cg_part.setMomentum(ev[cg_evt_->cepgenId(*ip - offset_)].momentum());
        }
        for (const auto& moth_id : p.motherList()) {
          if (moth_id <= offset_)
            continue;
          if (const unsigned short moth_cg_id = cg_evt_->cepgenId(moth_id - offset_);
              moth_cg_id != Pythia8::CepGenEvent::INVALID_ID)
            cg_part.addMother(ev[moth_cg_id]);
          else
            cg_part.addMother(addParticle(ev, pythia_->event[moth_id], p.p(), role));
          if (!p.isFinal()) {
            if (p.isResonance() || !p.daughterList().empty())
              cg_part.setStatus(Particle::Status::Resonance);
            else
              cg_part.setStatus(Particle::Status::Undefined);
          }
        }
      }
    }
  }

  unsigned short Pythia8Hadroniser::findRole(const Event& ev, const Pythia8::Particle& p) const {
    for (const auto& par_id : p.motherList()) {
      if (par_id == 1 && offset_ > 0)
        return (unsigned short)Particle::OutgoingBeam1;
      if (par_id == 2 && offset_ > 0)
        return (unsigned short)Particle::OutgoingBeam2;
      if (const unsigned short par_cg_id = cg_evt_->cepgenId(par_id - offset_);
          par_cg_id != Pythia8::CepGenEvent::INVALID_ID)
        return (unsigned short)ev(par_cg_id).role();
      if (par_id != Pythia8::CepGenEvent::INVALID_ID)
        return findRole(ev, pythia_->event[par_id]);
    }
    return (unsigned short)Particle::UnknownRole;
  }

  ParametersDescription Pythia8Hadroniser::description() {
    auto desc = Hadroniser::description();
    desc.setDescription("Interface to the Pythia 8 string hadronisation/fragmentation algorithm");
    desc.add<bool>("correctCentralSystem", false)
        .setDescription("Correct the kinematics of the central system whenever required");
    desc.add<bool>("debugLHEF", false).setDescription("Switch on the dump of each event into a debugging LHEF file");
    desc.add<std::string>("outputConfig", "last_pythia_config.cmd")
        .setDescription("Output filename for a backup of the last Pythia configuration snapshot");
    return desc;
  }
}  // namespace cepgen::hadr
// register hadroniser
using cepgen::hadr::Pythia8Hadroniser;
REGISTER_MODIFIER("pythia8", Pythia8Hadroniser);
