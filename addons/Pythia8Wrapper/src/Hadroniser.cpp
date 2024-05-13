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
#include "CepGenPythia8/EventInterface.h"

namespace cepgen::pythia8 {
  /// Interface to the Pythia8 hadronisation algorithm
  /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
  class Hadroniser : public cepgen::hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& plist)
        : cepgen::hadr::Hadroniser(plist),
          pythia_(new Pythia8::Pythia),
          cg_evt_(new pythia8::EventInterface),
          correct_central_(steer<bool>("correctCentralSystem")),
          debug_lhef_(steer<bool>("debugLHEF")),
          output_config_(steer<std::string>("outputConfig")) {}

    inline virtual ~Hadroniser() {
      if (!output_config_.empty())
        pythia_->settings.writeFile(output_config_, false);
      if (debug_lhef_)
        cg_evt_->closeLHEF(true);
    }

    inline static ParametersDescription description() {
      auto desc = cepgen::hadr::Hadroniser::description();
      desc.setDescription("Interface to the Pythia 8 string hadronisation/fragmentation algorithm");
      desc.add<bool>("correctCentralSystem", false)
          .setDescription("correct any discrepancy of the kinematics of the central system?");
      desc.add<bool>("debugLHEF", false).setDescription("dump each event into a debugging LHEF file?");
      desc.add<std::string>("outputConfig", "last_pythia_config.cmd")
          .setDescription("Pythia configuration backup output filename");
      return desc;
    }

    inline void readString(const std::string& param) override {
      if (!pythia_->readString(param))
        throw CG_FATAL("pythia8:Hadroniser") << "The Pythia8 core failed to parse the following setting:\n\t" << param;
    }
    void initialise() override;
    bool run(Event& ev, double& weight, bool fast) override;

    inline void setCrossSection(const Value& cross_section) override { cg_evt_->setCrossSection(0, cross_section); }

  private:
    void* enginePtr() override { return (void*)pythia_.get(); }

    static constexpr unsigned short PYTHIA_STATUS_IN_BEAM = 12;
    static constexpr unsigned short PYTHIA_STATUS_IN_PARTON_KT = 61;

    pdgids_t min_ids_;
    std::unordered_map<short, short> py_cg_corresp_;

    const std::unique_ptr<Pythia8::Pythia> pythia_;          ///< Pythia 8 core to be wrapped
    const std::shared_ptr<pythia8::EventInterface> cg_evt_;  ///< Event interface between CepGen and Pythia

    const bool correct_central_;
    const bool debug_lhef_;
    const std::string output_config_;
    bool res_decay_{true};
    bool enable_hadr_{false};
    unsigned short offset_{0};
    bool first_evt_{true};
  };

  void Hadroniser::initialise() {
    cg_evt_->initialise(runParameters());
#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER < 8300
    pythia_->setLHAupPtr(cg_evt_.get());
#else
    pythia_->setLHAupPtr(cg_evt_);
#endif
    const auto& kin = runParameters().kinematics();

    pythia_->settings.flag("BeamRemnants:primordialKT", false);
    pythia_->settings.parm("Check:epTolErr", 1.);
    pythia_->settings.parm("Check:mTolErr", 1.);
    pythia_->settings.parm("Beams:idA", (long)kin.incomingBeams().positive().integerPdgId());
    pythia_->settings.parm("Beams:idB", (long)kin.incomingBeams().negative().integerPdgId());
    // specify we will be using a LHA input
    pythia_->settings.mode("Beams:frameType", 5);
    pythia_->settings.parm("Beams:eCM", kin.incomingBeams().sqrtS());
    //pythia_->settings.flag("Check:beams", false);  //FIXME
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
    pythia_->settings.flag("PartonLevel:ISR", false);
    pythia_->settings.flag("PartonLevel:FSR", false);
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
    CG_WARNING("pythia8:Hadroniser") << "Beam remnants framework for this version of Pythia "
                                     << "(" << utils::format("%.3f", pythia_->settings.parm("Pythia:versionNumber"))
                                     << ")\n\t"
                                     << "does not support mixing of unresolved hadron states.\n\t"
                                     << "The proton remnants output might hence be wrong.\n\t"
                                     << "Please update the Pythia version or disable this part.";
#endif
    if (correct_central_ && res_decay_)
      CG_WARNING("pythia8:Hadroniser") << "Central system's kinematics correction enabled while resonances are\n\t"
                                       << "expected to be decayed. Please check that this is fully intended.";

    if (!pythia_->init())
      throw CG_FATAL("pythia8:Hadroniser") << "Failed to initialise the Pythia8 core!\n\t"
                                           << "See the message above for more details.";

    if (debug_lhef_)
      cg_evt_->initLHEF();
  }

  bool Hadroniser::run(Event& ev, double& weight, bool fast) {
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

    cg_evt_->feedEvent(ev);  // convert our event into a custom LHA format
    if (debug_lhef_ && !fast)
      cg_evt_->eventLHEF();

    // launch the hadronisation / resonances decays, and update the event accordingly
    auto& num_hadr_trials = ev.metadata["pythia8:num_hadronisation_trials"];
    num_hadr_trials = 0;
    while (true) {  // run the hadronisation/fragmentation algorithm
      if (num_hadr_trials++ > max_trials_)
        return false;
      if (pythia_->next()) {        // hadronisation successful
        if (first_evt_ && !fast) {  // we build the association map between the CepGen and Pythia8 events
          for (unsigned short i = 1; i < pythia_->event.size(); ++i)
            if (pythia_->event[i].status() == -PYTHIA_STATUS_IN_BEAM)  // no incoming particles in later stages
              offset_++;
          first_evt_ = false;
        }
        break;
      }
    }
    CG_DEBUG("pythia8:Hadroniser") << "Pythia8 hadronisation performed successfully.\n\t"
                                   << "Number of trials: " << num_hadr_trials << "/" << max_trials_ << ".\n\t"
                                   << "Particles multiplicity: " << ev.particles().size() << " → "
                                   << pythia_->event.size() << ".\n\t"
                                   << "  indices offset: " << offset_ << ".";

    cg_evt_->updateEvent(pythia_->event, ev, weight);  // update the event content with Pythia's output
    CG_LOG << ev;
    return true;
  }
}  // namespace cepgen::pythia8
using Pythia8Hadroniser = cepgen::pythia8::Hadroniser;
REGISTER_MODIFIER("pythia8", Pythia8Hadroniser);
