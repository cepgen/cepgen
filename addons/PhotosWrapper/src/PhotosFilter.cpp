/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2026  Laurent Forthomme
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

#include <Photos/Log.h>
#include <Photos/Photos.h>
#include <Photos/PhotosEvent.h>
#include <Photos/PhotosHepMC3Event.h>

#include <set>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace Photospp;

namespace cepgen::photos {
  /// Interface to the Photos decay routine
  class PhotosFilter final : public EventModifier {
  public:
    explicit PhotosFilter(const ParametersList& params) : EventModifier(params) {
      if (steer<bool>("debug"))
        Log::LogAll(true);
      Photos::setMomentumConservationThreshold(steer<double>("momentumConservationThreshold"));
      if (seed_ >= 0)  // must happen before Photos::initialize(), which initialises the generator from these seeds
        Photos::setSeed(static_cast<int>(seed_ % 31328ll), static_cast<int>((seed_ / 31328ll) % 30081ll));
      Photos::setInterference(steer<bool>("interference"));
      Photos::setDoubleBrem(steer<bool>("doubleBrem"));
      Photos::setQuatroBrem(steer<bool>("quatroBrem"));
      Photos::setCorrectionWtForW(steer<bool>("correctionWtForW"));
      Photos::setExponentiation(steer<bool>("exponentiation"));
      Photos::setPairEmission(steer<bool>("pairEmission"));
      Photos::setPhotonEmission(steer<bool>("photonEmission"));
      Photos::setMeCorrectionWtForScalar(steer<bool>("meCorrectionWtForScalar"));
      Photos::setMeCorrectionWtForW(steer<bool>("meCorrectionWtForW"));
      Photos::setMeCorrectionWtForZ(steer<bool>("meCorrectionWtForZ"));
      Photos::setTopProcessRadiation(steer<bool>("topProcessRadiation"));
    }
    ~PhotosFilter() override { Log::SummaryAtExit(); }

    static ParametersDescription description() {
      auto desc = EventModifier::description();
      desc.add("debug", false).setDescription("log all debugging information?");
      desc.add("momentumConservationThreshold", 1.e-6).setDescription("numerical limit for momentum conservation");
      desc.add("maxWtInterference", -1.).setDescription("maximum interference weight (< 0: use the Photos++ default)");
      desc.add("infraredCutOff", -1.)
          .setDescription(
              "minimal energy (in units of decaying particle mass) for photons to be explicitly generated"
              " (<= 0: use the Photos++ default for the selected mode)");
      desc.add("interference", true).setDescription("key for interference, matrix element weight");
      desc.add("doubleBrem", true).setDescription("set double bremsstrahlung generation (ignored with exponentiation)");
      desc.add("quatroBrem", false)
          .setDescription("set bremsstrahlung generation up to multiplicity of 4 (ignored with exponentiation)");
      desc.add("correctionWtForW", true)
          .setDescription("key for partial effects of matrix element (in leptonic W decays)");
      desc.add("exponentiation", true).setDescription("set exponentiation mode");
      desc.add("pairEmission", false).setDescription("set pair emission");
      desc.add("photonEmission", true).setDescription("set photon emission");
      desc.add("meCorrectionWtForScalar", false)
          .setDescription("switch for complete effects of the matrix element (in scalar to two scalar decays)");
      desc.add("meCorrectionWtForW", false)
          .setDescription("switch for complete effects of matrix element (in leptonic W decays)");
      desc.add("meCorrectionWtForZ", false)
          .setDescription("switch for complete effects of matrix element (in leptonic Z decays)");
      desc.add("topProcessRadiation", true)
          .setDescription("set photon emission in top pair production in quark (gluon) pair annihilation");
      return desc;
    }

    void initialise() override {
      Photos::setMomentumUnit(Photos::GEV);
      Photos::setAlphaQED(constants::ALPHA_EM);
      Photos::initialize();
      // these are overwritten by Photos::initialize(), and hence can only be steered afterwards
      if (const auto max_wt_interference = steer<double>("maxWtInterference"); max_wt_interference >= 0.)
        Photos::maxWtInterference(max_wt_interference);
      if (const auto infrared_cut_off = steer<double>("infraredCutOff"); infrared_cut_off > 0.)
        Photos::setInfraredCutOff(infrared_cut_off);
    }
    bool run(Event& event, double& weight, bool /*fast*/) override {
      weight = 1.;
      // radiation off the beam -> parton + remnant lines would alter the partons feeding the central system
      for (const auto role : {Particle::Role::IncomingBeam1, Particle::Role::IncomingBeam2})
        if (event.hasRole(role))
          if (const auto beam_id = event.oneWithRole(role).integerPdgId(); suppressed_beams_.insert(beam_id).second)
            Photos::suppressBremForDecay(0, static_cast<int>(beam_id));
      CepGenEvent hepmc_event(event);  // conversion to a HepMC3 format
      PhotosHepMC3Event photos_event(&hepmc_event);
      photos_event.process();
      return hepmc_event.merge(event);  // merge everything back into the original event
    }

  private:
    std::set<long> suppressed_beams_;  ///< beam particles for which Photos++ radiation was already disabled
  };
}  // namespace cepgen::photos
using cepgen::photos::PhotosFilter;
REGISTER_MODIFIER("photos", PhotosFilter);
