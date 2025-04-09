/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2022  Laurent Forthomme
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
  class PhotosFilter : public EventModifier {
  public:
    explicit PhotosFilter(const ParametersList&);
    ~PhotosFilter();

    static ParametersDescription description();

    void initialise() override;
    bool run(Event& ev, double& weight, bool fast) override;
  };

  PhotosFilter::PhotosFilter(const ParametersList& params) : EventModifier(params) {
    if (steer<bool>("debug"))
      Log::LogAll(true);
    Photos::setMomentumConservationThreshold(1.e-10);
    Photos::maxWtInterference(steer<double>("maxWtInterference"));
    Photos::setInfraredCutOff(steer<double>("infraredCutOff"));
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

  PhotosFilter::~PhotosFilter() { Log::SummaryAtExit(); }

  void PhotosFilter::initialise() {
    Photos::setMomentumUnit(Photos::GEV);
    Photos::setAlphaQED(constants::ALPHA_EM);
    Photos::initialize();
  }

  bool PhotosFilter::run(Event& event, double& weight, bool) {
    weight = 1.;
    HepMC3::CepGenEvent hepmc_event(event);
    PhotosHepMC3Event photos_event(&hepmc_event);
    //event.dump();
    photos_event.process();
    hepmc_event.merge(event);
    event.dump();
    return true;
  }

  ParametersDescription PhotosFilter::description() {
    auto desc = EventModifier::description();
    desc.add("debug", false).setDescription("log all debugging information?");
    desc.add("maxWtInterference", 1.).setDescription("maximum interference weight");
    desc.add("infraredCutOff", 0.01)
        .setDescription("minimal energy (in units of decaying particle mass) for photons to be explicitly generated");
    desc.add("interference", true).setDescription("key for interference, matrix element weight");
    desc.add("doubleBrem", true).setDescription("set double bremsstrahlung generation");
    desc.add("quatroBrem", false).setDescription("set bremsstrahlung generation up to multiplicity of 4");
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
}  // namespace cepgen::photos

// register event modifier
using cepgen::photos::PhotosFilter;
REGISTER_MODIFIER("photos", PhotosFilter);
