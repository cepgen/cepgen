/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"
#include "CepGenPythia6/EventInterface.h"
#include "CepGenPythia6/Pythia6Interface.h"

namespace cepgen::pythia6 {
  /// Interface to the Pythia 6 algorithm
  /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
  class Hadroniser final : public cepgen::hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& params)
        : cepgen::hadr::Hadroniser(params),
          random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

    static ParametersDescription description() {
      auto desc = cepgen::hadr::Hadroniser::description();
      desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
      desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
          .setDescription("random number generator to use for the various intermediate computations");
      return desc;
    }

    void readString(const std::string& param) override { pygive(param); }
    void initialise() override {
      CG_WARNING("pythia6:Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                       << "You will have to specify manually the multiplication factor according\n\t"
                                       << "to your list of open channels.";
      kinematics_mode_ = runParameters().kinematics().incomingBeams().mode();
    }
    bool run(Event& event, double& weight, bool fast) override {
      weight = 1.;
      EventInterface pythia_event(
          event,
          fast ? mode::Kinematics::ElasticElastic /*do not treat beam remnants when running in fast mode*/
               : kinematics_mode_,
          random_generator_.get());
      pythia_event.prepareHadronisation();  // fill Pythia 6 common blocks

      CG_DEBUG_LOOP("pythia6:Hadroniser")
          << "Dump of the event before the hadronisation:" << event << "\n\t"
          << utils::s("string object", pythia_event.numStrings(), true) << " identified and constructed.";

      const auto old_particles_multiplicity = pyjets_.n;
      pythia_event.run();  // run the hadronisation/decay
      if (!fast && pyjets_.n == old_particles_multiplicity)
        return false;  // hadronisation failed
      return true;
    }

  private:
    mode::Kinematics kinematics_mode_{mode::Kinematics::ElasticElastic};
    const std::unique_ptr<utils::RandomGenerator> random_generator_;
  };
}  // namespace cepgen::pythia6
using Pythia6Hadroniser = cepgen::pythia6::Hadroniser;
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser);
