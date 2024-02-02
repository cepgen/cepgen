/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include "CepGenAddOns/Pythia6Wrapper/EventInterface.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

namespace cepgen {
  namespace hadr {
    /**
     * Interface to the Pythia 6 algorithm
     * \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     */
    class Pythia6Hadroniser : public Hadroniser {
    public:
      explicit Pythia6Hadroniser(const ParametersList& params)
          : Hadroniser(params),
            rnd_gen_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

      static inline ParametersDescription description() {
        auto desc = Hadroniser::description();
        desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
        desc.add<ParametersDescription>("randomGenerator", ParametersDescription().setName<std::string>("stl"))
            .setDescription("random number generator to use for the various intermediate computations");
        return desc;
      }

      inline void readString(const std::string& param) override { pythia6::pygive(param); }
      inline void initialise() override {
        CG_WARNING("Pythia6Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                        << "You will have to specify manually the multiplication factor according\n\t"
                                        << "to your list of open channels.";
        kin_mode_ = runParameters().kinematics().incomingBeams().mode();
      }
      inline bool run(Event& ev, double& weight, bool fast) override {
        weight = 1.;
        pythia6::EventInterface evt(
            ev,
            fast ? mode::Kinematics::ElasticElastic  // do not treat beam remnants when running in fast mode
                 : kin_mode_,
            rnd_gen_.get());
        evt.prepareHadronisation();  // fill Pythia 6 common blocks

        CG_DEBUG_LOOP("Pythia6Hadroniser")
            << "Dump of the event before the hadronisation:" << ev << "\n\t"
            << utils::s("string object", evt.numStrings(), true) << " identified and constructed.";

        const int old_npart = pyjets_.n;
        evt.run();  // run the hadronisation/decay
        if (!fast && pyjets_.n == old_npart)
          return false;  // hadronisation failed

        return true;
      }

    private:
      mode::Kinematics kin_mode_;
      std::unique_ptr<utils::RandomGenerator> rnd_gen_;
    };
  }  // namespace hadr
}  // namespace cepgen
using cepgen::hadr::Pythia6Hadroniser;
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser);
