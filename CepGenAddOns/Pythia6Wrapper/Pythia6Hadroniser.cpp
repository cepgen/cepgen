/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
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
      explicit Pythia6Hadroniser(const ParametersList& params) : Hadroniser(params) {}

      static inline ParametersDescription description() {
        auto desc = Hadroniser::description();
        desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
        return desc;
      }

      inline void readString(const std::string& param) override { pythia6::pygive(param); }
      inline void initialise() override {
        CG_WARNING("Pythia6Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                        << "You will have to specify manually the multiplication factor according\n\t"
                                        << "to your list of open channels.";
      }
      inline bool run(Event& ev, double& weight, bool full) override {
        weight = 1.;
        pythia6::EventInterface evt(ev);
        evt.prepareHadronisation();  // fill Pythia 6 common blocks

        CG_DEBUG_LOOP("Pythia6Hadroniser")
            << "Dump of the event before the hadronisation:" << ev << "\n\t"
            << utils::s("string object", evt.numStrings(), true) << " identified and constructed.";

        const int old_npart = pyjets_.n;
        evt.run();  // run the hadronisation/decay
        if (full && pyjets_.n == old_npart)
          return false;  // hadronisation failed

        return true;
      }
    };
  }  // namespace hadr
}  // namespace cepgen
using cepgen::hadr::Pythia6Hadroniser;
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser);
