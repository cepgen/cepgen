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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6EventInterface.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

namespace cepgen {
  namespace hadr {
    /**
     * Full interface to the Pythia 6 algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia 6 hadronisation algorithm
     */
    class Pythia6Hadroniser : public Hadroniser {
    public:
      explicit Pythia6Hadroniser(const ParametersList& params) : Hadroniser(params) {}

      static ParametersDescription description() {
        auto desc = Hadroniser::description();
        desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
        return desc;
      }

      inline void readString(const std::string& param) override { pythia6::pygive(param); }
      void initialise() override {
        CG_WARNING("Pythia6Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                        << "You will have to specify manually the multiplication factor according\n\t"
                                        << "to your list of open channels.";
      }
      bool run(Event& ev, double& weight, bool full) override;

    private:
      Pythia6EventInterface evt_;
    };

    bool Pythia6Hadroniser::run(Event& ev, double& weight, bool full) {
      weight = 1.;
      evt_.feedEvent(ev);  // fill Pythia 6 common blocks

      CG_DEBUG_LOOP("Pythia6Hadroniser") << "Dump of the event before the hadronisation:" << ev << "\n\t"
                                         << utils::s("string object", evt_.numStrings(), true)
                                         << " identified and constructed.";

      const int old_npart = pyjets_.n;

      evt_.run();  // run the hadronisation/decay

      if (full && pyjets_.n == old_npart)
        return false;  // hadronisation failed

      //--- update the event
      for (int p = old_npart; p < pyjets_.n; ++p) {
        // filter the first particles already present in the event
        const pdgid_t pdg_id = abs(pyjets_.k[1][p]);
        ParticleProperties prop;
        if (full)
          if (!PDG::get().has(pdg_id)) {
            const auto name = pythia6::pyname(pdg_id);
            prop.pdgid = pdg_id;
            prop.name = name;
            prop.descr = name;
            prop.colours = pythia6::pyk(p + 1, 12);  // colour factor
            prop.mass = pythia6::pymass(pdg_id);
            prop.width = -1.;                      //pmas( pdg_id, 2 ),
            prop.charge = pythia6::pyk(p + 1, 6);  // charge
            prop.fermion = false;
            PDG::get().define(prop);
          }

        const auto moth_id = pyjets_.k[2][p] - 1;
        const auto role = pyjets_.k[2][p] != 0 ? ev[moth_id].role()  // child particle inherits its mother's role
                                               : Particle::Role::UnknownRole;

        auto& pa = ev.addParticle(role).get();
        pa.setId(p);
        pa.setStatus(pyjets_.k[0][p]);
        pa.setPdgId((long)pyjets_.k[1][p]);
        pa.setMomentum(
            Momentum(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]).setMass(pyjets_.p[4][p]));
        auto& moth = ev[moth_id];
        if (role != Particle::Role::UnknownRole)
          moth.setStatus(role == Particle::Role::CentralSystem ? Particle::Status::Resonance
                                                               : Particle::Status::Fragmented);
        pa.addMother(moth);
      }
      return true;
    }
  }  // namespace hadr
}  // namespace cepgen
using cepgen::hadr::Pythia6Hadroniser;
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser);
