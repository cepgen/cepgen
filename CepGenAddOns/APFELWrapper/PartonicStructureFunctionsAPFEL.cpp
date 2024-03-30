/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include <APFEL/APFEL.h>

#include <cmath>

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/PartonicParameterisation.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace strfun {
    /// Generic partonic level perturbative structure functions built from an external PDFs grid
    class APFELPartonic final : public PartonicParameterisation {
    public:
      /// Quarks types
      enum class Mode { full = 0, valence = 1, sea = 2 };
      /// Build a calculator from its Parameters object
      explicit APFELPartonic(const ParametersList& params)
          : PartonicParameterisation(params), q_limits_(steer<Limits>("qLimits")) {
        const auto perturbative_order = steer<int>("perturbativeOrder");
        APFEL::SetPerturbativeOrder(perturbative_order);
        APFEL::InitializeAPFEL();
        APFEL::EvolveAPFEL(q_limits_.min(), q_limits_.max());
        APFEL::CachePDFsAPFEL(q_limits_.min());
        CG_INFO("APFELPartonic") << "Partonic structure functions evaluator successfully built.\n"
                                 << " * APFEL version: " << APFEL::GetVersion() << "\n"
                                 << " * number of flavours: " << num_flavours_ << "\n"
                                 << " * quarks mode: " << mode_ << "\n"
                                 << " * Q range: " << q_limits_ << "\n"
                                 << " * perturbative order: " << perturbative_order << ".";
      }

      static int index() { return 402; }

      static ParametersDescription description() {
        auto desc = PartonicParameterisation::description();
        desc.setDescription("APFEL (partonic)");
        desc.add<int>("perturbativeOrder", 2);
        desc.add<Limits>("qLimits", {1., 100.});
        return desc;
      }

    private:
      double evalxQ2(int flavour, double xbj, double q2) override {
        const auto q = std::sqrt(q2);
        if (!q_limits_.contains(q))
          return 0.;
        return APFEL::xPDFxQ(flavour, xbj, q);
      }
      const Limits q_limits_;
    };
  }  // namespace strfun
}  // namespace cepgen
using cepgen::strfun::APFELPartonic;
REGISTER_STRFUN("apfel", APFELPartonic);
