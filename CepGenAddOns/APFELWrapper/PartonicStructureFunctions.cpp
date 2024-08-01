/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

namespace cepgen::apfel {
  /// Generic partonic level perturbative structure functions built from an external PDFs grid
  class PartonicStructureFunctions final : public strfun::PartonicParameterisation {
  public:
    /// Quarks types
    enum class Mode { full = 0, valence = 1, sea = 2 };
    /// Build a calculator from its Parameters object
    explicit PartonicStructureFunctions(const ParametersList& params)
        : PartonicParameterisation(params), q_limits_(steer<Limits>("qrange")), xbj_min_(steer<double>("xBjmin")) {
      const auto perturbative_order = steer<int>("perturbativeOrder");
      APFEL::SetPerturbativeOrder(perturbative_order);
      APFEL::InitializeAPFEL();
      APFEL::EvolveAPFEL(q_limits_.min(), q_limits_.max());
      APFEL::CachePDFsAPFEL(q_limits_.min());
      CG_INFO("apfel:PartonicStructureFunctions") << "Partonic structure functions evaluator successfully built.\n"
                                                  << " * APFEL version: " << APFEL::GetVersion() << "\n"
                                                  << " * number of flavours: " << num_flavours_ << "\n"
                                                  << " * quarks mode: " << mode_ << "\n"
                                                  << " * Q range: " << q_limits_ << ", min xBj: " << xbj_min_ << "\n"
                                                  << " * perturbative order: " << perturbative_order << ".";
    }

    static ParametersDescription description() {
      auto desc = PartonicParameterisation::description();
      desc.setDescription("APFEL (partonic)");
      desc.add<int>("perturbativeOrder", 2);
      desc.add<Limits>("qrange", {1., 100.});
      desc.add<double>("xBjmin", 2.e-6).setDescription("minimum Bjorken-x reachable for this PDF set");
      return desc;
    }

  private:
    double evalxQ2(int flavour, double xbj, double q2) override {
      if (xbj < xbj_min_)
        return 0.;
      const auto q = std::sqrt(q2);
      if (!q_limits_.contains(q))
        return 0.;
      return APFEL::xPDFxQ(flavour, xbj, q);
    }
    const Limits q_limits_;
    const double xbj_min_;
  };
}  // namespace cepgen::apfel
using PartonicStructureFunctionsAPFEL = cepgen::apfel::PartonicStructureFunctions;
REGISTER_STRFUN("apfel", 402, PartonicStructureFunctionsAPFEL);
