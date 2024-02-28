/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Process/Central2to4PhaseSpaceGenerator.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsCollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

namespace cepgen {
  template <typename Tp, typename Tc>
  class FactorisedPhaseSpaceGenerator : public PhaseSpaceGenerator {
  public:
    explicit FactorisedPhaseSpaceGenerator(const ParametersList& params)
        : PhaseSpaceGenerator(params), part_psgen_(new Tp(params)), cent_psgen_(new Tc(params)) {}

    static ParametersDescription description() {
      auto desc = PhaseSpaceGenerator::description();
      desc.setDescription("Factorised parton/central phase space mapper (" + Tp::description().description() + "/" +
                          Tc::description().description() + ")");
      desc += Tp::description();
      desc += Tc::description();
      return desc;
    }

    bool ktFactorised() const override {
      CG_ASSERT(part_psgen_);
      return part_psgen_->ktFactorised();
    }

    void setCentralCuts(const cuts::Central& cuts) const override {
      CG_ASSERT(cent_psgen_);
      cent_psgen_->setCuts(cuts);
    }

    void initialise(proc::FactorisedProcess* process) override {
      CG_ASSERT(part_psgen_);
      CG_ASSERT(cent_psgen_);
      part_psgen_->initialise(process);
      cent_psgen_->initialise(process);
    }

    double generate() override {
      CG_ASSERT(part_psgen_);
      CG_ASSERT(cent_psgen_);
      if (!part_psgen_->generatePartonKinematics())
        return 0.;
      const auto cent_weight = cent_psgen_->generateKinematics();
      if (!utils::positive(cent_weight))
        return 0.;
      const auto fluxes_weight = part_psgen_->fluxes();
      if (!utils::positive(fluxes_weight))
        return 0.;
      return fluxes_weight * cent_weight;
    }

    pdgids_t partons() const override {
      CG_ASSERT(part_psgen_);
      return pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
    }
    pdgids_t central() const override {
      CG_ASSERT(cent_psgen_);
      return cent_psgen_->particles();
    }

  private:
    const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
    const std::unique_ptr<CentralPhaseSpaceGenerator> cent_psgen_;
  };
  typedef FactorisedPhaseSpaceGenerator<PartonsKTPhaseSpaceGenerator, Central2to4PhaseSpaceGenerator> KT2to4;
  typedef FactorisedPhaseSpaceGenerator<PartonsCollinearPhaseSpaceGenerator, Central2to4PhaseSpaceGenerator> Coll2to4;
}  // namespace cepgen

REGISTER_PSGEN("kt2to4", KT2to4);
REGISTER_PSGEN("coll2to4", Coll2to4);
