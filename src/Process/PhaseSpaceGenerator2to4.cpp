/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Process/CentralPhaseSpace.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;
using namespace std::string_literals;

/// Two-parton (with fluxes) to two-particle central system phase space generator
/// \note Transitional wrapper: composes a two-parton generator and a 2-to-2 central system phase space mapper
class PhaseSpaceGenerator2to4 : public PhaseSpaceGenerator {
public:
  explicit PhaseSpaceGenerator2to4(const ParametersList& params)
      : PhaseSpaceGenerator(params),
        part_psgen_(PartonsPhaseSpaceGeneratorFactory::get().build(steer<std::string>("partonsGenerator"), params_)),
        central_psgen_(std::make_unique<CentralPhaseSpace2to2>(params_)) {}

  static ParametersDescription description() {
    auto desc = PhaseSpaceGenerator::description();
    desc.setDescription("2-to-4 phase space mapper");
    desc.add("partonsGenerator", ""s).setDescription("type of partons generator algorithm to use");
    return desc;
  }

  bool ktFactorised() const override {
    CG_ASSERT(part_psgen_);
    return part_psgen_->ktFactorised();
  }

  void initialise(proc::FactorisedProcess* process) override {
    proc_ = process;
    CG_ASSERT(part_psgen_);
    part_psgen_->initialise(process);  // registers the partons' variables first, to preserve dimensions ordering
    central_psgen_->initialise(process->kinematics(), particles_);
    defineVariables(*process, central_psgen_->integrationVariables());
  }

  bool generate() override {
    CG_ASSERT(part_psgen_);
    if (!part_psgen_->generatePartonKinematics())
      return false;
    central_weight_ = central_psgen_->generate((proc_->q1() + proc_->q2()).transverse(), central_momenta_);
    if (!utils::positive(central_weight_))
      return false;
    for (size_t i = 0; i < central_momenta_.size(); ++i)
      proc_->pc(i) = central_momenta_.at(i);
    return true;
  }

  double weight() const override {
    const auto fluxes_weight = part_psgen_->fluxes();
    if (!utils::positive(fluxes_weight))
      return 0.;
    return fluxes_weight * central_weight_;
  }

  spdgids_t partons() const override {
    CG_ASSERT(part_psgen_);
    return {part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
  }

  spdgids_t central() const override { return particles_; }

  void setCentral(const spdgids_t& central_pdgids) override {
    if (central_pdgids.size() != central_psgen_->multiplicity())
      throw CG_FATAL("PhaseSpaceGenerator2to4:setCentral")
          << "Invalid central particles multiplicity: expecting " << central_psgen_->multiplicity() << ", got "
          << central_pdgids << ".";
    particles_ = central_pdgids;
  }

private:
  const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
  const std::unique_ptr<CentralPhaseSpace> central_psgen_;
  spdgids_t particles_;  ///< Type of particles produced in the final state (integer values)

  proc::FactorisedProcess* proc_{nullptr};  //NOT owning

  std::vector<Momentum> central_momenta_;
  double central_weight_{0.};
};
REGISTER_PHASE_SPACE_GENERATOR("2to4", PhaseSpaceGenerator2to4);
