/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"

using namespace cepgen;

PartonsPhaseSpaceGenerator::PartonsPhaseSpaceGenerator(const ParametersList& params,
                                                       std::unique_ptr<PartonSampler> positive_sampler,
                                                       std::unique_ptr<PartonSampler> negative_sampler)
    : SteeredObject(params), samplers_{std::move(positive_sampler), std::move(negative_sampler)} {
  for (const auto& sampler : samplers_)
    if (!sampler)
      throw CG_FATAL("PartonsPhaseSpaceGenerator") << "Invalid single-beam parton sampler.";
}

void PartonsPhaseSpaceGenerator::initialise(proc::FactorisedProcess* process) {
  proc_ = process;
  const auto& kin = this->process().kinematics();
  samplers_.at(0)->initialise(kin, PartonSampler::Side::positive);
  samplers_.at(1)->initialise(kin, PartonSampler::Side::negative);

  // register the samplers' integration variables, interleaved between beams
  // (i.e. var1(+z), var1(-z), var2(+z), var2(-z), ...) to preserve the historical dimensions ordering
  const auto& initial_cuts = kin.cuts().initial;
  const auto pos_vars = samplers_.at(0)->integrationVariables(initial_cuts),
             neg_vars = samplers_.at(1)->integrationVariables(initial_cuts);
  auto define = [this](const PartonSampler::IntegrationVariable& var) {
    this->process().defineVariable(var.value.get(), var.mapping, var.range, var.name, var.description);
  };
  for (size_t i = 0; i < std::max(pos_vars.size(), neg_vars.size()); ++i) {
    if (i < pos_vars.size())
      define(pos_vars.at(i));
    if (i < neg_vars.size())
      define(neg_vars.at(i));
  }
}

bool PartonsPhaseSpaceGenerator::ktFactorised() const {
  return samplers_.at(0)->ktFactorised() && samplers_.at(1)->ktFactorised();
}

bool PartonsPhaseSpaceGenerator::generatePartonKinematics() {
  process().q1() = samplers_.at(0)->partonMomentum();
  process().q2() = samplers_.at(1)->partonMomentum();
  return true;
}

double PartonsPhaseSpaceGenerator::fluxes() const {
  return samplers_.at(0)->weight(process().x1(), process().mX2()) *
         samplers_.at(1)->weight(process().x2(), process().mY2());
}

const PartonSampler& PartonsPhaseSpaceGenerator::sampler(PartonSampler::Side side) const {
  return *samplers_.at(static_cast<size_t>(side));
}

proc::FactorisedProcess& PartonsPhaseSpaceGenerator::process() {
  if (!proc_)
    throw CG_FATAL("PartonsPhaseSpaceGenerator:process") << "Consumer process was not yet set.";
  return *proc_;
}

const proc::FactorisedProcess& PartonsPhaseSpaceGenerator::process() const {
  if (!proc_)
    throw CG_FATAL("PartonsPhaseSpaceGenerator:process") << "Consumer process was not yet set.";
  return *proc_;
}

ParametersDescription PartonsPhaseSpaceGenerator::description() {
  auto desc = PartonSampler::description();
  desc.setDescription("Two-parton phase space mapper");
  return desc;
}
