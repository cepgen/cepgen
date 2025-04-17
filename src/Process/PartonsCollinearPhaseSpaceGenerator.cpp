/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"

using namespace cepgen;

/// Collinear factorisation phase space generator
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Jul 2023
class PartonsCollinearPhaseSpaceGenerator final : public PartonsPhaseSpaceGenerator {
public:
  explicit PartonsCollinearPhaseSpaceGenerator(const ParametersList& params)
      : PartonsPhaseSpaceGenerator(params), log_parton_virtuality_(steer<bool>("logPartonVirtuality")) {}

  static ParametersDescription description() {
    auto desc = PartonsPhaseSpaceGenerator::description();
    desc.setDescription("Collinear phase space mapper");
    desc.add("logPartonVirtuality", true).setDescription("generate linearly to log(Q^2) instead of Q^2?");
    return desc;
  }

  bool ktFactorised() const override { return false; }
  bool generatePartonKinematics() override {
    //TODO: gaussian smearing of kt can be introduced here
    process().q1() = Momentum::fromPtYPhiM(0., 0., 0., std::sqrt(m_t1_));
    process().q2() = Momentum::fromPtYPhiM(0., 0., 0., std::sqrt(m_t2_));
    // define a window in central system invariant mass
    if (const auto invariant_mass = (process().q1() + process().q2()).mass();
        !process().kinematics().cuts().central.mass_sum.contains(invariant_mass))
      return false;
    return true;
  }
  double fluxes() const override {
    return positiveFlux<CollinearFlux>().fluxQ2(process().x1(), m_t1_) / m_t1_ *
           negativeFlux<CollinearFlux>().fluxQ2(process().x2(), m_t2_) / m_t2_;
  }

private:
  void initialise() override;

  const bool log_parton_virtuality_;
  // mapped variables
  double m_t1_{0.};  ///< Virtuality of the first intermediate parton
  double m_t2_{0.};  ///< Virtuality of the second intermediate parton
};

void PartonsCollinearPhaseSpaceGenerator::initialise() {
  const auto& kin = process().kinematics();

  // pick a parton flux parameterisation for each beam
  auto set_flux_properties = [&kin](const Beam& beam, std::unique_ptr<PartonFlux>& flux) {
    auto params = beam.partonFluxParameters();
    const auto params_e_el = CollinearFluxFactory::get().describeParameters(
        "EPAFlux",
        ParametersList().set(
            "formFactors",
            ParametersList().setName("PointLikeFermion").set<pdgid_t>("pdgId", std::abs(beam.integerPdgId()))));
    const auto params_p_el = CollinearFluxFactory::get().describeParameters(
        "EPAFlux", ParametersList().set("formFactors", beam.formFactors()));
    const auto params_p_inel = CollinearFluxFactory::get().describeParameters(
        "EPAFlux",
        ParametersList().set("formFactors",
                             ParametersList()
                                 .setName("InelasticNucleon")
                                 .set("structureFunctions", kin.incomingBeams().structureFunctions())));
    const auto params_hi_el = CollinearFluxFactory::get().describeParameters(
        "EPAFlux", ParametersList().set("formFactors", ParametersList().setName("HeavyIonDipole")));
    if (params.name().empty()) {
      if (beam.elastic()) {
        if (HeavyIon::isHI(beam.integerPdgId()))
          params = params_hi_el.validate(params);
        else if (const auto pdg_id = std::abs(beam.integerPdgId());
                 pdg_id == PDG::electron || pdg_id == PDG::muon || pdg_id == PDG::tau)
          params = params_e_el.validate(params);
        else
          params = params_p_el.validate(params);
      } else
        params = params_p_inel.validate(params);
      //TODO: fermions/pions
    }
    if (flux = CollinearFluxFactory::get().build(params); !flux)
      throw CG_FATAL("PartonsCollinearPhaseSpaceGenerator:init")
          << "Failed to initiate a parton flux object with properties: " << params << ".";
    if (flux->ktFactorised())
      throw CG_FATAL("PartonsCollinearPhaseSpaceGenerator:init")
          << "Invalid incoming parton flux modelling: " << flux->name() << ".";
  };
  set_flux_properties(kin.incomingBeams().positive(), pos_flux_);
  set_flux_properties(kin.incomingBeams().negative(), neg_flux_);

  // register the incoming partons' virtualities range
  const auto lim_q2_1 = kin.cuts().initial.q2.at(0).truncate(Limits{1.e-10, 5.}),
             lim_q2_2 = kin.cuts().initial.q2.at(1).truncate(Limits{1.e-10, 5.});
  if (log_parton_virtuality_)
    process()
        .defineVariable(
            m_t1_, proc::Process::Mapping::exponential, lim_q2_1.compute(std::log), "Positive-z parton virtuality")
        .defineVariable(
            m_t2_, proc::Process::Mapping::exponential, lim_q2_2.compute(std::log), "Negative-z parton virtuality");
  else
    process()
        .defineVariable(m_t1_, proc::Process::Mapping::linear, lim_q2_1, "Positive-z parton virtuality")
        .defineVariable(m_t2_, proc::Process::Mapping::linear, lim_q2_2, "Negative-z parton virtuality");
}
REGISTER_PARTONS_PHASE_SPACE_GENERATOR("coll", PartonsCollinearPhaseSpaceGenerator);
