/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2024  Laurent Forthomme
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
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"

namespace cepgen {
  PartonsKTPhaseSpaceGenerator::PartonsKTPhaseSpaceGenerator(const ParametersList& params)
      : PartonsPhaseSpaceGenerator(params), log_parton_virtuality_(steer<bool>("logPartonVirtuality")) {}

  void PartonsKTPhaseSpaceGenerator::initialise() {
    const auto& kin = process().kinematics();

    // pick a parton flux parameterisation for each beam
    auto set_flux_properties = [](const Beam& beam, std::unique_ptr<PartonFlux>& flux) {
      auto params = beam.partonFluxParameters();
      const auto params_p_el = KTFluxFactory::get().describeParameters("BudnevElastic");
      const auto params_p_inel = KTFluxFactory::get().describeParameters("BudnevInelastic");
      const auto params_hi_el = KTFluxFactory::get().describeParameters("ElasticHeavyIon");
      if (params.name().empty()) {
        if (beam.elastic()) {
          if (HeavyIon::isHI(beam.integerPdgId()))
            params = params_hi_el.validate(params);
          else
            params = params_p_el.validate(params);
        } else
          params = params_p_inel.validate(params);
        //TODO: fermions/pions
      }
      flux = KTFluxFactory::get().build(params);
      if (!flux)
        throw CG_FATAL("PartonsKTPhaseSpaceGenerator:init")
            << "Failed to initiate a parton flux object with properties: " << params << ".";
      if (!flux->ktFactorised())
        throw CG_FATAL("PartonsKTPhaseSpaceGenerator:init")
            << "Invalid incoming parton flux modelling: " << flux->name() << ".";
    };
    set_flux_properties(kin.incomingBeams().positive(), pos_flux_);
    set_flux_properties(kin.incomingBeams().negative(), neg_flux_);

    // register the incoming partons' transverse virtualities range
    if (log_parton_virtuality_) {
      const auto log_lim_kt = kin.cuts().initial.qt.compute(std::log).truncate(Limits{-10., 10.});
      process()
          .defineVariable(
              m_qt1_, proc::Process::Mapping::exponential, log_lim_kt, "qt1", "Positive-z parton virtuality")
          .defineVariable(
              m_qt2_, proc::Process::Mapping::exponential, log_lim_kt, "qt2", "Negative-z parton virtuality");
    } else {
      const auto lim_kt = kin.cuts().initial.qt.truncate(Limits{1.e-5, 1.e3});
      process()
          .defineVariable(m_qt1_, proc::Process::Mapping::linear, lim_kt, "qt1", "Positive-z parton virtuality")
          .defineVariable(m_qt2_, proc::Process::Mapping::linear, lim_kt, "qt2", "Negative-z parton virtuality");
    }

    // register the incoming partons' azimuthal angles range
    const auto lim_phi = kin.cuts().initial.phi.truncate(Limits{0., 2. * M_PI});
    process()
        .defineVariable(
            m_phi_qt1_, proc::Process::Mapping::linear, lim_phi, "phi_qt1", "Positive-z parton azimuthal angle")
        .defineVariable(
            m_phi_qt2_, proc::Process::Mapping::linear, lim_phi, "phi_qt2", "Negative-z parton azimuthal angle");
  }

  bool PartonsKTPhaseSpaceGenerator::generatePartonKinematics() {
    // set the fully transverse kinematics (eta = 0) of initial partons
    process().q1() = Momentum::fromPtEtaPhiE(m_qt1_, 0., m_phi_qt1_);
    process().q2() = Momentum::fromPtEtaPhiE(m_qt2_, 0., m_phi_qt2_);
    return true;
  }

  double PartonsKTPhaseSpaceGenerator::fluxes() const {
    // factors 1/pi due to integration over d^2(kt1) d^2(kt2) instead of d(kt1^2) d(kt2^2)
    return (positiveFlux<KTFlux>().fluxMX2(process().x1(), m_qt1_ * m_qt1_, process().mX2()) * M_1_PI * m_qt1_) *
           (negativeFlux<KTFlux>().fluxMX2(process().x2(), m_qt2_ * m_qt2_, process().mY2()) * M_1_PI * m_qt2_);
  }

  ParametersDescription PartonsKTPhaseSpaceGenerator::description() {
    auto desc = PartonsPhaseSpaceGenerator::description();
    desc.setDescription("KT-dependent phase space mapper");
    desc.add<bool>("logPartonVirtuality", true);
    return desc;
  }
}  // namespace cepgen
