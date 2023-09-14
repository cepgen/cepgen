/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2023  Laurent Forthomme
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
#include "CepGen/Process/KTPhaseSpaceGenerator.h"
#include "CepGen/Process/Process.h"

namespace cepgen {
  namespace proc {
    KTPhaseSpaceGenerator::KTPhaseSpaceGenerator(Process* proc) : PhaseSpaceGenerator(proc) {}

    void KTPhaseSpaceGenerator::init() {
      auto& kin = process().kinematics();
      auto set_flux_properties = [](Beam& beam, std::unique_ptr<PartonFlux>& flux) {
        auto params = beam.partonFluxParameters();
        if (params.name<std::string>().empty()) {
          if (beam.elastic()) {
            if (HeavyIon::isHI(beam.pdgId()))
              params = PartonFluxFactory::get().describeParameters("ElasticHeavyIonKT").validate(params);
            else
              params = PartonFluxFactory::get().describeParameters("BudnevElasticKT").validate(params);
          } else
            params = PartonFluxFactory::get().describeParameters("BudnevInelasticKT").validate(params);
          //TODO: fermions/pions
        }
        flux = std::move(PartonFluxFactory::get().build(params));
        if (!flux)
          throw CG_FATAL("KTPhaseSpaceGenerator:init")
              << "Failed to initiate a parton flux object with properties: " << params << ".";
      };
      set_flux_properties(kin.incomingBeams().positive(), pos_flux_);
      set_flux_properties(kin.incomingBeams().negative(), neg_flux_);

      if (!pos_flux_->ktFactorised() || !neg_flux_->ktFactorised())
        throw CG_FATAL("KTPhaseSpaceGenerator:init")
            << "Invalid incoming parton fluxes: " << std::vector<std::string>{pos_flux_->name(), neg_flux_->name()}
            << ".";

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      const auto log_lim_kt = kin.cuts().initial.qt.compute(std::log).truncate(Limits{-10., 10.});
      process().defineVariable(m_qt1_, Process::Mapping::exponential, log_lim_kt, "First incoming parton virtuality");
      process().defineVariable(m_qt2_, Process::Mapping::exponential, log_lim_kt, "Second incoming parton virtuality");

      const auto lim_phi_kt = kin.cuts().initial.phi_qt.truncate(Limits{0., 2. * M_PI});
      process().defineVariable(
          m_phi_qt1_, Process::Mapping::linear, lim_phi_kt, "First incoming parton azimuthal angle");
      process().defineVariable(
          m_phi_qt2_, Process::Mapping::linear, lim_phi_kt, "Second incoming parton azimuthal angle");
    }

    bool KTPhaseSpaceGenerator::generatePartonKinematics() {
      // set the transverse kinematics of initial partons
      process().q1() = Momentum::fromPtEtaPhiE(m_qt1_, 0., m_phi_qt1_);
      process().q2() = Momentum::fromPtEtaPhiE(m_qt2_, 0., m_phi_qt2_);
      return true;
    }

    double KTPhaseSpaceGenerator::fluxes() const {
      const auto& flux1 = dynamic_cast<const KTFlux&>(*pos_flux_);
      const auto& flux2 = dynamic_cast<const KTFlux&>(*neg_flux_);

      // factors 1/2pi and 1/2pi due to integration over d^2(kt1) d^2(kt2) instead of d(kt1^2) d(kt2^2)
      return (flux1.fluxMX2(process().x1(), m_qt1_ * m_qt1_, process().mX2()) * M_1_PI) *
             (flux2.fluxMX2(process().x2(), m_qt2_ * m_qt2_, process().mY2()) * M_1_PI);
    }
  }  // namespace proc
}  // namespace cepgen
