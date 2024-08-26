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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/PartonFluxes/PartonFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGenGammaUPC/GammaUPCInterface.h"

using namespace cepgen;

namespace gammaUPC {
  class TrivialPartonFlux final : public PartonFlux, public SteeredObject<TrivialPartonFlux> {
  public:
    explicit TrivialPartonFlux(const ParametersList& params)
        : PartonFlux(params),
          SteeredObject<TrivialPartonFlux>(params),
          beam_id_(SteeredObject<TrivialPartonFlux>::steer<int>("beamPdgId")),
          parton_id_(SteeredObject<TrivialPartonFlux>::steer<int>("partonPdgId")),
          mass2_(std::pow(PDG::get().mass(beam_id_), 2)) {}

    static ParametersDescription description() {
      auto desc = cepgen::PartonFlux::description();
      desc.add<int>("beamPdgId", PDG::proton);
      desc.add<int>("partonPdgId", PDG::photon);
      return desc;
    }

    bool fragmenting() const override { return false; }
    pdgid_t partonPdgId() const override { return parton_id_; }
    double mass2() const override { return mass2_; }

  private:
    const spdgid_t beam_id_, parton_id_;
    const double mass2_;
  };

  class PartonsPhaseSpaceGenerator final : public cepgen::PartonsPhaseSpaceGenerator {
  public:
    explicit PartonsPhaseSpaceGenerator(const ParametersList& params)
        : cepgen::PartonsPhaseSpaceGenerator(params),
          log_parton_virtuality_(steer<bool>("logPartonVirtuality")),
          force_p_nonhadr_(steer<bool>("forcePNonHadr")),
          heavy_ion_mode_(steerAs<int, gammaUPC::HeavyIonMode>("heavyIonMode")) {}
    static ParametersDescription description() {
      auto desc = cepgen::PartonsPhaseSpaceGenerator::description();
      desc.setDescription("Two-photon central phase space generator w/ gammaUPC fluxes");
      desc.add<bool>("logPartonVirtuality", false);
      desc.add<bool>("forcePNonHadr", true)
          .setDescription("disallow the hadronisation of beam remnants after parton emission?");
      desc.addAs<int, HeavyIonMode>("heavyIonMode", HeavyIonMode::HardSphere)
          .setDescription("type of modelling for the heavy ion parton emission")
          .allow(static_cast<int>(gammaUPC::HeavyIonMode::HardSphere), "hard sphere model")
          .allow(static_cast<int>(gammaUPC::HeavyIonMode::WoodsSaxon), "Woods-Saxon model");
      return desc;
    }

    void initialise() override {
      pos_flux_.reset(new TrivialPartonFlux(params_));
      neg_flux_.reset(new TrivialPartonFlux(params_));
      // register the incoming partons' virtuality
      const auto lim_q2_1 = process().kinematics().cuts().initial.q2.at(0).truncate(Limits{1.e-10, 5.}),
                 lim_q2_2 = process().kinematics().cuts().initial.q2.at(1).truncate(Limits{1.e-10, 5.});
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
      const auto hi_a = HeavyIon::isHI(process().kinematics().incomingBeams().positive().integerPdgId()),
                 hi_b = HeavyIon::isHI(process().kinematics().incomingBeams().negative().integerPdgId());
      mass_number_prod_ = 1;
      if (hi_a) {  // determine beams mode, and mass numbers for HI cases
        mass_number_prod_ *= HeavyIon::fromPdgId(process().kinematics().incomingBeams().positive().integerPdgId()).A;
        beams_mode_ = hi_b ? BeamsMode::AB : BeamsMode::Ap;
      } else
        beams_mode_ = hi_b ? BeamsMode::pB : BeamsMode::pp;
      if (hi_b)
        mass_number_prod_ *= HeavyIon::fromPdgId(process().kinematics().incomingBeams().negative().integerPdgId()).A;
      to_collider_.ebeam[0] = process().kinematics().incomingBeams().positive().momentum().energy();
      to_collider_.ebeam[1] = process().kinematics().incomingBeams().negative().momentum().energy();
      if (beams_mode_ == BeamsMode::Ap)  // special case for asymmetric beams
        std::swap(to_collider_.ebeam[0], to_collider_.ebeam[1]);
    }
    bool ktFactorised() const override { return false; }
    bool generatePartonKinematics() override {
      static Limits x_limits{0., 1.};
      if (!x_limits.contains(process().x1(), false) || !x_limits.contains(process().x2(), false))
        return false;
      process().q1() = Momentum::fromPtYPhiM(0., 0., 0., std::sqrt(m_t1_));
      process().q2() = Momentum::fromPtYPhiM(0., 0., 0., std::sqrt(m_t2_));
      return true;
    }
    double fluxes() const override {
      const auto prefactor = mass_number_prod_ * process().x1() * process().x2();
      switch (beams_mode_) {
        case BeamsMode::pp:
          return prefactor * twoPhotonFluxPP(process().x1(), process().x2(), force_p_nonhadr_);
        case BeamsMode::pB:
          return prefactor * twoPhotonFluxPA(process().x1(), process().x2(), force_p_nonhadr_, heavy_ion_mode_);
        case BeamsMode::Ap:
          return prefactor * twoPhotonFluxPA(process().x2(), process().x1(), force_p_nonhadr_, heavy_ion_mode_);
        case BeamsMode::AB:
          return prefactor * twoPhotonFluxAB(process().x1(), process().x2(), force_p_nonhadr_, heavy_ion_mode_);
        default:
          throw CG_FATAL("gammaUPC:fluxes") << "Unknown beams mode: " << static_cast<int>(beams_mode_) << ".";
      }
    }

  private:
    const bool log_parton_virtuality_;
    const bool force_p_nonhadr_;
    const HeavyIonMode heavy_ion_mode_;

    enum struct BeamsMode { pp, pB, Ap, AB } beams_mode_{BeamsMode::pp};  ///< type of beam-beam system
    unsigned short mass_number_prod_{1};                                  ///< product of beams mass numbers

    // mapped variables in process kinematics definition
    double m_t1_{0.};  ///< virtuality of positive-z beam parton
    double m_t2_{0.};  ///< virtuality of negative-z beam parton
  };
}  // namespace gammaUPC
using gammaUPC_PPSG = gammaUPC::PartonsPhaseSpaceGenerator;
REGISTER_PARTONS_PHASE_SPACE_GENERATOR("gammaUPC", gammaUPC_PPSG);
