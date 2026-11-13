/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/PartonSampler.h"

using namespace cepgen;

//--------------------------------------------------------------------------------------------------
// generic single-beam sampler
//--------------------------------------------------------------------------------------------------

PartonSampler::PartonSampler(const ParametersList& params)
    : SteeredObject(params), log_parton_virtuality_(steer<bool>("logPartonVirtuality")) {}

void PartonSampler::initialise(const Kinematics& kin, Side side) {
  side_ = side;
  const auto& beam = side_ == Side::positive ? kin.incomingBeams().positive() : kin.incomingBeams().negative();
  auto params = beam.partonFluxParameters();
  if (params.name().empty())  // no flux steered by the user; pick the default one for this beam type
    params = defaultFluxParameters(beamCategory(beam), beam, kin).validate(params);
  if (flux_ = buildFlux(params); !flux_)
    throw CG_FATAL("PartonSampler:initialise")
        << "Failed to initiate a parton flux object with properties: " << params << ".";
  flux_->initialise(beam);
  if (flux_->ktFactorised() != ktFactorised())
    throw CG_FATAL("PartonSampler:initialise")
        << "Invalid incoming parton flux modelling for a " << (ktFactorised() ? "kT-factorised" : "collinear")
        << " sampler: " << flux_->name() << ".";
  fluxInitialised();
}

const PartonFlux& PartonSampler::flux() const {
  if (!flux_)
    throw CG_FATAL("PartonSampler:flux") << "Parton flux was not yet initialised for this sampler.";
  return *flux_;
}

PartonSampler::BeamCategory PartonSampler::beamCategory(const Beam& beam) {
  if (!beam.elastic())
    return BeamCategory::inelasticHadron;
  if (HeavyIon::isHI(beam.integerPdgId()))
    return BeamCategory::elasticHeavyIon;
  if (const auto pdg_id = std::abs(beam.integerPdgId());
      pdg_id == PDG::electron || pdg_id == PDG::muon || pdg_id == PDG::tau)
    return BeamCategory::elasticPointLike;
  return BeamCategory::elasticHadron;  //TODO: fermions/pions
}

std::string PartonSampler::variableName(const std::string& base) const {
  return base + (side_ == Side::positive ? "1" : "2");
}

std::string PartonSampler::sideLabel() const { return side_ == Side::positive ? "Positive-z" : "Negative-z"; }

ParametersDescription PartonSampler::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Single-beam incoming parton sampler");
  desc.add("logPartonVirtuality", true).setDescription("generate linearly to log(Q^2) instead of Q^2?");
  return desc;
}

//--------------------------------------------------------------------------------------------------
// kT-factorised sampler
//--------------------------------------------------------------------------------------------------

ParametersDescription KTPartonSampler::defaultFluxParameters(BeamCategory category,
                                                             const Beam&,
                                                             const Kinematics&) const {
  switch (category) {
    case BeamCategory::elasticPointLike:
      return KTFluxFactory::get().describeParameters("BudnevElasticLepton");
    case BeamCategory::elasticHadron:
      return KTFluxFactory::get().describeParameters("BudnevElastic");
    case BeamCategory::inelasticHadron:
      return KTFluxFactory::get().describeParameters("BudnevInelastic");
    case BeamCategory::elasticHeavyIon:
      return KTFluxFactory::get().describeParameters("ElasticHeavyIon");
  }
  throw CG_FATAL("KTPartonSampler:defaultFluxParameters") << "Invalid beam category.";
}

std::unique_ptr<PartonFlux> KTPartonSampler::buildFlux(const ParametersList& params) const {
  return KTFluxFactory::get().build(params);
}

void KTPartonSampler::fluxInitialised() { kt_flux_ = &flux<KTFlux>(); }

std::vector<IntegrationVariable> KTPartonSampler::integrationVariables(const cuts::Initial& cuts) {
  std::vector<IntegrationVariable> vars;
  // parton transverse virtuality
  if (log_parton_virtuality_)
    vars.push_back({qt_,
                    proc::Process::Mapping::exponential,
                    cuts.qt.compute(std::log).truncate(Limits{-10., 10.}),
                    variableName("qt"),
                    sideLabel() + " parton virtuality"});
  else
    vars.push_back({qt_,
                    proc::Process::Mapping::linear,
                    cuts.qt.truncate(Limits{1.e-5, 1.e3}),
                    variableName("qt"),
                    sideLabel() + " parton virtuality"});
  // parton azimuthal angle
  vars.push_back({phi_qt_,
                  proc::Process::Mapping::linear,
                  cuts.phi.truncate(Limits{0., 2. * M_PI}),
                  variableName("phi_qt"),
                  sideLabel() + " parton azimuthal angle"});
  return vars;
}

Momentum KTPartonSampler::partonMomentum() const {
  return Momentum::fromPtEtaPhiE(qt_, 0., phi_qt_);  // fully transverse kinematics (eta = 0)
}

double KTPartonSampler::weight(double x, double mx2) const {
  // factor 1/pi due to integration over d^2(kt) instead of d(kt^2)
  return kt_flux_->fluxMX2(x, qt_ * qt_, mx2) * M_1_PI * qt_;
}

//--------------------------------------------------------------------------------------------------
// collinear sampler
//--------------------------------------------------------------------------------------------------

ParametersDescription CollinearPartonSampler::defaultFluxParameters(BeamCategory category,
                                                                    const Beam& beam,
                                                                    const Kinematics& kin) const {
  auto epa_flux = [](const ParametersList& form_factors) {
    return CollinearFluxFactory::get().describeParameters("EPAFlux", ParametersList().set("formFactors", form_factors));
  };
  switch (category) {
    case BeamCategory::elasticPointLike:
      return epa_flux(
          ParametersList().setName("PointLikeFermion").set<pdgid_t>("pdgId", std::abs(beam.integerPdgId())));
    case BeamCategory::elasticHadron:
      return epa_flux(beam.formFactors());
    case BeamCategory::inelasticHadron:
      return epa_flux(ParametersList()
                          .setName("InelasticNucleon")
                          .set("structureFunctions", kin.incomingBeams().structureFunctions()));
    case BeamCategory::elasticHeavyIon:
      return epa_flux(ParametersList().setName("HeavyIonDipole"));
  }
  throw CG_FATAL("CollinearPartonSampler:defaultFluxParameters") << "Invalid beam category.";
}

std::unique_ptr<PartonFlux> CollinearPartonSampler::buildFlux(const ParametersList& params) const {
  return CollinearFluxFactory::get().build(params);
}

void CollinearPartonSampler::fluxInitialised() { coll_flux_ = &flux<CollinearFlux>(); }

std::vector<IntegrationVariable> CollinearPartonSampler::integrationVariables(const cuts::Initial& cuts) {
  const auto lim_q2 = cuts.q2.at(static_cast<size_t>(side_)).truncate(Limits{1.e-10, 5.});
  if (log_parton_virtuality_)
    return {{t_,
             proc::Process::Mapping::exponential,
             lim_q2.compute(std::log),
             variableName("t"),
             sideLabel() + " parton virtuality"}};
  return {{t_, proc::Process::Mapping::linear, lim_q2, variableName("t"), sideLabel() + " parton virtuality"}};
}

Momentum CollinearPartonSampler::partonMomentum() const {
  //TODO: gaussian smearing of kt can be introduced here
  return Momentum::fromPtYPhiM(0., 0., 0., std::sqrt(t_));
}

double CollinearPartonSampler::weight(double x, double) const { return coll_flux_->fluxQ2(x, t_) / t_; }
