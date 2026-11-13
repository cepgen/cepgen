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
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/CentralPhaseSpaceFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::proc;
using namespace std::string_literals;

FactorisedProcess::FactorisedProcess(const ParametersList& params)
    : Process(params),
      symmetrise_(steer<bool>("symmetrise")),
      store_alphas_(steer<bool>("storeAlphas")),
      randomise_charge_(steer<bool>("randomiseCharge")),
      partons_generator_(partonsGeneratorType(params)),
      sampler_params_(partonSamplerParameters(params)) {}

FactorisedProcess::FactorisedProcess(const ParametersList& params, const spdgids_t& central_particles)
    : FactorisedProcess(params) {
  setCentral(central_particles);
}

FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
    : Process(proc),
      symmetrise_(proc.symmetrise_),
      store_alphas_(proc.store_alphas_),
      randomise_charge_(proc.randomise_charge_),
      partons_generator_(proc.partons_generator_),
      sampler_params_(proc.sampler_params_),
      central_particles_(proc.central_particles_) {}  // samplers and central mapper are built in prepareKinematics

std::string FactorisedProcess::partonsGeneratorType(const ParametersList& params) {
  auto type = params.has<std::string>("partonsGenerator") ? params.get<std::string>("partonsGenerator") : ""s;
  if (type.empty() && params.has<ParametersList>("kinematicsGenerator")) {  // deprecated steering
    const auto& kin_gen = params.get<ParametersList>("kinematicsGenerator");
    if (const auto tokens = utils::split(kin_gen.name(), ':'); tokens.size() >= 2)
      type = tokens.at(0);
    else if (kin_gen.has<std::string>("partonsGenerator"))
      type = kin_gen.get<std::string>("partonsGenerator");
    CG_WARNING("FactorisedProcess") << "Key 'kinematicsGenerator' is deprecated. Please use 'partonsGenerator = \""
                                    << (type.empty() ? "kt"s : type) << "\"' instead.";
  }
  if (type.empty())
    type = "kt";
  if (type != "kt" && type != "coll")
    throw CG_FATAL("FactorisedProcess") << "Invalid partons generator type: '" << type << "'. "
                                        << "Allowed values: 'kt' (kT-factorised), 'coll' (collinear).";
  return type;
}

ParametersList FactorisedProcess::partonSamplerParameters(const ParametersList& params) {
  auto log_parton_virtuality = params.has<bool>("logPartonVirtuality") ? params.get<bool>("logPartonVirtuality") : true;
  if (params.has<ParametersList>("kinematicsGenerator")) {  // deprecated steering
    if (const auto& kin_gen = params.get<ParametersList>("kinematicsGenerator");
        kin_gen.has<bool>("logPartonVirtuality"))
      log_parton_virtuality = kin_gen.get<bool>("logPartonVirtuality");
  }
  return ParametersList().set("logPartonVirtuality", log_parton_virtuality);
}

std::unique_ptr<PartonSampler> FactorisedProcess::buildPartonSampler(PartonSampler::Side side) const {
  const auto& beam = side == PartonSampler::Side::positive ? kinematics().incomingBeams().positive()
                                                           : kinematics().incomingBeams().negative();
  auto kt_factorised = partons_generator_ == "kt";
  // a parton flux explicitly steered for this beam fixes the sampler type
  if (const auto& flux_name = beam.partonFluxParameters().name(); !flux_name.empty()) {
    const auto is_kt_flux = KTFluxFactory::get().has(flux_name),
               is_coll_flux = CollinearFluxFactory::get().has(flux_name);
    if (is_kt_flux != is_coll_flux && is_kt_flux != kt_factorised) {
      CG_INFO("FactorisedProcess") << (side == PartonSampler::Side::positive ? "Positive" : "Negative")
                                   << "-z beam parton flux '" << flux_name << "' is "
                                   << (is_kt_flux ? "kT-factorised" : "collinear") << ": using a "
                                   << (is_kt_flux ? "kT-factorised" : "collinear")
                                   << " parton sampler for this beam (partons generator: '" << partons_generator_
                                   << "').";
      kt_factorised = is_kt_flux;
    }
  }
  if (kt_factorised)
    return std::make_unique<KTPartonSampler>(sampler_params_);
  return std::make_unique<CollinearPartonSampler>(sampler_params_);
}

std::unique_ptr<CentralPhaseSpace> FactorisedProcess::buildCentralPhaseSpace() const {
  // central system phase space mappers are indexed by the multiplicity they handle: "2to<N>"
  const auto mapper_name = "2to" + std::to_string(central_particles_.size());
  if (!CentralPhaseSpaceFactory::get().has(mapper_name))
    throw CG_FATAL("FactorisedProcess:buildCentralPhaseSpace")
        << "No central system phase space mapper available for a " << central_particles_.size()
        << "-particle central system (" << central_particles_ << "). "
        << "Registered mappers: " << CentralPhaseSpaceFactory::get().modules() << ".";
  return CentralPhaseSpaceFactory::get().build(mapper_name);
}

spdgids_t FactorisedProcess::partons() const {
  return {partonSampler(PartonSampler::Side::positive).flux().partonPdgId(),
          partonSampler(PartonSampler::Side::negative).flux().partonPdgId()};
}

const PartonSampler& FactorisedProcess::partonSampler(PartonSampler::Side side) const {
  const auto& sampler = parton_samplers_.at(static_cast<size_t>(side));
  if (!sampler)
    throw CG_FATAL("FactorisedProcess:partonSampler")
        << "Parton samplers are not yet built. They are only available once the process kinematics is prepared.";
  return *sampler;
}

bool FactorisedProcess::partonsKtFactorised() const {
  return partonSampler(PartonSampler::Side::positive).ktFactorised() &&
         partonSampler(PartonSampler::Side::negative).ktFactorised();
}

bool FactorisedProcess::generatePhaseSpace() {
  q1() = parton_samplers_[0]->partonMomentum();
  q2() = parton_samplers_[1]->partonMomentum();
  central_weight_ = central_phase_space_->generate((q1() + q2()).transverse(), central_momenta_);
  if (!utils::positive(central_weight_))
    return false;
  for (size_t i = 0; i < central_momenta_.size(); ++i)
    pc(i) = central_momenta_[i];
  return true;
}

double FactorisedProcess::phaseSpaceWeight() const {
  const auto fluxes_weight = parton_samplers_[0]->weight(x1(), mX2()) * parton_samplers_[1]->weight(x2(), mY2());
  if (!utils::positive(fluxes_weight))
    return 0.;
  return fluxes_weight * central_weight_;
}

void FactorisedProcess::setCentral(const spdgids_t& central_particles) {
  central_particles_ = central_particles;
  CG_DEBUG("FactorisedProcess:setCentral") << "Setting the central particles system to: " << central_particles_ << ".";
}

void FactorisedProcess::addEventContent() {
  setEventContent({{Particle::Role::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::CentralSystem, central_particles_}});
}

void FactorisedProcess::prepareKinematics() {
  random_generator_ = RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"));
  if (central_particles_.empty())
    CG_WARNING("FactorisedProcess:prepareKinematics") << "Empty central particles content.";

  // incoming partons: one single-beam sampler (flux and its own variables) per beam
  parton_samplers_ = {buildPartonSampler(PartonSampler::Side::positive),
                      buildPartonSampler(PartonSampler::Side::negative)};
  parton_samplers_[0]->initialise(kinematics(), PartonSampler::Side::positive);
  parton_samplers_[1]->initialise(kinematics(), PartonSampler::Side::negative);
  {  // register the samplers' variables, interleaved between beams (var1(+z), var1(-z), var2(+z), var2(-z), ...)
    const auto pos_vars = parton_samplers_[0]->integrationVariables(kinematics().cuts().initial),
               neg_vars = parton_samplers_[1]->integrationVariables(kinematics().cuts().initial);
    for (size_t i = 0; i < std::max(pos_vars.size(), neg_vars.size()); ++i) {
      if (i < pos_vars.size())
        cepgen::defineVariable(*this, pos_vars.at(i));
      if (i < neg_vars.size())
        cepgen::defineVariable(*this, neg_vars.at(i));
    }
  }
  // central system: phase space mapper chosen from the central system multiplicity
  if (central_phase_space_ = buildCentralPhaseSpace(); !central_phase_space_)
    throw CG_FATAL("FactorisedProcess:prepareKinematics") << "Failed to build the central system phase space mapper.";
  central_phase_space_->initialise(kinematics(), central_particles_);
  cepgen::defineVariables(*this, central_phase_space_->integrationVariables());

  event().oneWithRole(Particle::Role::Parton1).setIntegerPdgId(partons().at(0));
  event().oneWithRole(Particle::Role::Parton2).setIntegerPdgId(partons().at(1));
  event().map()[Particle::Role::CentralSystem].resize(central_particles_.size());

  CG_DEBUG("FactorisedProcess:prepareKinematics")
      << "Partons: " << partons() << " (" << (parton_samplers_[0]->ktFactorised() ? "kT" : "collinear") << " / "
      << (parton_samplers_[1]->ktFactorised() ? "kT" : "collinear") << "), "
      << "central system: " << central_particles_ << ". " << event();

  // register all process-dependent variables
  prepareFactorisedPhaseSpace();

  // register the outgoing remnants' variables
  if (!kinematics().incomingBeams().positive().elastic())
    defineVariable(mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive-z beam remnant squared mass");
  if (!kinematics().incomingBeams().negative().elastic())
    defineVariable(mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative-z beam remnant squared mass");
  // symmetrisation of the phase space: factor 2 in Jacobian for single-dissociation
  kin_prefactor_ = symmetrise_ && (kinematics().incomingBeams().mode() == mode::Kinematics::ElasticInelastic ||
                                   kinematics().incomingBeams().mode() == mode::Kinematics::InelasticElastic)
                       ? 2.
                       : 1.;
}

double FactorisedProcess::computeWeight() {
  if (!generatePhaseSpace())
    return 0.;
  {  // compute and sanitise the momentum losses
    x1() = x2() = 0.;
    for (size_t i = 0; i < central_particles_.size(); ++i) {
      const auto energy = pc(i).energy(), pz = pc(i).pz();
      x1() += std::fabs(energy + pz);
      x2() += std::fabs(energy - pz);
    }
    x1() *= inverseSqrtS(), x2() *= inverseSqrtS();
    if (!x_validity_range_.contains(x1()) || !x_validity_range_.contains(x2()))
      return 0.;
  }
  computeBeamKinematics();
  if (!validatedBeamKinematics())
    return 0.;
  if (const auto cent_weight = computeFactorisedMatrixElement(); utils::positive(cent_weight))
    return cent_weight * phaseSpaceWeight() * kin_prefactor_;
  return 0.;
}

void FactorisedProcess::fillKinematics() {
  {  // central system
    auto central_system = event()[Particle::Role::CentralSystem];
    if (randomise_charge_ && central_system.size() == 2) {  // randomise the charge of outgoing system
      const auto sign = random_generator_->uniformInt(0, 1) == 1;
      central_system[0].get().setAntiparticle(sign);
      central_system[1].get().setAntiparticle(!sign);
    }
    for (auto& central_particle : central_system)
      central_particle.get().setStatus(Particle::Status::FinalState);
  }
  if (symmetrise_ && random_generator_->uniformInt(0, 1) == 1) {  // symmetrise the el-in and in-el cases
    std::swap(pX(), pY());
    std::swap(q1(), q2());
    std::swap(pc(0), pc(1));
    for (auto* momentum : {&q1(), &q2(), &pX(), &pY()})
      momentum->mirrorZ();
    for (size_t i = 0; i < central_particles_.size(); ++i)  // loop over all central particles
      pc(i).mirrorZ();
  }
  // parton systems
  event().oneWithRole(Particle::Role::Parton1).setMomentum(pA() - pX(), true);
  event().oneWithRole(Particle::Role::Parton2).setMomentum(pB() - pY(), true);
  if (store_alphas_) {  // add couplings to metadata
    const auto two_parton_mass = (q1() + q2()).mass();
    event().metadata["alphaEM"] = alphaEM(two_parton_mass);
    event().metadata["alphaS"] = alphaS(two_parton_mass);
  }
}

void FactorisedProcess::computeBeamKinematics() {  // compute the four-momenta of the outgoing protons (or remnants)
  const double norm = wCM() * wCM() * s(), inv_norm = 1. / norm, prefactor = 0.5 * std::sqrt(norm);

  // positive-z outgoing beam system
  const auto px_p = (1. - x1()) * pA().p(), px_m = (mX2() + q1().p2()) * 0.25 / px_p;
  pX() = (-q1()).setPz(px_p - px_m).setEnergy(px_p + px_m);
  if (!kinematics().incomingBeams().positive().elastic())
    pX().setMass2(mX2());
  const double tau1 = q1().p2() / x1() * inv_norm;
  q1().setPz(+prefactor * (x1() - tau1)).setEnergy(+prefactor * (x1() + tau1));

  // negative-z outgoing beam system
  const auto py_m = (1. - x2()) * pB().p(), py_p = (mY2() + q2().p2()) * 0.25 / py_m;
  pY() = (-q2()).setPz(py_p - py_m).setEnergy(py_p + py_m);
  if (!kinematics().incomingBeams().negative().elastic())
    pY().setMass2(mY2());
  const double tau2 = q2().p2() / x2() * inv_norm;
  q2().setPz(-prefactor * (x2() - tau2)).setEnergy(+prefactor * (x2() + tau2));

  CG_DEBUG_LOOP("FactorisedProcess:validatedBeamKinematics")
      << "px+ = " << px_p << " / px- = " << px_m << ", py+ = " << py_p << " / py- = " << py_m << ". Remnants:\n\t"
      << "- X (positive-z): pX = " << pX() << ", mX = " << pX().mass() << "\n\t"
      << "- Y (negative-z): pY = " << pY() << ", mY = " << pY().mass() << ".";
}

//----- utilities

utils::RandomGenerator& FactorisedProcess::randomGenerator() const {
  if (!random_generator_)
    throw CG_FATAL("FactorisedProcess:randomGenerator") << "Process-local random generator was not yet initialised.";
  return *random_generator_;
}

double FactorisedProcess::that() const {
  if (central_particles_.size() != 2)
    throw CG_FATAL("FactorisedProcess:that") << "t-hat is only defined for a two-particle central system.";
  return 0.5 * ((q1() - pc(0)).mass2() + (q2() - pc(1)).mass2());
}

double FactorisedProcess::uhat() const {
  if (central_particles_.size() != 2)
    throw CG_FATAL("FactorisedProcess:uhat") << "u-hat is only defined for a two-particle central system.";
  return 0.5 * ((q1() - pc(1)).mass2() + (q2() - pc(0)).mass2());
}

bool FactorisedProcess::validatedBeamKinematics() {
  const auto invariant_mass = (q1() + q2()).mass();
  // define a window in central system invariant mass
  if (!kinematics().cuts().central.mass_sum.contains(invariant_mass))
    return false;

  // impose additional conditions for energy-momentum conservation
  if (!kinematics().incomingBeams().positive().elastic() && x2() * s() - invariant_mass - q2().p2() <= mX2())
    return false;
  if (!kinematics().incomingBeams().negative().elastic() && x1() * s() - invariant_mass - q1().p2() <= mY2())
    return false;

  // impose numerical conditions on X/Y 4-momenta
  if (std::fabs(pX().mass2() - mX2()) > NUM_LIMITS) {
    CG_WARNING("FactorisedProcess:validatedBeamKinematics")
        << "Invalid X system squared mass: " << pX().mass2() << "/" << mX2() << ".";
    return false;
  }
  if (std::fabs(pY().mass2() - mY2()) > NUM_LIMITS) {
    CG_WARNING("FactorisedProcess:validatedBeamKinematics")
        << "Invalid Y system squared mass: " << pY().mass2() << "/" << mY2() << ".";
    return false;
  }
  CG_DEBUG_LOOP("FactorisedProcess:validatedBeamKinematics")
      << "Squared c.m. energy = " << s() << " GeV^2\n\t"
      << "pos-z parton: " << q1() << ", m^2=" << q1().mass2() << " GeV, x1=" << x1() << ", p=" << q1().p() << "GeV\n\t"
      << "neg-z parton: " << q2() << ", m^2=" << q2().mass2() << " GeV, x2=" << x2() << ", p=" << q2().p() << "GeV.";
  return true;
}

ParametersDescription FactorisedProcess::description() {
  auto desc = Process::description();
  desc.setDescription("Unnamed factorised process");
  desc.add("kinematics",
           Kinematics::description()
               .addParametersDescriptionVector(
                   "partonFluxes",
                   PartonFluxFactory::get().describeParameters("BudnevElastic"),
                   std::vector(2, PartonFluxFactory::get().describeParameters("BudnevElastic").parameters()))
               .setDescription("Parton fluxes modelling"));
  desc.add("partonsGenerator", ""s)
      .setDescription("incoming partons generator type: 'kt' (kT-factorised) or 'coll' (collinear); empty: 'kt'");
  desc.add("logPartonVirtuality", true).setDescription("generate linearly to log(Q^2) instead of Q^2?");
  desc.add("symmetrise", false).setDescription("Symmetrise along z the central system?");
  desc.add("randomiseCharge", true).setDescription("randomise the charges of the central system (if charged)?");
  desc.add("storeAlphas", false).setDescription("store electromagnetic & strong coupling constants in event content?");
  desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
      .setDescription("random number generator engine");
  return desc;
}
