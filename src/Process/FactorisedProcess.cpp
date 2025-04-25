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
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace cepgen::proc;

FactorisedProcess::FactorisedProcess(const ParametersList& params)
    : Process(params),
      phase_space_generator_(PhaseSpaceGeneratorFactory::get().build(steer<ParametersList>("kinematicsGenerator"))),
      symmetrise_(steer<bool>("symmetrise")),
      store_alphas_(steer<bool>("storeAlphas")) {}

FactorisedProcess::FactorisedProcess(const ParametersList& params, const spdgids_t& central_particles)
    : FactorisedProcess(params) {
  setCentral(central_particles);
}

FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
    : Process(proc),
      phase_space_generator_(PhaseSpaceGeneratorFactory::get().build(proc.phase_space_generator_->parameters())),
      symmetrise_(proc.symmetrise_),
      store_alphas_(proc.store_alphas_) {
  setCentral(spdgids_t(proc.central_particles_.begin(), proc.central_particles_.end()));
}

void FactorisedProcess::setCentral(const spdgids_t& central_particles) {
  central_particles_ = std::vector<int>(central_particles.begin(), central_particles.end());
  phase_space_generator_->setCentral(central_particles_);
}

void FactorisedProcess::addEventContent() {
  CG_ASSERT(phase_space_generator_);
  const auto central_pdg_ids = phase_space_generator_->central();
  setEventContent({{Particle::Role::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::CentralSystem, spdgids_t(central_pdg_ids.begin(), central_pdg_ids.end())}});
}

void FactorisedProcess::prepareKinematics() {
  random_generator_ = RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"));
  if (!phase_space_generator_)
    throw CG_FATAL("FactorisedProcess:prepareKinematics")
        << "Phase space generator not set. Please check your process initialisation procedure, as you might "
           "be doing something irregular.";
  phase_space_generator_->initialise(this);

  event().oneWithRole(Particle::Role::Parton1).setIntegerPdgId(phase_space_generator_->partons().at(0));
  event().oneWithRole(Particle::Role::Parton2).setIntegerPdgId(phase_space_generator_->partons().at(1));
  event().map()[Particle::Role::CentralSystem].resize(central_particles_.size());

  CG_DEBUG("FactorisedProcess:prepareKinematics")
      << "Partons: " << phase_space_generator_->partons() << ", "
      << "central system: " << phase_space_generator_->central() << ". " << event();

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
  if (!phase_space_generator_->generate())
    return 0.;
  {  // compute and sanitise the momentum losses
    x1() = x2() = 0.;
    for (size_t i = 0; i < phase_space_generator_->central().size(); ++i) {
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
    return cent_weight * phase_space_generator_->weight() * kin_prefactor_;
  return 0.;
}

void FactorisedProcess::fillKinematics() {
  // parton systems
  event().oneWithRole(Particle::Role::Parton1).setMomentum(pA() - pX(), true);
  event().oneWithRole(Particle::Role::Parton2).setMomentum(pB() - pY(), true);

  if (symmetrise_ && random_generator_->uniformInt(0, 1) == 1) {  // symmetrise the el-in and in-el cases
    std::swap(pX(), pY());
    std::swap(q1(), q2());
    std::swap(pc(0), pc(1));
    for (auto* momentum : {&q1(), &q2(), &pX(), &pY(), &pc(0), &pc(1)})
      momentum->mirrorZ();
  }
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

double FactorisedProcess::that() const { return phase_space_generator_->that(); }

double FactorisedProcess::uhat() const { return phase_space_generator_->uhat(); }

bool FactorisedProcess::validatedBeamKinematics() {
  const auto invariant_mass = (q1() + q2()).mass();
  // define a window in central system invariant mass
  if (!kinematics().cuts().central.mass_sum.contains(invariant_mass))
    return 0.;

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
  desc.add("kinematicsGenerator", PhaseSpaceGeneratorFactory::get().describeParameters("kt:2to4"));
  desc.add("symmetrise", false).setDescription("Symmetrise along z the central system?");
  desc.add("storeAlphas", false).setDescription("store electromagnetic & strong coupling constants in event content?");
  desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
      .setDescription("random number generator engine");
  return desc;
}
