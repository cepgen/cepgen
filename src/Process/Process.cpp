/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::proc;

auto compute_value = [](double in, const Process::Mapping& type) -> double {
  switch (type) {
    case Process::Mapping::linear:
    case Process::Mapping::power_law:  //FIXME
    default:
      return in;
    case Process::Mapping::square:
      return in * in;
    case Process::Mapping::exponential:
      return std::exp(in);
  }
};

Process::Process(const ParametersList& params)
    : NamedModule(params),
      mp_(PDG::get().mass(PDG::proton)),
      mp2_(mp_ * mp_),
      rnd_gen_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {
  if (const auto& kin = steer<ParametersList>("kinematics"); !kin.empty())
    kin_.setParameters(kin);
  if (steer<bool>("hasEvent"))
    event_ = std::make_unique<Event>();
}

Process::Process(const Process& proc) : Process(proc.parameters()) { *this = proc; }

Process& Process::operator=(const Process& proc) {
  s_ = proc.s_;
  sqs_ = proc.sqs_;
  inv_s_ = proc.inv_s_;
  inv_sqs_ = proc.inv_sqs_;
  mA2_ = proc.mA2_;
  mB2_ = proc.mB2_;
  mX2_ = proc.mX2_;
  mY2_ = proc.mY2_;
  point_coord_ = proc.point_coord_;
  base_jacobian_ = proc.base_jacobian_;
  if (proc.event_)
    event_ = std::make_unique<Event>(*proc.event_);
  CG_DEBUG("Process").log([&](auto& log) {
    log << "Process " << name_ << " cloned with " << utils::s("integration variable", mapped_variables_.size(), true)
        << ":";
    for (const auto& var : mapped_variables_)
      log << "\n\t" << var.index << ") " << var.description << " (type: " << var.type << ", limits: " << var.limits
          << ").";
    if (event_)
      log << "\n\t" << *event_;
  });
  kin_ = proc.kin_;
  return *this;
}

std::unique_ptr<Process> Process::clone() const {
  throw CG_FATAL("Process:clone") << "Process \"" << name_ << "\" has no cloning method implementation!";
}

Momentum& Process::pA() { return event().oneWithRole(Particle::Role::IncomingBeam1).momentum(); }

const Momentum& Process::pA() const { return event().oneWithRole(Particle::Role::IncomingBeam1).momentum(); }

Momentum& Process::pB() { return event().oneWithRole(Particle::Role::IncomingBeam2).momentum(); }

const Momentum& Process::pB() const { return event().oneWithRole(Particle::Role::IncomingBeam2).momentum(); }

Momentum& Process::pX() { return event().oneWithRole(Particle::Role::OutgoingBeam1).momentum(); }

const Momentum& Process::pX() const { return event().oneWithRole(Particle::Role::OutgoingBeam1).momentum(); }

Momentum& Process::pY() { return event().oneWithRole(Particle::Role::OutgoingBeam2).momentum(); }

const Momentum& Process::pY() const { return event().oneWithRole(Particle::Role::OutgoingBeam2).momentum(); }

Momentum& Process::q1() { return event().oneWithRole(Particle::Role::Parton1).momentum(); }

const Momentum& Process::q1() const { return event().oneWithRole(Particle::Role::Parton1).momentum(); }

Momentum& Process::q2() { return event().oneWithRole(Particle::Role::Parton2).momentum(); }

const Momentum& Process::q2() const { return event().oneWithRole(Particle::Role::Parton2).momentum(); }

Momentum& Process::pc(size_t i) {
  if (event()[Particle::Role::CentralSystem].size() <= i)
    throw CG_FATAL("Process:pc") << "Trying to retrieve central particle #" << i << " while only "
                                 << event()[Particle::Role::CentralSystem].size() << " is/are registered.";
  return event()[Particle::Role::CentralSystem].at(i).get().momentum();
}

const Momentum& Process::pc(size_t i) const {
  if (event()(Particle::Role::CentralSystem).size() <= i)
    throw CG_FATAL("Process:pc") << "Trying to retrieve central particle #" << i << " while only "
                                 << event()(Particle::Role::CentralSystem).size() << " is/are registered.";
  return event()(Particle::Role::CentralSystem).at(i).momentum();
}

double Process::shat() const { return (q1() + q2()).mass2(); }

void Process::clear() {
  addEventContent();
  //--- initialise the "constant" (wrt x) part of the Jacobian
  base_jacobian_ = 1.;
  mapped_variables_.clear();
  CG_DEBUG("Process:clear") << "Process event content, and integration variables cleared.";
}

void Process::dumpVariables(std::ostream* os) const {
  std::ostringstream ss;
  ss << "List of variables handled by this process:";
  for (const auto& var : mapped_variables_)
    ss << "\n\t(" << var.index << ") " << var.type << " mapping (" << var.description << ")"
       << " in range " << var.limits;
  if (os)
    (*os) << ss.str();
  else
    CG_LOG << ss.str();
}

Process& Process::defineVariable(
    double& out, const Mapping& type, const Limits& lim, const std::string& name, const std::string& descr) {
  if (lim.min() == lim.max()) {
    if (lim.hasMin()) {
      out = compute_value(lim.min(), type);
      CG_DEBUG("Process:defineVariable") << "Quantity " << descr << " is set to be constant with a value " << out
                                         << ".";
      return *this;
    } else
      throw CG_FATAL("Process:defineVariable")
          << "The limits for '" << descr << "' (" << lim << ") could not be retrieved from the user configuration.";
  }

  double jacob_weight = 1.;  // initialise the local weight for this variable
  switch (type) {
    case Mapping::linear:
    case Mapping::exponential:
      jacob_weight = lim.range();
      break;
    case Mapping::square:
      jacob_weight = 2. * lim.range();
      break;
    case Mapping::power_law:
      jacob_weight = log(lim.max() / lim.min());
      break;
  }
  const auto var_desc =
      (!descr.empty() ? descr : (!name.empty() ? name : utils::format("var%z", mapped_variables_.size())));
  mapped_variables_.emplace_back(MappingVariable{name, var_desc, lim, out, type, mapped_variables_.size()});
  point_coord_.emplace_back(0.);
  base_jacobian_ *= jacob_weight;
  CG_DEBUG("Process:defineVariable") << "\n\t" << descr << " has been mapped to variable " << mapped_variables_.size()
                                     << ".\n\t"
                                     << "Allowed range for integration: " << lim << ".\n\t"
                                     << "Variable integration mode: " << type << ".\n\t"
                                     << "Weight in the Jacobian: " << jacob_weight << ".";
  return *this;
}

double Process::variableValue(size_t i, double x) const {
  const auto& var = mapped_variables_.at(i);
  return compute_value(var.limits.x(x), var.type);
}

double Process::generateVariables() const {
  if (mapped_variables_.size() == 0)
    throw CG_FATAL("Process:vars") << "No variables are mapped for this process!";
  if (base_jacobian_ == 0.)
    throw CG_FATAL("Process:vars") << "Point-independent component of the Jacobian for this "
                                   << "process is null.\n\t"
                                   << "Please check the validity of the phase space!";

  double jacobian = 1.;
  for (const auto& var : mapped_variables_) {
    if (!var.limits.valid())
      continue;
    if (var.index >= point_coord_.size())
      throw CG_FATAL("Process:x") << "Failed to retrieve coordinate " << var.index << " from "
                                  << "a dimension-" << ndim() << " process!";
    const auto& xv = point_coord_.at(var.index);  // between 0 and 1
    var.value = compute_value(var.limits.x(xv), var.type);
    switch (var.type) {
      case Mapping::linear: {
        // jacobian *= 1
      } break;
      case Mapping::exponential: {
        jacobian *= var.value;
      } break;
      case Mapping::square: {
        jacobian *= var.limits.x(xv);
      } break;
      case Mapping::power_law: {
        const double y = var.limits.max() / var.limits.min();
        var.value = var.limits.min() * std::pow(y, xv);
        jacobian *= var.value;
      } break;
    }
    CG_DEBUG_LOOP("Process:vars") << "\n\tvariable " << var.index << std::left << std::setw(60)
                                  << (!var.description.empty() ? " (" + var.description + ")" : "") << " in range "
                                  << std::setw(20) << var.limits << " has value " << std::setw(20) << var.value
                                  << " (x=" << point_coord_.at(var.index) << std::right << ")";
  }
  return jacobian;
}

double Process::weight(const std::vector<double>& x) {
  point_coord_ = x;

  //--- generate and initialise all variables and generate auxiliary
  //    (x-dependent) part of the Jacobian for this phase space point.
  const auto aux_jacobian = generateVariables();

  CG_DEBUG_LOOP("Process:weight").log([&](auto& log) {
    log << "Jacobian: " << base_jacobian_ << " * " << aux_jacobian << " = " << (base_jacobian_ * aux_jacobian)
        << ".\n\t";
    dumpPoint(&log.stream());
  });

  if (!utils::positive(aux_jacobian))
    return 0.;

  //--- compute the integrand
  const auto me_integrand = computeWeight();
  CG_DEBUG_LOOP("Process:weight") << "Integrand = " << me_integrand << "\n\t"
                                  << "(Process-specific integrand = " << me_integrand
                                  << ") * (Jacobian (excl. global Jacobian) = " << aux_jacobian
                                  << ") = " << (base_jacobian_ * aux_jacobian) << ".";
  if (!utils::positive(me_integrand))
    return 0.;

  //--- combine every component into a single weight for this point
  return (base_jacobian_ * aux_jacobian) * me_integrand * constants::GEVM2_TO_PB;
}

void Process::clearEvent() {
  if (event_)
    event_->restore();
}

const Event& Process::event() const {
  if (!event_)
    throw CG_FATAL("Process:event") << "Process does not have an event object!";
  return *event_;
}

Event& Process::event() {
  if (!event_)
    throw CG_FATAL("Process:event") << "Process does not have an event object!";
  return *event_;
}

Event* Process::eventPtr() {
  if (!event_)
    throw CG_FATAL("Process:event") << "Process does not have an event object!";
  return event_.get();
}

void Process::initialise() {
  CG_DEBUG("Process:initialise") << "Preparing to set the kinematics parameters. Input parameters: "
                                 << ParametersDescription(kin_.parameters()) << ".";

  clear();  // also resets the "first run" flag

  // build the coupling objects
  if (const auto& alpha_em = steer<ParametersList>("alphaEM"); !alpha_em.empty())
    alpha_em_ = AlphaEMFactory::get().build(alpha_em);
  if (const auto& alpha_s = steer<ParametersList>("alphaS"); !alpha_s.empty())
    alpha_s_ = AlphaSFactory::get().build(alpha_s);

  const auto& p1 = kin_.incomingBeams().positive().momentum();
  const auto& p2 = kin_.incomingBeams().negative().momentum();
  //--- define incoming system
  if (event_) {
    auto& ib1 = event_->oneWithRole(Particle::Role::IncomingBeam1);
    ib1.setIntegerPdgId(kin_.incomingBeams().positive().integerPdgId());
    ib1.setMomentum(p1);
    auto& ib2 = event_->oneWithRole(Particle::Role::IncomingBeam2);
    ib2.setIntegerPdgId(kin_.incomingBeams().negative().integerPdgId());
    ib2.setMomentum(p2);
    auto& ob1 = event_->oneWithRole(Particle::Role::OutgoingBeam1);
    ob1.setIntegerPdgId(kin_.incomingBeams().positive().integerPdgId());
    ob1.setStatus(kin_.incomingBeams().positive().elastic() ? Particle::Status::FinalState
                                                            : Particle::Status::Unfragmented);
    auto& ob2 = event_->oneWithRole(Particle::Role::OutgoingBeam2);
    ob2.setIntegerPdgId(kin_.incomingBeams().negative().integerPdgId());
    ob2.setStatus(kin_.incomingBeams().negative().elastic() ? Particle::Status::FinalState
                                                            : Particle::Status::Unfragmented);
    for (auto& cp : (*event_)[Particle::Role::CentralSystem])
      cp.get().setIntegerPdgId(cp.get().integerPdgId());
  }
  s_ = kin_.incomingBeams().s();
  sqs_ = std::sqrt(s_);
  inv_s_ = 1. / s_;
  inv_sqs_ = 1. / sqs_;

  mA2_ = mX2_ = p1.mass2();
  mB2_ = mY2_ = p2.mass2();
  wcm_ = 0.5 * (1. + std::sqrt(1. - 4. * std::sqrt(mA2_ * mB2_) * inv_s_));

  prepareKinematics();

  if (event_) {
    CG_DEBUG("Process:initialise").log([this, &p1, &p2](auto& log) {
      log << "Kinematics successfully set!\n"
          << "  sqrt(s) = " << sqs_ * 1.e-3 << " TeV,\n"
          << "  p1=" << p1 << ",\tmass=" << p1.mass() << " GeV\n"
          << "  p2=" << p2 << ",\tmass=" << p2.mass() << " GeV.\n";
      dumpVariables(&log.stream());
    });
    clearEvent();
  }
}

utils::RandomGenerator& Process::randomGenerator() const {
  if (!rnd_gen_)
    throw CG_FATAL("Process:randomGenerator") << "Process-local random generator was not yet initialised.";
  return *rnd_gen_;
}

double Process::alphaEM(double q) const {
  if (!alpha_em_)
    throw CG_FATAL("Process:alphaEM")
        << "Trying to compute the electromagnetic running coupling while it is not initialised.";
  return (*alpha_em_)(q);
}

double Process::alphaS(double q) const {
  if (!alpha_s_)
    throw CG_FATAL("Process:alphaS") << "Trying to compute the strong running coupling while it is not initialised.";
  return (*alpha_s_)(q);
}

void Process::dumpPoint(std::ostream* os) const {
  std::ostringstream oss;
  oss << "Number of integration parameters: " << mapped_variables_.size() << ", point: {"
      << utils::merge(point_coord_, ", ") << "}.";
  if (!os)
    CG_INFO("Process") << oss.str();
  else
    (*os) << oss.str();
}

void Process::setEventContent(const std::unordered_map<Particle::Role, spdgids_t>& part_ids) {
  if (!event_)
    return;
  if (part_ids.count(Particle::Role::CentralSystem) == 0)
    throw CG_FATAL("Process") << "The central system was not specified for this process.";

  *event_ = Event::minimal(part_ids.at(Particle::Role::CentralSystem).size());
  for (const auto& role_vs_parts : part_ids) {
    auto evt_parts = (*event_)[role_vs_parts.first];
    if (evt_parts.size() != role_vs_parts.second.size())
      throw CG_FATAL("Process") << "Invalid number of '" << role_vs_parts.first << "' given. "
                                << "Expecting " << evt_parts.size() << ", got " << role_vs_parts.second.size() << ".";
    for (size_t i = 0; i < evt_parts.size(); ++i) {
      auto& evt_part = evt_parts.at(i).get();
      const auto user_evt_part_pdgid = role_vs_parts.second.at(i);
      evt_part.setIntegerPdgId(user_evt_part_pdgid);
      evt_part.momentum().setMass(PDG::get().mass(user_evt_part_pdgid));
    }
  }
  event_->freeze();  // freeze the event as it is
}

void Process::setKinematics() {
  fillKinematics();
  if (event().hasRole(Particle::Role::Intermediate)) {
    Momentum intermediate_momentum;
    for (size_t i = 0; i < event()[Particle::Role::CentralSystem].size(); ++i)
      intermediate_momentum += pc(i);
    event().oneWithRole(Particle::Role::Intermediate).setMomentum(intermediate_momentum, true);
  }
}

ParametersDescription Process::description() {
  auto desc = ParametersDescription();
  desc.add("alphaEM", AlphaEMFactory::get().describeParameters("fixed"))
      .setDescription("e-m coupling evolution algorithm");
  desc.add("alphaS", AlphaSFactory::get().describeParameters("pegasus"))
      .setDescription("strong coupling evolution algorithm");
  desc.add("hasEvent", true).setDescription("does the process carry an event definition");
  desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
      .setDescription("random number generator engine");
  desc.add("kinematics", Kinematics::description());
  return desc;
}

namespace cepgen::proc {
  std::ostream& operator<<(std::ostream& os, const Process::Mapping& type) {
    switch (type) {
      case Process::Mapping::linear:
        return os << "linear";
      case Process::Mapping::exponential:
        return os << "exponential";
      case Process::Mapping::square:
        return os << "squared";
      case Process::Mapping::power_law:
        return os << "power law";
    }
    return os;
  }
}  // namespace cepgen::proc
