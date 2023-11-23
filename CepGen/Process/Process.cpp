/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace proc {
    Process::Process(const ParametersList& params)
        : NamedModule(params),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          rnd_gen_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {
      const auto& kin_params = steer<ParametersList>("kinematics");
      if (!kin_params.empty())
        kinematics().setParameters(kin_params);
      if (steer<bool>("hasEvent"))
        event_.reset(new Event);
    }

    Process::Process(const Process& proc)
        : NamedModule(proc),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          rnd_gen_(RandomGeneratorFactory::get().build(proc.rnd_gen_->parameters())),
          s_(proc.s_),
          sqs_(proc.sqs_),
          mA2_(proc.mA2_),
          mB2_(proc.mB2_),
          mapped_variables_(proc.mapped_variables_),
          point_coord_(proc.point_coord_),
          base_jacobian_(proc.base_jacobian_) {
      if (proc.event_)
        event_.reset(new Event(*proc.event_));
      CG_DEBUG("Process").log([&](auto& log) {
        log << "Process " << name_ << " cloned with "
            << utils::s("integration variable", mapped_variables_.size(), true) << ":";
        for (const auto& var : mapped_variables_)
          log << "\n\t" << var.index << ") " << var.description << " (type: " << var.type << ", limits: " << var.limits
              << ").";
        if (event_)
          log << "\n\t" << *event_;
      });
    }

    std::unique_ptr<Process> Process::clone() const {
      throw CG_FATAL("Process:clone") << "Process \"" << name_ << "\" has no cloning method implementation!";
    }

    Momentum& Process::pA() { return event().oneWithRole(Particle::IncomingBeam1).momentum(); }

    const Momentum& Process::pA() const { return event().oneWithRole(Particle::IncomingBeam1).momentum(); }

    Momentum& Process::pB() { return event().oneWithRole(Particle::IncomingBeam2).momentum(); }

    const Momentum& Process::pB() const { return event().oneWithRole(Particle::IncomingBeam2).momentum(); }

    Momentum& Process::pX() { return event().oneWithRole(Particle::OutgoingBeam1).momentum(); }

    const Momentum& Process::pX() const { return event().oneWithRole(Particle::OutgoingBeam1).momentum(); }

    Momentum& Process::pY() { return event().oneWithRole(Particle::OutgoingBeam2).momentum(); }

    const Momentum& Process::pY() const { return event().oneWithRole(Particle::OutgoingBeam2).momentum(); }

    Momentum& Process::q1() { return event().oneWithRole(Particle::Parton1).momentum(); }

    const Momentum& Process::q1() const { return event().oneWithRole(Particle::Parton1).momentum(); }

    Momentum& Process::q2() { return event().oneWithRole(Particle::Parton2).momentum(); }

    const Momentum& Process::q2() const { return event().oneWithRole(Particle::Parton2).momentum(); }

    Momentum& Process::pc(size_t i) {
      if (event()[Particle::CentralSystem].size() <= i)
        throw CG_FATAL("Process:pc") << "Trying to retrieve central particle #" << i << " while only "
                                     << event()[Particle::CentralSystem].size() << " is/are registered.";
      return event()[Particle::CentralSystem].at(i).get().momentum();
    }

    const Momentum& Process::pc(size_t i) const {
      if (event()(Particle::CentralSystem).size() <= i)
        throw CG_FATAL("Process:pc") << "Trying to retrieve central particle #" << i << " while only "
                                     << event()(Particle::CentralSystem).size() << " is/are registered.";
      return event()(Particle::CentralSystem).at(i).momentum();
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

    Process& Process::defineVariable(double& out, const Mapping& type, const Limits& lim, const std::string& descr) {
      if (!lim.valid())
        throw CG_FATAL("Process:defineVariable")
            << "The limits for '" << descr << "' (" << lim << ") could not be retrieved from the user configuration.";

      double jacob_weight = 1.;  // initialise the local weight for this variable
      switch (type) {
        case Mapping::linear:
          jacob_weight = lim.range();
          break;
        case Mapping::square:
          jacob_weight = 2. * lim.range();
          break;
        case Mapping::exponential:
          jacob_weight = lim.range();
          break;
        case Mapping::power_law:
          jacob_weight = log(lim.max() / lim.min());
          break;
      }
      const auto var_desc = descr.empty() ? utils::format("var%z", mapped_variables_.size()) : descr;
      mapped_variables_.emplace_back(MappingVariable{var_desc, lim, out, type, mapped_variables_.size()});
      point_coord_.emplace_back(0.);
      base_jacobian_ *= jacob_weight;
      CG_DEBUG("Process:defineVariable") << "\n\t" << descr << " has been mapped to variable "
                                         << mapped_variables_.size() << ".\n\t"
                                         << "Allowed range for integration: " << lim << ".\n\t"
                                         << "Variable integration mode: " << type << ".\n\t"
                                         << "Weight in the Jacobian: " << jacob_weight << ".";
      return *this;
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
        switch (var.type) {
          case Mapping::linear: {
            var.value = var.limits.x(xv);
            // jacobian *= 1
          } break;
          case Mapping::exponential: {
            var.value = std::exp(var.limits.x(xv));  // transform back to linear
            jacobian *= var.value;
          } break;
          case Mapping::square: {
            const auto val = var.limits.x(xv);
            var.value = val * val;  // transform to square
            jacobian *= val;
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
                                      << "Proc.-specific integrand * Jacobian (excl. global Jacobian) = "
                                      << (me_integrand * aux_jacobian) << ".";
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
                                     << ParametersDescription(kin_.parameters(false)) << ".";

      clear();  // also resets the "first run" flag

      // build the coupling objects
      const auto& alpha_em = steer<ParametersList>("alphaEM");
      if (!alpha_em.empty())
        alphaem_ = AlphaEMFactory::get().build(alpha_em);
      const auto& alpha_s = steer<ParametersList>("alphaS");
      if (!alpha_s.empty())
        alphas_ = AlphaSFactory::get().build(alpha_s);

      const auto& p1 = kin_.incomingBeams().positive().momentum();
      const auto& p2 = kin_.incomingBeams().negative().momentum();
      //--- define incoming system
      if (event_) {
        auto& ib1 = event_->oneWithRole(Particle::IncomingBeam1);
        ib1.setPdgId(kin_.incomingBeams().positive().pdgId());
        ib1.setMomentum(p1);
        auto& ib2 = event_->oneWithRole(Particle::IncomingBeam2);
        ib2.setPdgId(kin_.incomingBeams().negative().pdgId());
        ib2.setMomentum(p2);
        auto& ob1 = event_->oneWithRole(Particle::OutgoingBeam1);
        ob1.setPdgId(kin_.incomingBeams().positive().pdgId());
        ob1.setStatus(kin_.incomingBeams().positive().elastic() ? Particle::Status::FinalState
                                                                : Particle::Status::Unfragmented);
        auto& ob2 = event_->oneWithRole(Particle::OutgoingBeam2);
        ob2.setPdgId(kin_.incomingBeams().negative().pdgId());
        ob2.setStatus(kin_.incomingBeams().negative().elastic() ? Particle::Status::FinalState
                                                                : Particle::Status::Unfragmented);
        for (auto& cp : (*event_)[Particle::CentralSystem])
          cp.get().setPdgId(cp.get().pdgId());
      }
      s_ = kin_.incomingBeams().s();
      sqs_ = std::sqrt(s_);

      mA2_ = p1.mass2();
      mB2_ = p2.mass2();
      wcm_ = 0.5 * (1. + std::sqrt(1. - 4. * std::sqrt(mA2_ * mB2_) / s_));

      prepareKinematics();

      if (event_) {
        CG_DEBUG("Process:initialise") << "Kinematics successfully set!\n"
                                       << "  sqrt(s) = " << sqs_ * 1.e-3 << " TeV,\n"
                                       << "  p1=" << p1 << ",\tmass=" << p1.mass() << " GeV\n"
                                       << "  p2=" << p2 << ",\tmass=" << p2.mass() << " GeV.";
        clearEvent();
      }
    }

    double Process::alphaEM(double q) const {
      if (!alphaem_)
        throw CG_FATAL("Process:alphaEM")
            << "Trying to compute the electromagnetic running coupling while it is not initialised.";
      return (*alphaem_)(q);
    }

    double Process::alphaS(double q) const {
      if (!alphas_)
        throw CG_FATAL("Process:alphaS")
            << "Trying to compute the strong running coupling while it is not initialised.";
      return (*alphas_)(q);
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

    void Process::setEventContent(const std::unordered_map<Particle::Role, pdgids_t>& part_ids) {
      if (!event_)
        return;
      if (part_ids.count(Particle::Role::CentralSystem) == 0)
        throw CG_FATAL("Process") << "The central system was not specified for this process.";

      *event_ = Event::minimal(part_ids.at(Particle::Role::CentralSystem).size());
      for (const auto& role_vs_parts : part_ids) {
        auto evt_parts = (*event_)[role_vs_parts.first];
        if (evt_parts.size() != role_vs_parts.second.size())
          throw CG_FATAL("Process") << "Invalid number of '" << role_vs_parts.first << "' given. "
                                    << "Expecting " << evt_parts.size() << ", got " << role_vs_parts.second.size()
                                    << ".";
        for (size_t i = 0; i < evt_parts.size(); ++i) {
          auto& evt_part = evt_parts.at(i).get();
          const auto user_evt_part_pdgid = role_vs_parts.second.at(i);
          if (HeavyIon::isHI(user_evt_part_pdgid)) {
            evt_part.setPdgId(user_evt_part_pdgid);
            evt_part.momentum().setMass(HeavyIon::fromPdgId(user_evt_part_pdgid).mass());
          } else {
            const auto& part_info = PDG::get()(user_evt_part_pdgid);
            evt_part.setPdgId(user_evt_part_pdgid, part_info.charge / 3.);
            evt_part.momentum().setMass(part_info.mass);
          }
        }
      }
      event_->freeze();  // freeze the event as it is
    }

    ParametersDescription Process::description() {
      auto desc = ParametersDescription();
      desc.add<ParametersDescription>("alphaEM", AlphaEMFactory::get().describeParameters("fixed"))
          .setDescription("electromagnetic coupling evolution algorithm");
      desc.add<ParametersDescription>("alphaS", AlphaSFactory::get().describeParameters("pegasus"))
          .setDescription("strong coupling evolution algorithm");
      desc.add<bool>("hasEvent", true).setDescription("does the process carry an event definition");
      desc.add<ParametersDescription>("randomGenerator", ParametersDescription().setName<std::string>("stl"))
          .setDescription("random number generator engine");
      return desc;
    }

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
  }  // namespace proc
}  // namespace cepgen
