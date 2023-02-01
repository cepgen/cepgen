/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace proc {
    Process::Process(const ParametersList& params)
        : NamedModule(params), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {
      if (steer<bool>("hasEvent"))
        event_.reset(new Event);
    }

    Process::Process(const Process& proc)
        : NamedModule(proc),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          mapped_variables_(proc.mapped_variables_),
          point_coord_(proc.point_coord_),
          base_jacobian_(proc.base_jacobian_),
          s_(proc.s_),
          sqs_(proc.sqs_),
          mA2_(proc.mA2_),
          mB2_(proc.mB2_) {
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

    Momentum& Process::pc(size_t i) { return event()[Particle::CentralSystem].at(i).get().momentum(); }

    const Momentum& Process::pc(size_t i) const { return event()(Particle::CentralSystem).at(i).momentum(); }

    double Process::shat() const { return (q1() + q2()).mass2(); }

    void Process::clear() {
      addEventContent();
      //--- initialise the "constant" (wrt x) part of the Jacobian
      base_jacobian_ = 1.;
      mapped_variables_.clear();
      CG_DEBUG("Process:clear") << "Process event content, and integration variables cleared.";
    }

    void Process::dumpVariables(std::ostream* os) const {
      if (!os)
        os = &CG_LOG.stream();
      (*os) << "List of variables handled by this process:";
      for (const auto& var : mapped_variables_)
        (*os) << "\n\t(" << var.index << ") " << var.type << " mapping (" << var.description << ")"
              << " in range " << var.limits;
    }

    Process& Process::defineVariable(
        double& out, const Mapping& type, Limits in, const Limits& default_limits, const std::string& descr) {
      if (!in.valid()) {
        CG_DEBUG("Process:defineVariable") << descr << " could not be retrieved from the user configuration. "
                                           << "Setting it to the default value: " << default_limits << ".";
        in = default_limits;
      }

      Limits lim = in;
      out = 0.;                  // reset the variable
      double jacob_weight = 1.;  // initialise the local weight for this variable

      switch (type) {
        case Mapping::linear:
          jacob_weight = lim.range();
          break;
        case Mapping::square:
          jacob_weight = 2. * lim.range();
          break;
        case Mapping::exponential: {
          lim = {// limits already stored as log(limits)
                 (!lim.hasMin() || lim.min() == 0.) ? -10. : std::max(log(lim.min()), -10.),
                 (!lim.hasMax() || lim.max() == 0.) ? +10. : std::min(log(lim.max()), +10.)};
          jacob_weight = lim.range();  // use the linear version
        } break;
        case Mapping::power_law:
          jacob_weight = log(lim.max() / lim.min());
          break;
      }
      mapped_variables_.emplace_back(
          MappingVariable{descr.empty() ? utils::format("var%z", mapped_variables_.size()) : descr,
                          lim,
                          out,
                          type,
                          (unsigned short)mapped_variables_.size()});
      point_coord_.emplace_back(0.);
      base_jacobian_ *= jacob_weight;
      CG_DEBUG("Process:defineVariable") << "\n\t" << descr << " has been mapped to variable "
                                         << mapped_variables_.size() << ".\n\t"
                                         << "Allowed range for integration: " << in << " (" << lim << ").\n\t"
                                         << "Variable integration mode: " << type << ".\n\t"
                                         << "Weight in the Jacobian: " << jacob_weight << ".";
      return *this;
    }

    void Process::generateVariables() const {
      if (mapped_variables_.size() == 0)
        throw CG_FATAL("Process:vars") << "No variables are mapped for this process!";
      if (base_jacobian_ == 0.)
        throw CG_FATAL("Process:vars") << "Point-independent component of the Jacobian for this "
                                       << "process is null.\n\t"
                                       << "Please check the validity of the phase space!";

      for (const auto& var : mapped_variables_) {
        if (!var.limits.valid())
          continue;
        if (var.index >= point_coord_.size())
          throw CG_FATAL("Process:x") << "Failed to retrieve coordinate " << var.index << " from "
                                      << "a dimension-" << ndim() << " process!";
        const double xv = point_coord_.at(var.index);  // between 0 and 1
        switch (var.type) {
          case Mapping::linear: {
            var.value = var.limits.x(xv);
          } break;
          case Mapping::exponential: {          // limits already logarithmic
            var.value = exp(var.limits.x(xv));  // transform back to linear
          } break;
          case Mapping::square: {
            var.value = pow(var.limits.x(xv), 2);
          } break;
          case Mapping::power_law: {
            const double y = var.limits.max() / var.limits.min();
            var.value = var.limits.min() * pow(y, xv);
          } break;
        }
      }
      CG_DEBUG_LOOP("Process:vars").log([&](auto& dbg) {
        dbg << "Dump of all variables values:";
        for (const auto& var : mapped_variables_) {
          double value = 0.;
          switch (var.type) {
            case Mapping::linear:
            case Mapping::exponential:
            case Mapping::power_law:
              value = var.value;
              break;
            case Mapping::square:
              value = sqrt(var.value);
              break;
          }
          dbg << "\n\tvariable " << var.index << std::left << std::setw(60)
              << (!var.description.empty() ? " (" + var.description + ")" : "") << " in range " << std::setw(20)
              << var.limits << " has value " << std::setw(20) << value << " (x=" << point_coord_.at(var.index)
              << std::right << ")";
        }
      });
    }

    double Process::jacobian() const {
      double jac = 1.;
      for (const auto& var : mapped_variables_) {
        if (!var.limits.valid())
          continue;
        switch (var.type) {
          case Mapping::linear:
            break;
          case Mapping::square:
            jac *= sqrt(var.value);
            break;
          case Mapping::exponential:
          case Mapping::power_law:
            jac *= var.value;
            break;
        }
      }
      return jac;
    }

    double Process::weight(const std::vector<double>& x) {
      point_coord_ = x;
      if (CG_LOG_MATCH("Process:dumpPoint", debugInsideLoop))
        dumpPoint();

      //--- generate and initialise all variables
      generateVariables();

      //--- compute the integrand
      const auto me_integrand = computeWeight();
      if (me_integrand <= 0.)
        return 0.;

      //--- generate auxiliary (x-dependent) part of the Jacobian for
      //    this phase space point.
      const double aux_jacobian = jacobian();
      if (aux_jacobian <= 0.)
        return 0.;

      //--- combine every component into a single weight for this point
      const auto weight = (base_jacobian_ * aux_jacobian) * me_integrand;

      CG_DEBUG_LOOP("Process:weight") << "Jacobian: " << base_jacobian_ << " * " << aux_jacobian << " = "
                                      << (base_jacobian_ * aux_jacobian) << ".\n\t"
                                      << "Integrand = " << me_integrand << "\n\t"
                                      << "Proc.-specific integrand * Jacobian (excl. global Jacobian) = "
                                      << (me_integrand * aux_jacobian) << "\n\t"
                                      << "Point weight = " << weight << ".";

      return weight;
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

      kin_.incomingBeams().initialise();

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
        ob1.setStatus(kin_.incomingBeams().positive().fragmented() ? Particle::Status::Unfragmented
                                                                   : Particle::Status::FinalState);
        auto& ob2 = event_->oneWithRole(Particle::OutgoingBeam2);
        ob2.setPdgId(kin_.incomingBeams().negative().pdgId());
        ob2.setStatus(kin_.incomingBeams().negative().fragmented() ? Particle::Status::Unfragmented
                                                                   : Particle::Status::FinalState);
        for (auto& cp : (*event_)[Particle::CentralSystem])
          cp.get().setPdgId(cp.get().pdgId());
      }
      s_ = kin_.incomingBeams().s();
      sqs_ = std::sqrt(s_);

      mA2_ = p1.mass2();
      mB2_ = p2.mass2();

      prepareKinematics();

      if (event_) {
        CG_DEBUG("Process:initialise") << "Kinematics successfully set!\n"
                                       << "  sqrt(s) = " << sqs_ * 1.e-3 << " TeV,\n"
                                       << "  p1=" << p1 << ",\tmass=" << p1.mass() << " GeV\n"
                                       << "  p2=" << p2 << ",\tmass=" << p2.mass() << " GeV.";
        clearEvent();
      }
    }

    void Process::dumpPoint() const {
      CG_INFO("Process").log([&](auto& info) {
        info << "Number of integration parameters: " << mapped_variables_.size();
        for (unsigned short i = 0; i < point_coord_.size(); ++i)
          info << utils::format("\n\t  x(%2d) = %8.6f", i, point_coord_[i]);
        info << ".";
      });
    }

    void Process::setEventContent(const IncomingState& ini, const OutgoingState& fin) {
      if (!event_)
        return;

      event_->clear();
      //----- add the particles in the event

      //--- incoming state
      for (const auto& ip : ini) {
        auto& p = event_->addParticle(ip.first).get();
        const auto& part_info = PDG::get()(ip.second);
        p.setPdgId(ip.second, part_info.charge / 3.);
        p.setMass(part_info.mass);
        if (ip.first == Particle::IncomingBeam1 || ip.first == Particle::IncomingBeam2)
          p.setStatus(Particle::Status::PrimordialIncoming);
        if (ip.first == Particle::Parton1 || ip.first == Particle::Parton2)
          p.setStatus(Particle::Status::Incoming);
      }
      //--- central system (if not already there)
      const auto& central_system = ini.find(Particle::CentralSystem);
      if (central_system == ini.end()) {
        auto& p = event_->addParticle(Particle::Intermediate).get();
        p.setPdgId((pdgid_t)PDG::invalid);
        p.setStatus(Particle::Status::Propagator);
      }
      //--- outgoing state
      for (const auto& opl : fin) {  // pair(role, list of PDGids)
        for (const auto& pdg : opl.second) {
          auto& p = event_->addParticle(opl.first).get();
          const auto& part_info = PDG::get()(pdg);
          p.setPdgId(pdg, part_info.charge / 3.);
          p.setMass(part_info.mass);
        }
      }

      //----- define the particles parentage

      const Particles parts = event_->particles();
      for (const auto& p : parts) {
        Particle& part = (*event_)[p.id()];
        switch (part.role()) {
          case Particle::OutgoingBeam1:
          case Particle::Parton1:
            part.addMother(event_->oneWithRole(Particle::IncomingBeam1));
            break;
          case Particle::OutgoingBeam2:
          case Particle::Parton2:
            part.addMother(event_->oneWithRole(Particle::IncomingBeam2));
            break;
          case Particle::Intermediate:
            part.addMother(event_->oneWithRole(Particle::Parton1));
            part.addMother(event_->oneWithRole(Particle::Parton2));
            break;
          case Particle::CentralSystem:
            part.addMother(event_->oneWithRole(Particle::Intermediate));
            break;
          default:
            break;
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
