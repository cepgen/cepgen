/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2023  Laurent Forthomme
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

#include <algorithm>

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace cuts {
    //--------------------------------------------------------------------
    // physics system kinematic properties
    //--------------------------------------------------------------------

    Central::Central() : SteeredObject(ParametersList{}) {}

    Central::Central(const ParametersList& params) : SteeredObject(params) {
      (*this)
          .add("pt", pt_single)
          .add("eta", eta_single)
          .add("rapidity", rapidity_single)
          .add("energy", energy_single)
          .add("mass", mass_single)
          .add("ptsum", pt_sum)
          .add("etasum", eta_sum)
          .add("energysum", energy_sum)
          .add("invmass", mass_sum)
          .add("ptdiff", pt_diff)
          .add("dphi", phi_diff)
          .add("rapiditydiff", rapidity_diff);
    }

    ParametersDescription Central::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("pt", Limits{0.}).setDescription("Single particle pt (GeV/c)");
      desc.add<Limits>("eta", Limits{}).setDescription("Single particle eta");
      desc.add<Limits>("rapidity", Limits{}).setDescription("Single particle rapidity");
      desc.add<Limits>("energy", Limits{}).setDescription("Single particle energy (GeV)");
      desc.add<Limits>("mass", Limits{}).setDescription("Single particle mass (GeV/c^2)");
      desc.add<Limits>("ptsum", Limits{}).setDescription("System pt (GeV/c)");
      desc.add<Limits>("etasum", Limits{}).setDescription("System eta");
      desc.add<Limits>("energysum", Limits{}).setDescription("System energy (GeV)");
      desc.add<Limits>("invmass", Limits{}).setDescription("System mass (GeV/c^2)");
      desc.add<Limits>("ptdiff", Limits{}).setDescription("System D(pt) (GeV/c)");
      desc.add<Limits>("dphi", Limits{}).setDescription("System D(phi) (rad)");
      desc.add<Limits>("rapiditydiff", Limits{}).setDescription("System D(Y)");
      return desc;
    }

    bool Central::contain(const Particles& parts, const Event*) const {
      Momentum mom_sum;
      for (const auto& part : parts) {
        const auto& mom = part.momentum();
        if (!pt_single.contains(mom.pt()) || !eta_single.contains(mom.eta()) ||
            !rapidity_single.contains(mom.rapidity()) || !energy_single.contains(mom.energy()) ||
            !mass_single.contains(mom.mass()))
          return false;
        mom_sum += mom;
      }
      if (!pt_sum.contains(mom_sum.pt()) || !eta_sum.contains(mom_sum.eta()) ||
          !energy_sum.contains(mom_sum.energy()) || !mass_sum.contains(mom_sum.mass()))
        return false;
      if (parts.size() > 1) {  // look at correlations
        const auto &mom1 = parts.at(0).momentum(), &mom2 = parts.at(1).momentum();
        if (!pt_diff.contains(fabs(mom1.pt() - mom2.pt())) || !phi_diff.contains(mom1.deltaPhi(mom2)) ||
            !rapidity_diff.contains(fabs(mom1.rapidity() - mom2.rapidity())))
          return false;
      }
      return true;
    }

    Initial::Initial(const ParametersList& params) : SteeredObject(params) {
      (*this).add("q2", q2).add("qt", qt).add("phi", phi);
    }

    ParametersDescription Initial::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("q2", Limits{0., 1.e5}).setDescription("Virtuality (GeV^2)");
      desc.add<Limits>("qt", Limits{}).setDescription("Transverse virtuality (GeV)");
      desc.add<Limits>("phi", Limits{}).setDescription("Partons D(phi) (rad)");
      return desc;
    }

    bool Initial::contain(const Particles& parts, const Event*) const {
      for (const auto& part : parts) {
        const auto& mom = part.momentum();
        if (!qt.contains(mom.pt()))
          return false;
        if (!q2.contains(mom.mass2()))
          return false;
      }
      if (parts.size() == 2 && phi.valid() && !phi.contains(parts.at(0).momentum().deltaPhi(parts.at(1).momentum())))
        return false;
      return true;
    }

    Remnants::Remnants(const ParametersList& params) : SteeredObject(params) {
      (*this).add("mx", mx).add("yj", yj).add("xi", xi);
    }

    ParametersDescription Remnants::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("mx", Limits{Remnants::MX_MIN, 1.e3}).setDescription("Diffractive mass (GeV/c^2)");
      desc.add<Limits>("yj", Limits{}).setDescription("Diffractive jet rapidity");
      desc.add<Limits>("xi", Limits{}).setDescription("Longit.fract.mom. loss (\"xi\")");
      return desc;
    }

    bool Remnants::contain(const Particles& parts, const Event* evt) const {
      for (const auto& part : parts) {
        if (part.status() != Particle::Status::FinalState)
          continue;
        if (evt && xi.valid() &&
            !xi.contains(1. - part.momentum().pz() / (*evt)(*part.mothers().begin()).momentum().pz()))
          return false;
        if (!yj.contains(fabs(part.momentum().rapidity())))
          return false;
      }
      return true;
    }
  }  // namespace cuts

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  CutsList::CutsList(const ParametersList& params)
      : SteeredObject(params), initial(params_), central(params_), remnants(params_) {}

  const ParametersList& CutsList::parameters() const {
    params_ += initial.parameters() + central.parameters() + remnants.parameters();
    for (const auto& cuts_vs_part : central_particles)  // per-PDGid selection
      params_.operator[]<ParametersList>("cuts").set<ParametersList>(std::to_string(cuts_vs_part.first),
                                                                     cuts_vs_part.second.parameters());
    return SteeredObject::parameters();
  }

  void CutsList::setParameters(const ParametersList& params) {
    SteeredObject::setParameters(params);
    initial.setParameters(params_);
    central.setParameters(params_);
    remnants.setParameters(params_);
    if (params_.has<ParametersList>("cuts")) {  // per-particle cuts
      const auto& per_parts = steer<ParametersList>("cuts");
      for (const auto& part : per_parts.keys())
        central_particles[(pdgid_t)stoi(part)].setParameters(per_parts.get<ParametersList>(part));
    }

    // sanity checks
    if (params_.has<Limits>("phiptdiff")) {
      CG_WARNING("CutsList:setParameters") << "\"phiptdiff\" parameter is deprecated! "
                                           << "Please use \"phidiff\" instead.";
      params_.fill<Limits>("phiptdiff", central.phi_diff);  //legacy
    }
    if (initial.q2.max() <= 0.) {
      CG_WARNING("CutsList:setParameters") << "Maximum parton virtuality (" << initial.q2 << ") is invalid. "
                                           << "It is now set to " << 1.e4 << " GeV^2.";
      initial.q2.max() = 1.e4;
    }
    if (remnants.mx.min() < cuts::Remnants::MX_MIN) {
      CG_WARNING("CutsList:setParameters") << "Minimum diffractive mass range (" << remnants.mx << ") is invalid. "
                                           << "It is now set to " << cuts::Remnants::MX_MIN << " GeV/c^2.";
      remnants.mx.min() = cuts::Remnants::MX_MIN;
    }

    CG_DEBUG("CutsList:setParameters") << "User specified the following cuts list:\n" << *this << ".";
  }

  std::ostream& operator<<(std::ostream& os, const CutsList& cl) {
    auto dump_cuts = [&os](const auto& obj) {
      std::string sep;
      for (const auto& lim : obj.parameters().template keysOf<Limits>()) {
        const auto& limit = obj.parameters().template get<Limits>(lim);
        if (limit.valid() && obj.description().has(lim))
          os << sep << obj.description().get(lim).description() << ": " << limit, sep = ";";
      }
    };
    os << "init.system{";
    dump_cuts(cl.initial);
    os << "}, cent.system{";
    dump_cuts(cl.central);
    os << "}, remnants{";
    dump_cuts(cl.remnants);
    return os << "}";
  }
}  // namespace cepgen
