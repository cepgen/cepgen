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
#include "CepGen/Utils/String.h"

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
      if (params.has<Limits>("phiptdiff")) {
        CG_WARNING("Central") << "\"phiptdiff\" parameter is deprecated! "
                              << "Please use \"phidiff\" instead.";
        params.fill("phiptdiff", phi_diff);  // legacy
      }
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

    //------------------------------------------------------------------

    Initial::Initial(const ParametersList& params) : SteeredObject(params) {
      (*this).add("q21", q2_1).add("q22", q2_2).add("qt", qt).add("phi", phi);
    }

    void Initial::setParameters(const ParametersList& params) {
      if (params.empty())
        return;
      SteeredObject::setParameters(params);
      if (const auto q2lim = params.get<Limits>("q2"); q2lim.valid()) {  // symmetric Q^2 cut specified
        q2_1 = q2lim;
        q2_2 = q2lim;
      }
      for (auto* q2 : {&q2_1, &q2_2})
        if (q2->max() <= 0.) {
          CG_WARNING("Initial:setParameters") << "Maximum parton virtuality (" << *q2 << ") is invalid. "
                                              << "It is now set to " << 1.e4 << " GeV^2.";
          q2->max() = 1.e4;
        }
    }

    bool Initial::contain(const Particles& parts, const Event*) const {
      for (const auto& part : parts) {
        const auto& mom = part.momentum();
        if (!qt.contains(mom.pt()))
          return false;
      }
      if (parts.size() == 2) {
        if (!q2_1.contains(parts.at(0).momentum().mass2()))
          return false;
        if (!q2_2.contains(parts.at(1).momentum().mass2()))
          return false;
        if (phi.valid() && !phi.contains(parts.at(0).momentum().deltaPhi(parts.at(1).momentum())))
          return false;
      }
      return true;
    }

    ParametersDescription Initial::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("q21", Limits{0., 1.e5}).setDescription("Positive-z virtuality (GeV^2)");
      desc.add<Limits>("q22", Limits{0., 1.e5}).setDescription("Negative-z virtuality (GeV^2)");
      desc.add<Limits>("qt", Limits{}).setDescription("Transverse virtuality (GeV)");
      desc.add<Limits>("phi", Limits{}).setDescription("Partons D(phi) (rad)");
      return desc;
    }

    //------------------------------------------------------------------

    Remnants::Remnants(const ParametersList& params) : SteeredObject(params) {
      (*this).add("mx", mx).add("yj", yj).add("xi", xi);
    }

    void Remnants::setParameters(const ParametersList& params) {
      if (params.empty())
        return;
      SteeredObject::setParameters(params);
      if (mx.min() < MX_MIN) {
        CG_WARNING("CutsList:setParameters") << "Minimum diffractive mass range (" << mx << ") is invalid. "
                                             << "It is now set to " << MX_MIN << " GeV/c^2.";
        mx.min() = MX_MIN;
      }
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

    ParametersDescription Remnants::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("mx", Limits{Remnants::MX_MIN, 1.e3}).setDescription("Diffractive mass (GeV/c^2)");
      desc.add<Limits>("yj", Limits{}).setDescription("Diffractive jet rapidity");
      desc.add<Limits>("xi", Limits{}).setDescription("Longit.fract.mom. loss (\"xi\")");
      return desc;
    }
  }  // namespace cuts

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  CutsList::CutsList(const ParametersList& params)
      : SteeredObject(params),
        initial(steer<ParametersList>("initial")),
        central(steer<ParametersList>("central")),
        remnants(steer<ParametersList>("remnants")) {}

  void CutsList::setParameters(const ParametersList& params) {
    if (params.empty())
      return;
    initial.setDescribedParameters(params);
    central.setDescribedParameters(params);
    remnants.setDescribedParameters(params);
    if (params.has<ParametersList>("cuts")) {  // per-particle cuts
      const auto per_part_cuts = params.get<ParametersList>("cuts");
      for (const auto& part : per_part_cuts.keys())
        central_particles[(pdgid_t)std::stoi(part)].setDescribedParameters(per_part_cuts.get<ParametersList>(part));
    }

    // override the parameters from sub-parameters content
    params_.set("initial", initial.parameters())
        .set("central", central.parameters())
        .set("remnants", remnants.parameters());
    for (const auto& cuts_vs_part : central_particles)  // per-PDGid selection
      params_.operator[]<ParametersList>("cuts").set<ParametersList>(std::to_string(cuts_vs_part.first),
                                                                     cuts_vs_part.second.parameters());
    CG_DEBUG("CutsList:setParameters") << "User specified the following cuts list:\n" << *this << ".";
    for (const auto& key : params_.keys()) {
      if (key == "initial" || key == "central" || key == "remnants" || key == "cuts") {
        auto& cuts = params_.operator[]<ParametersList>(key);
        for (const auto& lim_key : cuts.keysOf<Limits>())
          if (auto& lim = cuts.operator[]<Limits>(lim_key); lim.min() == 0. && lim.max() == 0.) {
            CG_WARNING("CutsList:setParameters")
                << "Unset the range for '" << key << "/" << lim_key << "' from " << lim << ".";
            lim = Limits{};
          }
      } else
        params_.erase(key);
    }
    initial.setParameters(steer<ParametersList>("initial"));
    central.setParameters(steer<ParametersList>("central"));
    remnants.setParameters(steer<ParametersList>("remnants"));
  }

  ParametersDescription CutsList::description() {
    auto desc = ParametersDescription();
    desc.add("initial", cuts::Initial::description());
    desc.add("central", cuts::Central::description());
    desc.add("remnants", cuts::Remnants::description());
    return desc;
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
