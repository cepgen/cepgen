/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    const std::regex EventBrowser::rgx_select_id_("([a-zA-Z]+)\\(([0-9]+)\\)", std::regex_constants::extended);
    const std::regex EventBrowser::rgx_select_id2_("([a-zA-Z]+)\\(([0-9]+),([0-9]+)\\)",
                                                   std::regex_constants::extended);
    const std::regex EventBrowser::rgx_select_role_("([a-zA-Z]+)\\(([a-z]+[0-9]?)\\)", std::regex_constants::extended);
    const std::regex EventBrowser::rgx_select_role2_("([a-zA-Z]+)\\(([a-z]+[0-9]?),([a-z]+[0-9]?)\\)",
                                                     std::regex_constants::extended);

    double EventBrowser::get(const Event& ev, const std::string& var) const {
      std::smatch sm;
      //--- particle-level variables (indexed by integer id)
      if (std::regex_match(var, sm, rgx_select_id_)) {
        const auto& var_name = sm[1].str();
        const auto& part = ev(std::stoul(sm[2].str()));
        return variable(ev, part, var_name);
      }
      if (std::regex_match(var, sm, rgx_select_id2_)) {
        const auto& var_name = sm[1].str();
        const auto& part1 = ev(std::stoul(sm[2].str()));
        const auto& part2 = ev(std::stoul(sm[3].str()));
        return variable(ev, part1, part2, var_name);
      }
      //--- particle-level variables (indexed by role)
      const auto check_role = [&](const std::string& role, const std::string& var) -> bool {
        bool ret = role_str_.count(role) > 0;
        if (!ret)
          CG_WARNING("EventBrowser") << "Invalid particle role retrieved from configuration: \"" << role << "\".\n\t"
                                     << "Skipping the variable \"" << var << "\" in the output module.";
        return ret;
      };
      if (std::regex_match(var, sm, rgx_select_role_)) {
        const auto& var_name = sm[1].str();
        const auto& str_role = sm[2].str();
        if (!check_role(str_role, var))
          return INVALID_OUTPUT;
        const auto& part = ev(role_str_.at(str_role))[0];
        return variable(ev, part, var_name);
      }
      if (std::regex_match(var, sm, rgx_select_role2_)) {
        const auto& var_name = sm[1].str();
        const auto& str_role1 = sm[2].str();
        const auto& str_role2 = sm[3].str();
        if (!check_role(str_role1, var) || !check_role(str_role2, var))
          return INVALID_OUTPUT;
        const auto& part1 = ev(role_str_.at(str_role1))[0];
        const auto& part2 = ev(role_str_.at(str_role2))[0];
        return variable(ev, part1, part2, var_name);
      }
      //--- event-level variables
      return variable(ev, var);
    }

    double EventBrowser::variable(const Event& ev, const Particle& part, const std::string& var) const {
      if (m_mom_str_.count(var)) {
        const auto& meth = m_mom_str_.at(var);
        return (part.momentum().*meth)();
      }
      if (var == "xi") {
        const auto& moth = part.mothers();
        if (moth.empty()) {
          CG_WARNING("EventBrowser") << "Failed to retrieve parent particle to compute xi "
                                     << "for the following particle:\n"
                                     << part;
          return INVALID_OUTPUT;
        }
        return 1. - part.momentum().energy() / ev(int(*moth.begin())).momentum().energy();
      }
      if (var == "pdg")
        return (double)part.integerPdgId();
      if (var == "charge")
        return part.charge();
      if (var == "status")
        return (double)part.status();
      throw CG_ERROR("EventBrowser") << "Failed to retrieve variable \"" << var << "\".";
    }

    double EventBrowser::variable(const Event&,
                                  const Particle& part1,
                                  const Particle& part2,
                                  const std::string& var) const {
      if (m_two_mom_str_.count(var)) {
        const auto& meth = m_two_mom_str_.at(var);
        return (part1.momentum().*meth)(part2.momentum());
      }
      if (m_mom_str_.count(var)) {
        const auto& meth = m_mom_str_.at(var);
        return ((part1.momentum() + part2.momentum()).*meth)();
      }
      if (var == "acop")
        return 1. - fabs(part1.momentum().deltaPhi(part2.momentum()) * M_1_PI);
      throw CG_ERROR("EventBrowser") << "Failed to retrieve variable \"" << var << "\".";
    }

    double EventBrowser::variable(const Event& ev, const std::string& var) {
      if (var == "np")
        return (double)ev.size();
      //if ( var == "nev" )
      //  return (double)num_evts_+1;
      if (var == "nob1" || var == "nob2") {
        const auto& bparts = ev(var == "nob1" ? Particle::Role::OutgoingBeam1 : Particle::Role::OutgoingBeam2);
        return (double)std::count_if(
            bparts.begin(), bparts.end(), [](const auto& part) { return (int)part.status() > 0; });
      }
      if (var == "met")
        return ev.missingMomentum().pt();
      if (var == "mephi")
        return ev.missingMomentum().phi();
      if (utils::startsWith(var, "meta:"))
        return ev.metadata(var.substr(5));
      throw CG_ERROR("EventBrowser") << "Failed to retrieve the event-level variable \"" << var << "\".";
    }
  }  // namespace utils
}  // namespace cepgen
