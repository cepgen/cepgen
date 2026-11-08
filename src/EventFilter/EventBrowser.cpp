/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Utils/String.h"

using namespace cepgen::utils;
using namespace std::string_literals;

const std::regex EventBrowser::rgx_select_id_("([a-zA-Z0-9]+)\\(([0-9]+)\\)", std::regex_constants::extended);
const std::regex EventBrowser::rgx_select_id2_("([a-zA-Z0-9]+)\\(([0-9]+),([0-9]+)\\)", std::regex_constants::extended);
const std::regex EventBrowser::rgx_select_role_("([a-zA-Z0-9]+)\\(([a-z]+[0-9]?)\\)", std::regex_constants::extended);
const std::regex EventBrowser::rgx_select_role2_("([a-zA-Z0-9]+)\\(([a-z]+[0-9]?),([a-z]+[0-9]?)\\)",
                                                 std::regex_constants::extended);

double EventBrowser::get(const Event& event, const std::string& variable_name) const {
  std::smatch sm;
  if (std::regex_match(variable_name, sm, rgx_select_id_))  // particle-level variables (indexed by integer id)
    return variable(event, event(std::stoul(sm[2].str())), sm[1].str());
  if (std::regex_match(variable_name, sm, rgx_select_id2_))  // particle-level variables (indexed by integer id)
    return variable(event, event(std::stoul(sm[2].str())), event(std::stoul(sm[3].str())), sm[1].str());
  const auto check_role = [&](const std::string& role, const std::string& variable) -> bool {
    if (role_str_.count(role) > 0)
      return true;
    CG_WARNING("EventBrowser") << "Invalid particle role retrieved from configuration: \"" << role << "\".\n\t"
                               << "Skipping the variable \"" << variable << "\" in the output module.";
    return false;
  };
  if (std::regex_match(variable_name, sm, rgx_select_role_)) {  // particle-level variables (indexed by role)
    if (const auto& str_role = sm[2].str(); check_role(str_role, variable_name))
      return variable(event, event(role_str_.at(str_role))[0], sm[1].str());
    return INVALID_OUTPUT;
  }
  if (std::regex_match(variable_name, sm, rgx_select_role2_)) {  // particle-level variables (indexed by role)
    if (const auto &str_role1 = sm[2].str(), &str_role2 = sm[3].str();
        check_role(str_role1, variable_name) && check_role(str_role2, variable_name))
      return variable(event, event(role_str_.at(str_role1))[0], event(role_str_.at(str_role2))[0], sm[1].str());
    return INVALID_OUTPUT;
  }
  return variable(event, variable_name);  // event-level variables
}

double EventBrowser::variable(const Event& event, const Particle& particle, const std::string& variable_name) const {
  if (m_mom_str_.count(variable_name)) {
    const auto& meth = m_mom_str_.at(variable_name);
    return (particle.momentum().*meth)();
  }
  if (variable_name == "xi") {
    if (const auto& moth = particle.mothers(); !moth.empty())
      return 1. - particle.momentum().energy() / event(*moth.begin()).momentum().energy();
    CG_WARNING("EventBrowser") << "Failed to retrieve parent particle to compute xi "
                               << "for the following particle:\n"
                               << particle;
    return INVALID_OUTPUT;
  }
  if (variable_name == "pdg")
    return static_cast<double>(particle.integerPdgId());
  if (variable_name == "charge")
    return particle.charge();
  if (variable_name == "status")
    return static_cast<double>(particle.status());
  throw CG_ERROR("EventBrowser") << "Failed to retrieve variable \"" << variable_name << "\".";
}

double EventBrowser::variable(const Event&,
                              const Particle& particle1,
                              const Particle& particle2,
                              const std::string& variable_name) const {
  if (m_two_mom_str_.count(variable_name)) {
    const auto& meth = m_two_mom_str_.at(variable_name);
    return (particle1.momentum().*meth)(particle2.momentum());
  }
  if (m_mom_str_.count(variable_name)) {
    const auto& meth = m_mom_str_.at(variable_name);
    return ((particle1.momentum() + particle2.momentum()).*meth)();
  }
  if (variable_name == "acop"s)  // two-particle acoplanarity
    return 1. - fabs(particle1.momentum().deltaPhi(particle2.momentum()) * M_1_PI);
  throw CG_ERROR("EventBrowser") << "Failed to retrieve variable \"" << variable_name << "\".";
}

double EventBrowser::variable(const Event& event, const std::string& variable_name) {
  if (variable_name == "np")  // number of particles in event (whatever the status)
    return static_cast<double>(event.size());
  if (variable_name == "nob1" ||
      variable_name == "nob2") {  // number of (stable) particles in outgoing diffractive system
    const auto& beam_particles =
        event(variable_name == "nob1" ? Particle::Role::OutgoingBeam1 : Particle::Role::OutgoingBeam2);
    return static_cast<double>(std::count_if(beam_particles.begin(), beam_particles.end(), [](const auto& particle) {
      return static_cast<int>(particle.status()) > 0;
    }));
  }
  if (variable_name == "met")  // missing transverse energy
    return event.missingMomentum().pt();
  if (variable_name == "mephi"s)  // azimuthal component of the missing transverse energy
    return event.missingMomentum().phi();
  if (variable_name == "cmEnergy")  // two-beam centre-of-mass energy
    return event.cmEnergy();
  if (startsWith(variable_name, "meta:"))  // metadata field
    return event.metadata(variable_name.substr(5));
  throw CG_ERROR("EventBrowser") << "Failed to retrieve the event-level variable \"" << variable_name << "\".";
}
