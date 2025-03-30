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

#ifndef CepGen_EventFilter_EventBrowser_h
#define CepGen_EventFilter_EventBrowser_h

#include <regex>

#include "CepGen/Event/Particle.h"

namespace cepgen {
  class Event;
}
namespace cepgen::utils {
  /// User-friendly browser for the Event content
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class EventBrowser {
  public:
    EventBrowser() = default;
    double get(const Event& ev, const std::string& var) const;  ///< Get/compute a variable value

  private:
    /// Retrieve a particle named variable
    double variable(const Event&, const Particle&, const std::string&) const;
    /// Retrieve a two-particle system named variable
    double variable(const Event&, const Particle&, const Particle&, const std::string&) const;
    /// Retrieve a whole event named variable
    static double variable(const Event&, const std::string&);

    static const std::regex rgx_select_id_, rgx_select_id2_, rgx_select_role_, rgx_select_role2_;
    static constexpr double INVALID_OUTPUT = -999.;

    //--- auxiliary helper maps
    const std::unordered_map<std::string, Particle::Role> role_str_ = {{"ib1", Particle::Role::IncomingBeam1},
                                                                       {"ib2", Particle::Role::IncomingBeam2},
                                                                       {"ob1", Particle::Role::OutgoingBeam1},
                                                                       {"ob2", Particle::Role::OutgoingBeam2},
                                                                       {"pa1", Particle::Role::Parton1},
                                                                       {"pa2", Particle::Role::Parton2},
                                                                       {"cs", Particle::Role::CentralSystem},
                                                                       {"int", Particle::Role::Intermediate}};
    typedef double (Momentum::*pMethod)() const;
    /// Mapping of string variables to momentum getter methods
    const std::unordered_map<std::string, pMethod> m_mom_str_ = {
        {"px", &Momentum::px},        {"py", &Momentum::py},      {"pz", &Momentum::pz},
        {"pt", &Momentum::pt},        {"pt2", &Momentum::pt2},    {"eta", &Momentum::eta},
        {"phi", &Momentum::phi},      {"m", &Momentum::mass},     {"m2", &Momentum::mass2},
        {"mt", &Momentum::massT},     {"mt2", &Momentum::massT2}, {"e", &Momentum::energy},
        {"e2", &Momentum::energy2},   {"et", &Momentum::energyT}, {"et2", &Momentum::energyT2},
        {"p", &Momentum::p},          {"p2", &Momentum::p2},      {"th", &Momentum::theta},
        {"y", &Momentum::rapidity},   {"beta", &Momentum::beta},  {"gamma", &Momentum::gamma},
        {"gamma2", &Momentum::gamma2}};
    typedef double (Momentum::*pMethodOth)(const Momentum&) const;
    const std::unordered_map<std::string, pMethodOth> m_two_mom_str_ = {{"deta", &Momentum::deltaEta},
                                                                        {"dphi", &Momentum::deltaPhi},
                                                                        {"dpt", &Momentum::deltaPt},
                                                                        {"dr", &Momentum::deltaR}};
  };
}  // namespace cepgen::utils

#endif
