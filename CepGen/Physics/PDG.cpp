/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const PDG::Id& pdg) { return os << PDG::get().name(pdg); }

  PDG::PDG() {
    // PDG id, name, description, colour, mass, width, charge, is fermion
    define(ParticleProperties(invalid, "invalid", "invalid", 0, -1., -1., 0, false));
    define(ParticleProperties(diffractiveProton, "diff_proton", "p\u002A", 0, 0., 0., 3, false));
    define(ParticleProperties(pomeron, "pomeron", "\u2119", 0, 0., 0., 0, false));
    define(ParticleProperties(reggeon, "reggeon", "\u211D", 0, 0., 0., 0, false));
  }

  PDG& PDG::get() {
    static PDG instance;
    return instance;
  }

  bool PDG::has(pdgid_t id) const { return particles_.count(id) > 0; }

  const ParticleProperties& PDG::operator()(pdgid_t id) const {
    auto it = particles_.find(id);
    if (it != particles_.end())
      return it->second;
    dump();
    throw CG_FATAL("PDG") << "No particle with PDG id " << id << " in the catalogue.";
  }

  ParticleProperties& PDG::operator[](pdgid_t id) { return particles_[id]; }

  void PDG::define(const ParticleProperties& props) {
    CG_DEBUG("PDG:define").log([&](auto& log) {
      if (has(props.pdgid))
        log << "Updating the properties of a particle with PDG id=" << props.pdgid << ".\n\t"
            << "Old properties: " << operator()(props.pdgid) << ",\n\t"
            << "New properties: " << props << ".";
      else
        log << "Adding a new particle with PDG id=" << std::setw(8) << props.pdgid << ", properties: " << props << ".";
    });
    particles_[props.pdgid] = props;
  }

  const std::vector<pdgid_t> PDG::particles() const {
    std::vector<pdgid_t> out;
    std::transform(
        particles_.begin(), particles_.end(), std::back_inserter(out), [](const auto& pt) { return pt.first; });
    return out;
  }

  const std::string& PDG::name(pdgid_t id) const {
    const auto& descr = operator()(id).descr;
    if (!descr.empty())
      return descr;
    return operator()(id).name;
  }

  double PDG::colours(pdgid_t id) const { return operator()(id).colours; }

  double PDG::mass(pdgid_t id) const { return operator()(id).mass; }

  double PDG::width(pdgid_t id) const { return operator()(id).width; }

  double PDG::charge(pdgid_t id) const { return operator()(id).charge / 3.; }

  size_t PDG::size() const { return particles_.size(); }

  void PDG::dump() const {
    //--- first build a sorted vector out of the (unsorted) map
    std::vector<std::pair<pdgid_t, ParticleProperties> > tmp;
    std::transform(particles_.begin(), particles_.end(), std::back_inserter(tmp), [](const auto& prt) {
      return std::pair<pdgid_t, ParticleProperties>{prt.first, prt.second};
    });
    std::sort(tmp.begin(),
              tmp.end(),
              [](const std::pair<pdgid_t, ParticleProperties>& a, const std::pair<pdgid_t, ParticleProperties>& b) {
                return a.first < b.first;
              });
    //--- then the proper dump begins
    CG_INFO("PDG").log([&tmp](auto& info) {
      info << "List of particles registered:";
      for (const auto& prt : tmp)
        if (prt.first != PDG::invalid)
          info << utils::format(
              "\n%20s %-32s\tcharge: %2de, colour factor: %1d, mass: %8.4f GeV/c^2, width: %6.3f GeV.",
              utils::colourise(std::to_string(prt.second.pdgid), utils::Colour::none, utils::Modifier::italic).data(),
              (utils::boldify(prt.second.name) + " " + (prt.second.fermion ? "fermion" : "boson") + ":").data(),
              prt.second.charge / 3,
              prt.second.colours,
              prt.second.mass,
              prt.second.width);
    });
  }
}  // namespace cepgen
