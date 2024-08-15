/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const PDG::Id& pdg) {
    if (HeavyIon::isHI(pdg))
      return os << HeavyIon::fromPdgId(pdg);
    return os << PDG::get().name(pdg);
  }

  PDG::PDG() {
    // PDG id, name, description, colour, mass, width, charge, is fermion
    define(ParticleProperties(invalid, "invalid", "invalid", 0, -1., -1., {}, false));
    define(ParticleProperties(diffractiveProton, "diff_proton", "p\u002A", 0, 0., 0., {-3, 3}, false));
    define(ParticleProperties(pomeron, "pomeron", "\u2119", 0, 0., 0., {}, false));
    define(ParticleProperties(reggeon, "reggeon", "\u211D", 0, 0., 0., {}, false));
  }

  PDG& PDG::get() {
    static PDG instance;
    return instance;
  }

  bool PDG::has(spdgid_t id) const { return particles_.count(std::abs(id)) > 0; }

  const ParticleProperties& PDG::operator()(spdgid_t id) const {
    if (auto it = particles_.find(std::abs(id)); it != particles_.end())
      return it->second;
    CG_DEBUG("PDG").log([this](auto& log) {
      log << "List of particles registered in the PDG runtime database:\n";
      dump(&log.stream());
    });
    throw CG_ERROR("PDG") << "No particle with PDG id " << id << " in the catalogue.";
  }

  ParticleProperties& PDG::operator[](spdgid_t id) { return particles_[std::abs(id)]; }

  void PDG::define(const ParticleProperties& props) {
    if (props.pdgid == PDG::invalid && props.name != "invalid")
      throw CG_FATAL("PDG:define") << "Trying to define a particle with invalid PDG id: " << props << ".";
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

  pdgids_t PDG::particles() const {
    pdgids_t out;
    std::transform(
        particles_.begin(), particles_.end(), std::back_inserter(out), [](const auto& pt) { return pt.first; });
    return out;
  }

  const std::string& PDG::name(spdgid_t id) const {
    if (const auto& human_name = operator()(id).human_name; !human_name.empty())
      return human_name;
    return operator()(id).name;
  }

  double PDG::colours(spdgid_t id) const { return operator()(id).colours; }

  double PDG::mass(spdgid_t id) const {
    if (HeavyIon::isHI(id))
      return HeavyIon::fromPdgId(id).mass();
    return operator()(id).mass;
  }

  double PDG::width(spdgid_t id) const { return operator()(id).width; }

  double PDG::charge(spdgid_t id) const { return operator()(id).integerCharge() * (id / std::abs(id)) * 1. / 3.; }

  std::vector<double> PDG::charges(spdgid_t id) const {
    std::vector<double> chs;
    for (const auto& ch : operator()(id).charges)
      chs.emplace_back(ch * 1. / 3);
    return chs;
  }

  size_t PDG::size() const { return particles_.size(); }

  void PDG::dump(std::ostream* os) const {
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
    std::ostringstream oss;
    oss << "List of particles registered:";
    for (const auto& prt : tmp)
      if (prt.first != PDG::invalid)
        oss << utils::format(
            "\n%16s %-32s\tcharges: {%6s}, colour factor: %1d, mass: %8.4f GeV/c^2, width: %6.3f GeV.",
            utils::colourise(std::to_string(prt.second.pdgid), utils::Colour::none, utils::Modifier::italic).data(),
            (utils::boldify(prt.second.name) + " " + (prt.second.fermion ? "fermion" : "boson") + ":").data(),
            utils::merge(prt.second.charges, ",").data(),
            prt.second.colours,
            prt.second.mass,
            prt.second.width);
    if (os)
      (*os) << oss.str();
    else
      CG_INFO("PDG") << oss.str();
  }
}  // namespace cepgen
