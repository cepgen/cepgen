/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#include <iostream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  ParticleProperties::ParticleProperties(const ParametersList& params) : SteeredObject(params) {
    (*this)
        .add("pdgid", pdgid)
        .add("name", name)
        .add("description", descr)
        .add("colours", colours)
        .add("mass", mass)
        .add("width", width)
        .add("charges", charges)
        .add("fermion", fermion);
  }

  ParticleProperties::ParticleProperties(pdgid_t ppdgid,
                                         const std::string& pname,
                                         const std::string& pdescr,
                                         int pcolours,
                                         double pmass,
                                         double pwidth,
                                         const std::vector<int>& pcharges,
                                         bool pfermion)
      : ParticleProperties(ParametersList()
                               .set("pdgid", ppdgid)
                               .set("name", pname)
                               .set("description", pdescr)
                               .set("colours", pcolours)
                               .set("mass", pmass)
                               .set("width", pwidth)
                               .set("charges", pcharges)
                               .set("fermion", pfermion)) {}

  short ParticleProperties::integerCharge() const {
    if (charges.empty())
      return 0;
    if (charges.size() > 2)
      throw CG_ERROR("ParticleProperties:integerCharge")
          << "Multiple charges are possible for the given particle: " << charges << ".";
    return charges.at(0);
  }

  ParametersDescription ParticleProperties::description() {
    auto desc = ParametersDescription();
    desc.add<pdgid_t>("pdgid", 0).setDescription("PDG unique identifier");
    desc.add<std::string>("name", "n/a").setDescription("particle computer-readable name");
    desc.add<std::string>("description", "n/a").setDescription("particle human-readable name");
    desc.add<int>("colours", 0).setDescription("colour factor");
    desc.add<double>("mass", 0.).setDescription("particle mass (in GeV/c^2)");
    desc.add<double>("width", 0.).setDescription("particle width (in GeV)");
    desc.add<std::vector<int> >("charges", {}).setDescription("possible electric charges (in units of e)");
    desc.add<bool>("fermion", false).setDescription("is the particle following the Fermi-Dirac statistics?");
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const ParticleProperties& prop) {
    return os << (prop.name.empty() ? "unnamed" : prop.name) << "{"
              << "pdgid=" << prop.pdgid << ",desc=" << prop.descr << ",colours=" << prop.colours
              << ",mass=" << prop.mass << ",width=" << prop.width << ",charges={" << utils::merge(prop.charges, ", ")
              << "}" << (prop.fermion ? ",fermion" : "") << "}";
  }
}  // namespace cepgen
