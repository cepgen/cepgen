/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2022  Laurent Forthomme
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

#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  ParticleProperties::ParticleProperties(pdgid_t ppdgid,
                                         const std::string& pname,
                                         const std::string& pdescr,
                                         int pcolours,
                                         double pmass,
                                         double pwidth,
                                         int pcharge,
                                         bool pfermion)
      : pdgid(ppdgid),
        name(pname),
        descr(pdescr),
        colours(pcolours),
        mass(pmass),
        width(pwidth),
        charge(pcharge),
        fermion(pfermion) {
    (*this)
        .add("pdgid", pdgid)
        .add("name", name)
        .add("description", descr)
        .add("colours", colours)
        .add("mass", mass)
        .add("width", width)
        .add("charge", charge)
        .add("fermion", fermion);
  }

  ParticleProperties::ParticleProperties(const ParametersList& params) : SteeredObject(params) {
    (*this)
        .add("pdgid", pdgid)
        .add("name", name)
        .add("description", descr)
        .add("colours", colours)
        .add("mass", mass)
        .add("width", width)
        .add("charge", charge)
        .add("fermion", fermion);
  }

  ParametersDescription ParticleProperties::description() {
    auto pdesc = ParametersDescription();
    pdesc.add<pdgid_t>("pdgid", 0).setDescription("PDG unique identifier");
    pdesc.add<std::string>("name", "n/a").setDescription("particle computer-readable name");
    pdesc.add<std::string>("description", "n/a").setDescription("particle human-readable name");
    pdesc.add<int>("colours", 0).setDescription("colour factor");
    pdesc.add<double>("mass", 0.).setDescription("particle mass (in GeV/c^2)");
    pdesc.add<double>("width", 0.).setDescription("particle width (in GeV)");
    pdesc.add<int>("charge", 0).setDescription("electric charge (in units of e)");
    pdesc.add<bool>("fermion", false).setDescription("is the particle following the Fermi-Dirac statistics?");
    return pdesc;
  }

  bool ParticleProperties::operator==(const ParticleProperties& oth) const {
    if (pdgid != oth.pdgid)
      return false;
    if (mass != oth.mass)
      return false;
    if (charge != oth.charge)
      return false;
    if (width != oth.width)
      return false;
    if (fermion != oth.fermion)
      return false;
    if (colours != oth.colours)
      return false;
    return true;
  }

  std::ostream& operator<<(std::ostream& os, const ParticleProperties& prop) {
    return os << (prop.name.empty() ? "unnamed" : prop.name) << "{"
              << "id=" << prop.pdgid << ",desc=" << prop.descr << ",colours=" << prop.colours << ",mass=" << prop.mass
              << ",width=" << prop.width << ",charge=" << prop.charge << (prop.fermion ? ",fermion" : "") << "}";
  }
}  // namespace cepgen
