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

#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <iosfwd>
#include <string>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  /// Alias for the integer-like particle PDG id
  typedef unsigned long long pdgid_t;
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties : SteeredObject<ParticleProperties> {
    explicit ParticleProperties(pdgid_t ppdgid = 0ull,  // PDG::invalid
                                const std::string& pname = "",
                                const std::string& pdescr = "",
                                int pcolours = -1,
                                double pmass = -1.,
                                double pwidth = -1.,
                                int pcharge = 0.,
                                bool pfermion = false)
        : ParticleProperties(ParametersList()
                                 .set("pdgid", ppdgid)
                                 .set("name", pname)
                                 .set("description", pdescr)
                                 .set("colours", pcolours)
                                 .set("mass", pmass)
                                 .set("width", pwidth)
                                 .set("charge", pcharge)
                                 .set("fermion", pfermion)) {}
    explicit ParticleProperties(const ParametersList& params) : SteeredObject(params) {
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
    static ParametersDescription description() {
      auto pdesc = ParametersDescription();
      pdesc.add<pdgid_t>("pdgid", 0).setDescription("PDG unique identifier");
      pdesc.add<std::string>("name", "").setDescription("particle computer-readable name");
      pdesc.add<std::string>("description", "").setDescription("particle human-readable name");
      pdesc.add<int>("colours", 0).setDescription("colour factor");
      pdesc.add<double>("mass", 0.).setDescription("particle mass (in GeV/c^2)");
      pdesc.add<double>("width", 0.).setDescription("particle width (in GeV)");
      pdesc.add<int>("charge", 0).setDescription("electric charge (in units of e)");
      pdesc.add<bool>("fermion", false).setDescription("is the particle following the Fermi-Dirac statistics?");
      return pdesc;
    }

    pdgid_t pdgid{0ull};  ///< PDG identifier
    std::string name;     ///< Particle name
    std::string descr;    ///< Human-readable name
    int colours{0};       ///< Colour factor
    double mass{0.};      ///< Mass, in GeV/c\f$^2\f$
    double width{0.};     ///< Decay width, in GeV/c\f$^2\f$
    int charge{0};        ///< Electric charge, in \f$e\f$/3
    bool fermion{false};  ///< Is the particle a fermion?
    friend std::ostream& operator<<(std::ostream&, const ParticleProperties&);
  };
}  // namespace cepgen

#endif
