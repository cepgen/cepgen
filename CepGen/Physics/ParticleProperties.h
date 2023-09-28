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

#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <iosfwd>
#include <string>
#include <vector>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  /// Alias for the integer-like particle PDG id
  typedef unsigned long long pdgid_t;
  typedef std::vector<pdgid_t> pdgids_t;
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties final : SteeredObject<ParticleProperties> {
    explicit ParticleProperties(const ParametersList&);
    explicit ParticleProperties(pdgid_t ppdgid = 0ull,  // PDG::invalid
                                const std::string& pname = "",
                                const std::string& pdescr = "",
                                int pcolours = -1,
                                double pmass = -1.,
                                double pwidth = -1.,
                                int pcharge = 0.,
                                bool pfermion = false);

    static ParametersDescription description();

    friend std::ostream& operator<<(std::ostream&, const ParticleProperties&);

    pdgid_t pdgid{0ull};  ///< PDG identifier
    std::string name{};   ///< Particle name
    std::string descr{};  ///< Human-readable name
    int colours{0};       ///< Colour factor
    double mass{0.};      ///< Mass, in GeV/c\f$^2\f$
    double width{0.};     ///< Decay width, in GeV/c\f$^2\f$
    int charge{0};        ///< Electric charge, in \f$e\f$/3
    bool fermion{false};  ///< Is the particle a fermion?
  };
}  // namespace cepgen

#endif
