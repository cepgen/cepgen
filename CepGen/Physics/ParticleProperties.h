/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
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

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  using pdgid_t = unsigned long long;       ///< Alias for the integer-like particle PDG id
  using pdgids_t = std::vector<pdgid_t>;    ///< Alias for a collection of particles PDG ids
  using spdgid_t = long long;               ///< Alias for a signed particle PDG id (adding charge information)
  using spdgids_t = std::vector<spdgid_t>;  ///< Alias for a collection of particles signed PDG ids
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties final : SteeredObject<ParticleProperties> {
    explicit ParticleProperties(const ParametersList&);
    explicit ParticleProperties(pdgid_t _pdg_id = 0ull,  // PDG::invalid
                                const std::string& _name = "",
                                const std::string& _description = "",
                                int colours = -1,
                                double _mass = -1.,
                                double _width = -1.,
                                const std::vector<int>& _charges = {},
                                bool _fermion = false);

    static ParametersDescription description();

    friend std::ostream& operator<<(std::ostream&, const ParticleProperties&);

    short integerCharge() const;  ///< Integer charge, in \f$e\f$/3

    pdgid_t pdgid{0ull};         ///< PDG identifier
    std::string name{};          ///< Particle name
    std::string human_name{};    ///< Human-readable name
    int colours{0};              ///< Colour factor
    double mass{0.};             ///< Mass, in GeV/c\f$^2\f$
    double width{0.};            ///< Decay width, in GeV/c\f$^2\f$
    std::vector<int> charges{};  ///< Electric charges, in \f$e\f$/3
    bool fermion{false};         ///< Is the particle a fermion?
  };
}  // namespace cepgen

#endif
