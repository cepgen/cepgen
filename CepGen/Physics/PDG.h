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

#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include <cstddef>  // size_t
#include <unordered_map>
#include <vector>

#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  /// A singleton holding all physics constants associated to particles
  class PDG {
  public:
    /// PDG ids of all known particles
    /// \note From \cite Beringer:1900zz :
    /// > The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.
    enum PdgIdEnum : pdgid_t {
      invalid = 0,
      down = 1,
      up = 2,
      electron = 11,
      muon = 13,
      tau = 15,
      gluon = 21,
      photon = 22,
      W = 24,
      pomeron = 990,
      reggeon = 110,
      piZero = 111,
      piPlus = 211,
      eta = 221,
      phi1680 = 100333,
      proton = 2212,
      diffractiveProton = 9902210
    };

    /// A class-in-the-middle PDG identifier for printout operations
    class Id {
    public:
      /// Construct an object from a PDG identifier
      Id(pdgid_t pdgid) : pdgid_(pdgid) {}
      /// Recasting operator to get back the pdgid_t number
      operator pdgid_t() const { return pdgid_; }
      /// Human-readable PDG name
      friend std::ostream& operator<<(std::ostream&, const Id&);

    private:
      pdgid_t pdgid_;
    };

    /// Retrieve a unique instance of this particles info collection
    static PDG& get();
    PDG(const PDG&) = delete;
    void operator=(const PDG&) = delete;
    /// Default destructor
    ~PDG() = default;

    /// Add a new particle definition to the library
    void define(const ParticleProperties& props);
    const std::vector<pdgid_t> particles() const;  ///< All particles ids in this library
    void dump() const;                             ///< Dump all particles in this library
    size_t size() const;                           ///< Number of particles defined in this library

    //--- per-particles information

    bool has(pdgid_t) const;                              ///< Is the particle defined for a given PDG id
    const ParticleProperties& operator()(pdgid_t) const;  ///< All physical properties for one particle
    const std::string& name(pdgid_t) const;               ///< Human-readable name for this particle
    double colours(pdgid_t) const;                        ///< Colour factor for this particle
    double mass(pdgid_t) const;                           ///< Particle mass (in GeV)
    double width(pdgid_t) const;                          ///< Resonance width (in GeV)
    double charge(pdgid_t) const;                         ///< Electric charge (in \f$e\f$) for this particle

  private:
    explicit PDG();
    /** \note Indexing variable: PDG id of particle */
    std::unordered_map<pdgid_t, ParticleProperties> particles_;
  };
}  // namespace cepgen

#endif
