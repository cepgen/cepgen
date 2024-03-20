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

#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include <cstddef>  // size_t

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
      neutron = 2112,
      proton = 2212,
      diffractiveProton = 9902210
    };

    /// A class-in-the-middle PDG identifier for printout operations
    class Id {
    public:
      inline Id(pdgid_t pdgid) : pdgid_(pdgid) {}                 ///< Construct an object from a PDG identifier
      inline operator pdgid_t() const { return pdgid_; }          ///< Recasting operator to get back the pdgid_t number
      friend std::ostream& operator<<(std::ostream&, const Id&);  ///< Human-readable PDG name

    private:
      pdgid_t pdgid_;
    };

    static PDG& get();  ///< Retrieve a unique instance of this particles info collection
    PDG(const PDG&) = delete;
    void operator=(const PDG&) = delete;
    ~PDG() = default;  ///< Default destructor

    void define(const ParticleProperties&);    ///< Add a new particle definition to the library
    pdgids_t particles() const;                ///< All particles ids in this library
    void dump(std::ostream* = nullptr) const;  ///< Dump all particles in this library
    size_t size() const;                       ///< Number of particles defined in this library

    //--- per-particles information

    bool has(spdgid_t) const;                              ///< Is the particle defined for a given PDG id
    const ParticleProperties& operator()(spdgid_t) const;  ///< All physical properties for one particle
    ParticleProperties& operator[](spdgid_t id);           /// Accessor for particle properties
    const std::string& name(spdgid_t) const;               ///< Human-readable name for this particle
    double colours(spdgid_t) const;                        ///< Colour factor for this particle
    double mass(spdgid_t) const;                           ///< Particle mass (in GeV)
    double width(spdgid_t) const;                          ///< Resonance width (in GeV)
    double charge(spdgid_t) const;                         ///< Electric charge (in \f$e\f$) for this particle
    std::vector<double> charges(spdgid_t) const;  ///< Electric charges (in \f$e\f$) for this particle and anti-particles

  private:
    explicit PDG();
    std::unordered_map<pdgid_t, ParticleProperties> particles_;  ///< Collection of properties, indexed by PDG id
  };
}  // namespace cepgen

#endif
