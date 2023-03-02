/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#ifndef CepGen_Physics_Beam_h
#define CepGen_Physics_Beam_h

#include <iosfwd>
#include <memory>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/FormFactors/FormFactors.h"
#include "CepGen/PartonFluxes/PartonFlux.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  /// Incoming beams characteristics
  class Beam : public SteeredObject<Beam> {
  public:
    explicit Beam(const ParametersList&);  ///< Default constructor

    static ParametersDescription description();

    /// Human-readable description of a beam particle/system
    friend std::ostream& operator<<(std::ostream&, const Beam&);

    /// Type of beam treatment
    enum class Mode {
      invalid = 0,
      ProtonElastic = 1,     ///< Elastic scattering from proton
      ProtonInelastic = 2,   ///< Inelastic scattering from proton (according to the proton structure functions set)
      PointLikeScalar = 3,   ///< Trivial, spin-0 emission
      PointLikeFermion = 4,  ///< Trivial, spin-1/2 emission
      CompositeScalar = 5,   ///< Composite pion emission
      HIElastic = 10,        ///< Elastic scattering from heavy ion
      Other = 6,             ///< Other beam type
    };
    /// Human-readable format of a beam mode (elastic/dissociative parts)
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    const Mode& mode() const { return mode_; }
    /// Initialise the fluxes evaluator object
    void initialise();
    /// Is the beam particle expected to be fragmented after emission?
    bool fragmented() const;

    /// Beam particle PDG id
    pdgid_t pdgId() const { return pdg_; }
    /// Set the beam particle PDG id
    Beam& setPdgId(pdgid_t pdg) {
      pdg_ = pdg;
      return *this;
    }
    /// Scattered parton PDG id
    pdgid_t daughterId() const;
    /// Beam particle 4-momentum
    const Momentum& momentum() const { return momentum_; }
    /// Set the beam particle 4-momentum
    Beam& setMomentum(const Momentum& mom) {
      momentum_ = mom;
      return *this;
    }

    /// Scalar parton flux modelling
    const PartonFlux& flux() const;
    /// Compute the scalar parton flux given its modelling
    double flux(double x, double q2, double mx2 = -1.) const;

  private:
    pdgid_t pdg_;                       ///< PDG identifier for the beam
    Momentum momentum_;                 ///< Incoming particle momentum
    Mode mode_;                         ///< Beam treatment mode
    std::unique_ptr<PartonFlux> flux_;  ///< Incoming parton flux evaluator
  };
}  // namespace cepgen

#endif
