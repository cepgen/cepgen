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
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PartonFlux.h"

namespace cepgen {
  /// Incoming beams characteristics
  class Beam : public SteeredObject<Beam> {
  public:
    explicit Beam(const ParametersList&);  ///< Default constructor

    static ParametersDescription description();

    /// Human-readable description of a beam particle/system
    friend std::ostream& operator<<(std::ostream&, const Beam&);

    /// Initialise the fluxes evaluator object
    void initialise();

    /// Does the beam remain on-shell after parton emission?
    bool elastic() const { return elastic_; }
    /// Specify if the beam remains on-shell after parton emission
    Beam& setElastic(bool elastic) {
      elastic_ = elastic;
      return *this;
    }
    /// Beam particle PDG id
    pdgid_t pdgId() const { return pdg_id_; }
    /// Set the beam particle PDG id
    Beam& setPdgId(pdgid_t pdg) {
      pdg_id_ = pdg;
      return *this;
    }
    /// Beam particle 4-momentum
    const Momentum& momentum() const { return momentum_; }
    /// Set the beam particle 4-momentum
    Beam& setMomentum(const Momentum& mom) {
      momentum_ = mom;
      return *this;
    }
    const ParametersList& partonFluxParameters() const { return flux_info_; }

  private:
    pdgid_t pdg_id_{0};         ///< PDG identifier for the beam
    Momentum momentum_;         ///< Incoming particle momentum
    ParametersList flux_info_;  ///< Incoming parton flux parameters
    bool elastic_;              ///< Elastic parton emission?
  };
}  // namespace cepgen

#endif
