/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include <memory>

#include "CepGen/PartonFluxes/PartonFlux.h"
#include "CepGen/Physics/Momentum.h"

namespace cepgen {
  /// Incoming beams characteristics
  class Beam final : public SteeredObject<Beam> {
  public:
    explicit Beam(const ParametersList&);  ///< Default constructor

    static ParametersDescription description();

    friend std::ostream& operator<<(std::ostream&, const Beam&);  ///< Human-readable description of beam system

    void initialise();  ///< Initialise the fluxes evaluator object

    inline bool elastic() const { return elastic_; }  ///< Does the beam remain on-shell after parton emission?
    /// Specify if the beam remains on-shell after parton emission
    inline Beam& setElastic(bool el) {
      elastic_ = el;
      return *this;
    }

    inline spdgid_t integerPdgId() const { return pdg_id_; }  ///< Beam particle PDG id
    /// Set the beam particle PDG id
    inline Beam& setIntegerPdgId(spdgid_t pdg) {
      pdg_id_ = pdg;
      return *this;
    }

    inline const Momentum& momentum() const { return momentum_; }  ///< Beam particle 4-momentum
    /// Set the beam particle 4-momentum
    inline Beam& setMomentum(const Momentum& mom) {
      momentum_ = mom;
      return *this;
    }

    inline const ParametersList& formFactors() const { return form_factors_; }        ///< Form factors parameters
    inline const ParametersList& partonFluxParameters() const { return flux_info_; }  ///< Parton flux modelling

  private:
    spdgid_t pdg_id_{0};           ///< PDG identifier for the beam
    Momentum momentum_;            ///< Incoming particle momentum
    ParametersList form_factors_;  ///< Form factors modelling parameters
    ParametersList flux_info_;     ///< Incoming parton flux parameters
    bool elastic_{true};           ///< Elastic parton emission?
  };
}  // namespace cepgen

#endif
