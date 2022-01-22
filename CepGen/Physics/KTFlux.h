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

#ifndef CepGen_Physics_KTFlux_h
#define CepGen_Physics_KTFlux_h

#include <iosfwd>

#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace formfac {
    class Parameterisation;
  }
  namespace strfun {
    class Parameterisation;
  }
  struct HeavyIon;
  /// Collection of fundamental constants for \f$k_{\rm T}\f$ fluxes definition
  struct KTFluxParameters {
    static const double kMinKTFlux;  ///< Minimal value taken for a \f$\k_{\rm T}\f$-factorised flux
  };
  /// Type of incoming partons fluxes
  enum class KTFlux {
    invalid = -1,                    ///< Invalid flux
    P_Photon_Elastic = 0,            ///< Elastic photon emission from proton
    P_Photon_Elastic_Budnev = 10,    ///< Elastic photon emission from proton (Budnev flux approximation)
    P_Photon_Inelastic = 1,          ///< Inelastic photon emission from proton
    P_Photon_Inelastic_Budnev = 11,  ///< Inelastic photon emission from proton (Budnev flux approximation)
    P_Gluon_KMR = 20,                ///< Inelastic gluon emission from proton (KMR flux modelling)
    HI_Photon_Elastic = 100          ///< Elastic photon emission from heavy ion (from Starlight \cite Klein:2016yzr)
  };
  /// Human version of the flux name
  std::ostream& operator<<(std::ostream&, const KTFlux&);
  /// \brief Compute the flux for a given parton \f$(x,k_{\rm T})\f$
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
  /// \param[in] ff Form factors evaluator
  /// \param[in] sf Structure functions evaluator
  /// \param[in] mi2 Incoming particle squared mass
  /// \param[in] mf2 Outgoing diffractive squared mass
  double ktFlux(const KTFlux& type,
                double x,
                double kt2,
                formfac::Parameterisation& ff,
                strfun::Parameterisation& sf,
                double mi2,
                double mf2);
  /// \brief Compute the flux (from heavy ion) for a given parton \f$(x,k_{\rm T})\f$
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
  /// \param[in] hi Heavy ion properties
  double ktFlux(const KTFlux& type, double x, double kt2, const HeavyIon& hi);
}  // namespace cepgen

#endif
