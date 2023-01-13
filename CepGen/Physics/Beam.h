/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/FormFactors/FormFactors.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  enum class KTFlux;
  namespace strfun {
    class Parameterisation;
  }
  namespace formfac {
    class Parameterisation;
  }
  struct HeavyIon;

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
      Other = 6,             ///< Other beam type
    };
    /// Human-readable format of a beam mode (elastic/dissociative parts)
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    const Mode& mode() const { return mode_; }
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
    friend std::ostream& operator<<(std::ostream&, const KTFlux&);
    const KTFlux& ktFlux() const { return kt_flux_; }

    /// Compute the electromagnetic form factors for the current beam mode
    formfac::FormFactors flux(double q2,
                              double mx2,
                              formfac::Parameterisation* ff = nullptr,
                              strfun::Parameterisation* sf = nullptr) const;

    /// Compute the scalar kT-dependent flux
    double ktFlux(double x,
                  double q2,
                  double mx2 = -1.,
                  formfac::Parameterisation* ff = nullptr,
                  strfun::Parameterisation* sf = nullptr) const;

    /// Compute the flux for a given parton \f$(x,k_{\rm T})\f$
    /// \param[in] type Flux modelling
    /// \param[in] x Parton momentum fraction
    /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
    /// \param[in] ff Form factors evaluator
    /// \param[in] sf Structure functions evaluator
    /// \param[in] mi2 Incoming particle squared mass
    /// \param[in] mf2 Outgoing diffractive squared mass
    static double ktFluxNucl(const KTFlux& type,
                             double x,
                             double kt2,
                             formfac::Parameterisation* ff = nullptr,
                             strfun::Parameterisation* sf = nullptr,
                             double mi2 = -1.,
                             double mf2 = -1.);
    /// Compute the flux (from heavy ion) for a given parton \f$(x,k_{\rm T})\f$
    /// \param[in] type Flux modelling
    /// \param[in] x Parton momentum fraction
    /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
    /// \param[in] hi Heavy ion properties
    static double ktFluxHI(const KTFlux& type, double x, double kt2, const HeavyIon& hi);

  private:
    static const double kMinKTFlux;  ///< Minimal value taken for a \f$\k_{\rm T}\f$-factorised flux
    pdgid_t pdg_;                    ///< PDG identifier for the beam
    Momentum momentum_;              ///< Incoming particle momentum
    Mode mode_;                      ///< Beam treatment mode
    KTFlux kt_flux_;                 ///< Type of \f$k_{\rm T}\f$-factorised flux to be considered (if any)
  };
}  // namespace cepgen

#endif
