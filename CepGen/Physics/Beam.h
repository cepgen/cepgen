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
#include "CepGen/Physics/Momentum.h"

namespace cepgen {
  enum class KTFlux;
  namespace strfun {
    class Parameterisation;
  }
  namespace formfac {
    class Parameterisation;
  }

  /// Incoming beams characteristics
  class Beam : public SteeredObject<Beam> {
  public:
    explicit Beam(const ParametersList& params = ParametersList());  ///< Default constructor

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
    };
    /// Human-readable format of a beam mode (elastic/dissociative parts)
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    const Mode& mode() const { return mode_; }
    /// Is the beam particle expected to be fragmented after emission?
    bool fragmented() const;

    pdgid_t pdgId() const { return pdg_; }
    const Momentum& momentum() const { return momentum_; }
    const KTFlux& ktFlux() const { return kt_flux_; }

    struct FormFactors {
      double FE, FM;
    };
    /// Compute the electromagnetic form factors
    FormFactors flux(double q2,
                     double mx2,
                     formfac::Parameterisation* ff = nullptr,
                     strfun::Parameterisation* sf = nullptr) const;
    /// Compute the scalar kT-dependent flux
    double ktFlux(double x,
                  double q2,
                  double mx2 = -1.,
                  formfac::Parameterisation* ff = nullptr,
                  strfun::Parameterisation* sf = nullptr) const;

  private:
    pdgid_t pdg_;        ///< PDG identifier for the beam
    Momentum momentum_;  ///< Incoming particle momentum
    Mode mode_;          ///< Beam treatment mode
    KTFlux kt_flux_;     ///< Type of \f$k_{\rm T}\f$-factorised flux to be considered (if any)
  };
}  // namespace cepgen

#endif
