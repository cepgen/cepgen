/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#ifndef CepGen_Process_PartonSampler_h
#define CepGen_Process_PartonSampler_h

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/PartonFluxes/PartonFlux.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class Beam;
  class Kinematics;
  class KTFlux;
  class CollinearFlux;
}  // namespace cepgen
namespace cepgen::cuts {
  struct Initial;
}  // namespace cepgen::cuts

namespace cepgen {
  /// Incoming parton kinematics sampler for a single beam
  /// \note A sampler owns the parton flux of one beam, and the integration variables this flux needs
  ///  (e.g. parton transverse momentum and azimuthal angle for kT-factorised fluxes, or parton virtuality for
  ///  virtuality-unintegrated collinear fluxes). It is evaluated in two steps:
  ///  1. partonMomentum() gives the parton 4-momentum seed (transverse part), before the momentum fraction is known;
  ///  2. weight(x, mx2) gives the flux (and its own variables' Jacobian) once x and the remnant mass are known.
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Nov 2026
  class PartonSampler : public SteeredObject<PartonSampler> {
  public:
    explicit PartonSampler(const ParametersList&);

    static ParametersDescription description();

    /// Beam side this sampler is attached to
    enum class Side { positive = 0, negative = 1 };

    /// Description of one integration variable required by this sampler
    struct IntegrationVariable {
      std::reference_wrapper<double> value;  ///< Reference to the variable to be mapped
      proc::Process::Mapping mapping;        ///< Mapping type
      Limits range;                          ///< Integration range (in mapped coordinates)
      std::string name;                      ///< Computer-readable variable name
      std::string description;               ///< Human-readable variable description
    };

    /// Build the parton flux for one beam (user-steered, or default for this beam type)
    void initialise(const Kinematics&, Side);

    virtual bool ktFactorised() const = 0;  ///< Does the parton carry a primordial kT?
    /// List of integration variables required by this sampler (to be registered by the consumer process)
    virtual std::vector<IntegrationVariable> integrationVariables(const cuts::Initial&) = 0;
    virtual Momentum partonMomentum() const = 0;  ///< Parton 4-momentum seed for the current phase space point
    /// Flux weight (including this sampler's own variables Jacobian) once x and the remnant mass are known
    virtual double weight(double x, double mx2) const = 0;

    const PartonFlux& flux() const;  ///< Parton flux modelling
    /// Type-casted parton flux modelling
    template <typename T>
    inline const T& flux() const {
      return dynamic_cast<const T&>(flux());
    }

    inline Side side() const { return side_; }  ///< Beam side this sampler is attached to

  protected:
    /// Beam categories used to pick a default parton flux
    enum class BeamCategory { elasticPointLike, elasticHadron, inelasticHadron, elasticHeavyIon };
    static BeamCategory beamCategory(const Beam&);
    /// Default flux parameters for this sampler's flux family, and a given beam
    virtual ParametersDescription defaultFluxParameters(BeamCategory, const Beam&, const Kinematics&) const = 0;
    /// Build (and type-check) a flux object of this sampler's family
    virtual std::unique_ptr<PartonFlux> buildFlux(const ParametersList&) const = 0;
    virtual void fluxInitialised() {}  ///< Hook called once the flux object is built
    /// Side-dependent variable name/description helpers
    std::string variableName(const std::string& base) const;
    std::string sideLabel() const;

    const bool log_parton_virtuality_;  ///< Generate linearly in log(Q^2) instead of Q^2?
    Side side_{Side::positive};
    std::unique_ptr<PartonFlux> flux_;
  };

  /// Sampler for kT-factorised (unintegrated) parton fluxes: two variables, parton qT and its azimuthal angle
  class KTPartonSampler final : public PartonSampler {
  public:
    using PartonSampler::PartonSampler;

    bool ktFactorised() const override { return true; }
    std::vector<IntegrationVariable> integrationVariables(const cuts::Initial&) override;
    Momentum partonMomentum() const override;
    double weight(double x, double mx2) const override;

  private:
    ParametersDescription defaultFluxParameters(BeamCategory, const Beam&, const Kinematics&) const override;
    std::unique_ptr<PartonFlux> buildFlux(const ParametersList&) const override;
    void fluxInitialised() override;

    const KTFlux* kt_flux_{nullptr};  ///< Type-casted flux (NOT owning)
    double qt_{0.};                   ///< Parton transverse virtuality
    double phi_qt_{0.};               ///< Parton transverse virtuality azimuthal angle
  };

  /// Sampler for collinear parton fluxes: one variable, the parton virtuality
  class CollinearPartonSampler final : public PartonSampler {
  public:
    using PartonSampler::PartonSampler;

    bool ktFactorised() const override { return false; }
    std::vector<IntegrationVariable> integrationVariables(const cuts::Initial&) override;
    Momentum partonMomentum() const override;
    double weight(double x, double mx2) const override;

  private:
    ParametersDescription defaultFluxParameters(BeamCategory, const Beam&, const Kinematics&) const override;
    std::unique_ptr<PartonFlux> buildFlux(const ParametersList&) const override;
    void fluxInitialised() override;

    const CollinearFlux* coll_flux_{nullptr};  ///< Type-casted flux (NOT owning)
    double t_{0.};                             ///< Parton virtuality
  };
}  // namespace cepgen

#endif
