/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2026  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
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

#ifndef CepGenProcesses_PPtoVV_h
#define CepGenProcesses_PPtoVV_h

#include <cmath>
#include <complex>
#include <memory>
#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/HelicityAmplitudes/PolarisationState.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  /// Base class for the \f$k_{\rm T}\f$-factorised \f$\gamma\gamma\rightarrow VV\f$ processes
  class PPtoVV : public proc::FactorisedProcess {
  public:
    static ParametersDescription description() {
      auto desc = proc::FactorisedProcess::description();
      desc.add("ktFactorised", true);
      desc.add("polarisationStates", PolarisationState::description());
      return desc;
    }

  protected:
    /// \param[in] params steering parameters
    /// \param[in] central central system particles
    /// \param[in] boson_pdg PDG id of the produced bosons (both are assumed to have the same mass)
    /// \param[in] symmetry_factor multiplicative factor applied on the squared amplitudes (1/2 for identical bosons)
    explicit PPtoVV(const ParametersList& params,
                    const spdgids_t& central,
                    pdgid_t boson_pdg,
                    double symmetry_factor = 1.)
        : proc::FactorisedProcess(params, central),
          mv_(PDG::get().mass(boson_pdg)),
          mv2_(mv_ * mv_),
          symmetry_factor_(symmetry_factor),
          polarisation_(polarisationParameters(params_)) {}

    /// \f$k_{\rm T}\f$-factorised matrix element, summed over the requested polarisation states of the two bosons
    double offShellME() const {
      const auto kin = amplitudes::Kinematics(mv2_, shat(), that(), uhat());
      const double p1 = q1().px() * q2().px() + q1().py() * q2().py(),
                   p2 = q1().px() * q2().py() - q1().py() * q2().px(),
                   p3 = q1().px() * q2().px() - q1().py() * q2().py(),
                   p4 = q1().px() * q2().py() + q1().py() * q2().px();

      using namespace std::complex_literals;
      double hel_mat_elem{0.};
      for (const auto& lam3 : polarisation_.polarisations().first)
        for (const auto& lam4 : polarisation_.polarisations().second) {
          // compute all photon helicity amplitudes
          const auto pp = amplitudes_->amplitude(kin, {+1, +1, lam3, lam4}),
                     mm = amplitudes_->amplitude(kin, {-1, -1, lam3, lam4}),
                     pm = amplitudes_->amplitude(kin, {+1, -1, lam3, lam4}),
                     mp = amplitudes_->amplitude(kin, {-1, +1, lam3, lam4});
          // add ME for this boson helicity to total ME
          hel_mat_elem += std::norm(p1 * (pp + mm) - 1i * p2 * (pp - mm) - p3 * (pm + mp) - 1i * p4 * (pm - mp));
        }
      return symmetry_factor_ * hel_mat_elem * std::pow(0.5 / q1().pt() / q2().pt() / shat(), 2);
    }

    const double mv_;   ///< boson mass (in GeV/c^2)
    const double mv2_;  ///< squared boson mass (in GeV^2/c^4)
    const double symmetry_factor_;
    const PolarisationState polarisation_;
    std::unique_ptr<amplitudes::HelicityAmplitudes> amplitudes_;

  private:
    /// Build the polarisation states steering, also accepting the integer mode used by steering cards
    static ParametersList polarisationParameters(const ParametersList& params) {
      if (params.has<int>("polarisationStates"))
        return ParametersList().set<int>("mode", params.get<int>("polarisationStates"));
      return params.get<ParametersList>("polarisationStates");
    }
  };
}  // namespace cepgen

#endif
