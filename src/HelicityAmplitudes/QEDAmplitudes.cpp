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

#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"

using namespace cepgen;
using namespace cepgen::amplitudes;

/// \file
/// Tree-level QED helicity amplitudes for \f$\gamma\gamma\to X\bar X\f$, with \f$X\f$ a fermion or a charged scalar.
/// All amplitudes are returned for a unit coupling and a unit charge (i.e. \f$e=Q=1\f$): the couplings
/// (running \f$\alpha\f$, charge, colour factors) are to be applied by the process using the amplitudes.
/// Photon helicities are \f$\pm1\f$ (photon 1 along +z, photon 2 along -z), and the amplitudes are given in the
/// scattering plane (azimuthal angle of the outgoing pair set to zero).

///////////////////////////////////////////////////////////////////////////////////////////////////

/// \f$\gamma\gamma\to S^+S^-\f$ in scalar QED
struct ScalarQED final : HelicityAmplitudes {
  explicit ScalarQED(const ParametersList& params) : HelicityAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Scalar QED (γγ → S⁺S⁻)");
    return desc;
  }
  bool hasAASS() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    if (hel[0] == 0 || hel[1] == 0 || hel[2] != 0 || hel[3] != 0)  // transverse photons
      return 0.;                                                   // -> spinless scalars only
    if (hel[0] == hel[1])                                          // J_z = 0
      return 2. * kin.inv_gamma2 * kin.invA;
    return 2. * kin.beta2 * kin.sin_theta2 * kin.invA;  // |J_z| = 2
  }
};
REGISTER_HELICITY_AMPLITUDES("scalarQED", ScalarQED);

///////////////////////////////////////////////////////////////////////////////////////////////////

/// \f$\gamma\gamma\to f\bar f\f$ in QED, with \f$m_f\neq 0\f$
struct FermionQED final : HelicityAmplitudes {
  explicit FermionQED(const ParametersList& params) : HelicityAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Fermion QED (γγ → f f̄)");
    return desc;
  }
  bool hasAAFF() const override { return true; }
  /// \param[in] hel Photons helicities (±1) followed by the fermion and anti-fermion helicities (±1, i.e. 2λ)
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    if (hel[0] == 0 || hel[1] == 0 || hel[2] == 0 || hel[3] == 0)
      return 0.;
    const double lam1 = hel[0], sig = hel[2];
    if (hel[0] == hel[1])  // J_z = 0: only same-helicity fermions
      return hel[2] != hel[3] ? 0. : -2. * kin.inv_gamma * kin.invA * (lam1 - sig * kin.beta);
    if (hel[2] == hel[3])  // |J_z| = 2 along the beam, no spin projection along the fermion axis
      return -2. * sig * kin.beta * kin.inv_gamma * kin.sin_theta2 * kin.invA;
    return -2. * kin.beta * kin.sin_theta * kin.invA * (lam1 * sig + kin.cos_theta);
  }
};
REGISTER_HELICITY_AMPLITUDES("fermionQED", FermionQED);
