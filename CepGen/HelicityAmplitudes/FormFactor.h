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

#ifndef CepGen_HelicityAmplitudes_FormFactor_h
#define CepGen_HelicityAmplitudes_FormFactor_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen::amplitudes {
  /// \brief Form factor taming the growth of anomalous couplings with the two-photon invariant mass \f$W=\sqrt{\hat s}\f$
  /// \note It is meant to multiply the anomalous couplings (not the squared matrix element), such that the interference
  ///   with the SM amplitude scales as \f$F\f$ and the pure anomalous contribution as \f$F^2\f$.
  struct FormFactor : NamedModule<FormFactor> {
    using NamedModule<FormFactor>::NamedModule;
    static std::unique_ptr<FormFactor> build(const ParametersList&);  ///< Build a form factor modelling

    static ParametersDescription description();

    /// Form factor value at two-parton squared invariant mass (in GeV^2/c^4)
    virtual double operator()(double /*shat*/) const { return 1.; }
  };
}  // namespace cepgen::amplitudes

#endif
