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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/HelicityAmplitudes/FormFactor.h"

using namespace cepgen;
using namespace cepgen::amplitudes;

/// Dipole form factor, as described for \f$n=2\f$ in arXiv:0807.1121
/// \note Each anomalous coupling is replaced by \f$a\to a/(1+W^2/\Lambda^2)^n\f$.
class DipoleFormFactor : public FormFactor {
public:
  explicit DipoleFormFactor(const ParametersList& params)
      : FormFactor(params), scale2_(std::pow(steer<double>("scale"), 2)), n_(steer<double>("exponent")) {
    if (scale2_ == 0.)
      throw CG_FATAL("DipoleFormFactor")
          << "Strictly positive cut-off scale (in GeV) is required for dipole form factor.";
    if (n_ <= 0.)
      throw CG_FATAL("DipoleFormFactor") << "Dipole exponent must be strictly positive, got " << n_ << ".";
  }

  static ParametersDescription description() {
    auto desc = FormFactor::description();
    desc.setDescription("Dipole (1+W^2/scale^2)^-n form factor");
    desc.add("scale", 500.).setDescription("cut-off scale, in GeV");
    desc.add("exponent", 2.).setDescription("dipole exponent");
    return desc;
  }

  double operator()(double shat) const override { return std::pow(1. + shat / scale2_, -n_); }

private:
  const double scale2_;
  const double n_;
};

class SharpFormFactor : public FormFactor {
public:
  explicit SharpFormFactor(const ParametersList& params)
      : FormFactor(params), scale2_(std::pow(steer<double>("scale"), 2)) {
    if (scale2_ == 0.)
      throw CG_FATAL("SharpFormFactor")
          << "Strictly positive cut-off scale (in GeV) is required for sharp form factor.";
  }

  static ParametersDescription description() {
    auto desc = FormFactor::description();
    desc.setDescription("Dipole (1+W^2/scale^2)^-n form factor");
    desc.add("scale", 500.).setDescription("cut-off scale, in GeV");
    return desc;
  }

  double operator()(double shat) const override { return shat > scale2_ ? 0. : 1.; }

private:
  const double scale2_;
};

ParametersDescription FormFactor::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Trivial form factor");
  return desc;
}

std::unique_ptr<FormFactor> FormFactor::build(const ParametersList& params) {
  if (const auto& name = params.name(); name == "none")
    return std::make_unique<FormFactor>(params);
  else if (name == "dipole")
    return std::make_unique<DipoleFormFactor>(params);
  else if (name == "sharp")
    return std::make_unique<SharpFormFactor>(params);
  else
    throw CG_FATAL("FormFactor:build") << "Failed to build a taming form factor with name '" << name << "'.";
}
