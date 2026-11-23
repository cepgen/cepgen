/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2026  Laurent Forthomme
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

#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"

using namespace cepgen;
using namespace cepgen::amplitudes;

HelicityAmplitudes::HelicityAmplitudes(const ParametersList& params) : NamedModule(params), s1_(steer<double>("s1")) {}

ParametersDescription HelicityAmplitudes::description() {
  auto desc = ParametersDescription();
  desc.setDescription("generic helicity model");
  desc.add("s1", 0.);
  return desc;
}

std::complex<double> HelicityAmplitudes::amplitudeLinearPolarisations(const Kinematics& kin,
                                                                      short lam3,
                                                                      short lam4,
                                                                      const std::complex<double>& e1,
                                                                      const std::complex<double>& e2) const {
  const auto pp = amplitude(kin, {+1, +1, lam3, lam4}), mm = amplitude(kin, {-1, -1, lam3, lam4}),
             pm = amplitude(kin, {+1, -1, lam3, lam4}), mp = amplitude(kin, {-1, +1, lam3, lam4});
  return 0.5 * (e1 * std::conj(e2) * pp + std::conj(e1) * e2 * mm - e1 * e2 * pm - std::conj(e1 * e2) * mp);
}

double HelicityAmplitudes::s1() const { return s1_; }

double HelicityAmplitudes::c1() const { return std::sqrt(1. - s1_ * s1_); }
