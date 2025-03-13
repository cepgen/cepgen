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

#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Physics/PDG.h"

using namespace cepgen;
using namespace cepgen::formfac;

Parameterisation::Parameterisation(const ParametersList& params)
    : NamedModule(params),
      pdg_id_(steer<pdgid_t>("pdgId")),
      mass2_(std::pow(PDG::get().mass(pdg_id_), 2)),
      mp_(PDG::get().mass(PDG::proton)),
      mp2_(mp_ * mp_) {}

double Parameterisation::tau(double q2) const { return 0.25 * q2 / mass2_; }

const FormFactors& Parameterisation::operator()(double q2) {
  if (q2 < 0.)
    ff_ = FormFactors{};
  else if (q2 != q2_) {
    q2_ = q2;
    eval();
  }
  return ff_;
}

ParametersDescription Parameterisation::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed form factors parameterisation");
  desc.add<pdgid_t>("pdgId", PDG::proton);
  return desc;
}

void Parameterisation::setFEFM(double fe, double fm) {
  ff_.FE = fe;
  ff_.FM = fm;
  ff_.GM = std::sqrt(ff_.FM);
  const auto ta = tau(q2_);
  ff_.GE = std::sqrt((1. + ta) * ff_.FE - ta * ff_.FM);
}

void Parameterisation::setGEGM(double ge, double gm) {
  ff_.GE = ge;
  ff_.GM = gm;
  ff_.FM = ff_.GM * ff_.GM;
  const auto ta = tau(q2_);
  ff_.FE = (ff_.GE * ff_.GE + ta * ff_.FM) / (1. + ta);
}

//------------------------------------------------------------------

namespace cepgen::formfac {
  std::ostream& operator<<(std::ostream& os, const Parameterisation& ff) {
    os << ff.name();
    if (ff.q2_ >= 0.)
      os << "(Q^2=" << ff.q2_ << " GeV^2): " << ff.ff_;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const FormFactors& ff) {
    return os << "FE=" << ff.FE << ", FM=" << ff.FM << ", GE=" << ff.GE << ", GM=" << ff.GM;
  }
}  // namespace cepgen::formfac
