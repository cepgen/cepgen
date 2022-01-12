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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"

namespace cepgen {
  /// General definition for a resonance
  ResonanceObject::ResonanceObject(const ParametersList& params)
      : SteeredObject(params),
        br_(steer<ParametersList>("branchingRatios")),
        ang_mom_(steer<int>("angularMomentum")),
        x0_(steer<double>("x0")),
        mass_(steer<double>("mass")),
        width_(steer<double>("width")),
        mp_(PDG::get().mass(PDG::proton)),
        mp2_(mp_ * mp_),
        mpi2_(std::pow(PDG::get().mass(PDG::piZero), 2)),
        meta2_(std::pow(PDG::get().mass(PDG::eta), 2)),
        x02_(x0_ * x0_) {}

  ParametersDescription ResonanceObject::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Set of physical properties for one resonance");
    desc.add<ParametersDescription>("branchingRatios", BranchingRatios::description());
    desc.add<int>("angularMomentum", 0).setDescription("meson angular momentum");
    desc.add<double>("x0", 0.).setDescription("damping parameter");
    desc.add<double>("mass", 0.).setDescription("mass, in GeV/c^2");
    desc.add<double>("width", 0.).setDescription("full width, in GeV");
    return desc;
  }

  double ResonanceObject::ecmr(double m2) const { return mass_ == 0 ? 0. : utils::energyFromW(mass_, mp2_, m2); }

  ResonanceObject::BranchingRatios::BranchingRatios(const ParametersList& params)
      : SteeredObject(params),
        singlepi(steer<double>("singlePi")),
        doublepi(steer<double>("doublePi")),
        eta(steer<double>("eta")) {
    if (!valid())
      CG_WARNING("ResonanceObject:BranchingRatios")
          << "Invalid branching fractions. Sum = " << (singlepi + doublepi + eta) << " != 1.";
  }

  ParametersDescription ResonanceObject::BranchingRatios::description() {
    auto desc = ParametersDescription();
    desc.add<double>("singlePi", 0.).setDescription("branching fraction for a resonance decay into a single pion");
    desc.add<double>("doublePi", 0.).setDescription("branching fraction for a resonance decay into a pion pair");
    desc.add<double>("eta", 0.).setDescription("branching fraction for a resonance decay into an eta");
    return desc;
  }
}  // namespace cepgen
