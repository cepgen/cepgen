/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PolarisationState.h"

namespace cepgen {
  PolarisationState::PolarisationState(const ParametersList& params)
      : SteeredObject(params),
        mode_(params_.has<int>("mode") ? steerAs<int, Mode>("mode") : Mode::invalid),
        pol_(mode_ == Mode::invalid ? std::make_pair(steer<std::vector<int> >("W1"), steer<std::vector<int> >("W2"))
                                    : computePolarisations(mode_)) {}

  PolarisationState::Polarisations PolarisationState::computePolarisations(const Mode& mode) {
    switch (mode) {
      case Mode::LL:
        return std::make_pair(Polarisation{0}, Polarisation{0});
      case Mode::LT:
        return std::make_pair(Polarisation{0}, Polarisation{-1, 1});
      case Mode::TL:
        return std::make_pair(Polarisation{-1, 1}, Polarisation{0});
      case Mode::TT:
        return std::make_pair(Polarisation{-1, 1}, Polarisation{-1, 1});
      case Mode::full:
        return std::make_pair(Polarisation{-1, 0, 1}, Polarisation{-1, 0, 1});
      default:
        throw CG_FATAL("PolarisationState:computePolarisations")
            << "Invalid mode for polarisation states computation: " << mode << ".";
    }
  }

  ParametersDescription PolarisationState::description() {
    auto desc = ParametersDescription();
    desc.addAs<int, Mode>("mode", Mode::invalid);
    desc.add<std::vector<int> >("W1", {-1, 0, 1}).setDescription("First polarisation states");
    desc.add<std::vector<int> >("W2", {-1, 0, 1}).setDescription("Second polarisation states");
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const PolarisationState::Mode& mode) {
    switch (mode) {
      case PolarisationState::Mode::invalid:
      default:
        return os << "invalid";
      case PolarisationState::Mode::full:
        return os << "full";
      case PolarisationState::Mode::LT:
        return os << "long-trans";
      case PolarisationState::Mode::TL:
        return os << "trans-long";
      case PolarisationState::Mode::TT:
        return os << "trans-trans";
      case PolarisationState::Mode::LL:
        return os << "long-long";
    }
  }
}  // namespace cepgen
