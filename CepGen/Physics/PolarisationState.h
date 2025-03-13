/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGen_Physics_PolarisationState_h
#define CepGen_Physics_PolarisationState_h

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  class PolarisationState final : public SteeredObject<PolarisationState> {
  public:
    explicit PolarisationState(const ParametersList&);

    static ParametersDescription description();

    enum class Mode { invalid = -1, full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    typedef std::vector<int> Polarisation;
    typedef std::pair<Polarisation, Polarisation> Polarisations;

    const Mode& mode() const { return mode_; }
    const Polarisations& polarisations() const { return pol_; }

  private:
    static Polarisations computePolarisations(const Mode&);

    const Mode mode_;
    const Polarisations pol_;
  };
}  // namespace cepgen

#endif
