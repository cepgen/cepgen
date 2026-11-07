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

#ifndef CepGen_Physics_SurvivalFactor_h
#define CepGen_Physics_SurvivalFactor_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
}  // namespace cepgen

namespace cepgen {
  /// A generic (soft) gap survival factor \f$S^2\f$ modelling
  /// \note Evaluated once per phase space point at the integrand level, on the event content filled by the process,
  ///  before any event modification algorithm. Both the cross section and the (unweighted) events distributions
  ///  account for the rescaling.
  /// \date Nov 2026
  class SurvivalFactor : public NamedModule<SurvivalFactor> {
  public:
    explicit SurvivalFactor(const ParametersList& params) : NamedModule(params) {}
    ~SurvivalFactor() override = default;

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed gap survival factor modelling");
      return desc;
    }

    /// Compute the survival factor for a given event
    /// \return A value in [0, 1] (values outside this range are not prevented, but will be warned about once)
    virtual double operator()(const Event&) const = 0;
  };
}  // namespace cepgen

#endif
