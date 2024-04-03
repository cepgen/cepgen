/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#ifndef CepGenAddOns_CubaWrapper_CubaIntegrator_h
#define CepGenAddOns_CubaWrapper_CubaIntegrator_h

#include "CepGen/Integration/Integrator.h"

namespace cepgen {
  class Integrand;
  /// Interface objects to Cuba algorithms
  namespace cuba {
    /// Cuba integration algorithm
    class Integrator : public cepgen::Integrator {
    public:
      explicit Integrator(const ParametersList&);

      static ParametersDescription description();
      static Integrand* gIntegrand;

      Value integrate(Integrand&) override;

    protected:
      virtual Value integrate() = 0;

      int ncomp_, nvec_;
      double epsrel_, epsabs_;
      int mineval_, maxeval_;
      int verbose_;
    };

    int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/);
  }  // namespace cuba
}  // namespace cepgen

#endif
