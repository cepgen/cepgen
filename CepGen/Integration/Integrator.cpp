/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"

namespace cepgen {
  Integrator::Integrator(const ParametersList& params)
      : NamedModule(params),
        seed_(params.get<int>("seed", time(nullptr))),
        verbosity_(params.get<int>("verbose", 1)),
        result_(0.),
        err_result_(0.),
        initialised_(false),
        rnd_(0., 1.) {}

  void Integrator::setIntegrand(Integrand& integr) {
    integrand_ = &integr;

    //--- force the reinitialisation
    initialised_ = false;
  }

  //------------------------------------------------------------------------------------------------
  // helper / alias methods
  //------------------------------------------------------------------------------------------------

  size_t Integrator::size() const {
    if (!integrand_)
      throw CG_FATAL("Integrator:size") << "Trying to retrieve phase space size on an unitialised integrand!";
    return integrand_->size();
  }

  double Integrator::eval(const std::vector<double>& x) const {
    if (!integrand_)
      throw CG_FATAL("Integrator:eval") << "Trying to evaluate the weight on a phase space point "
                                        << "on an unitialised integrand!";
    return integrand_->eval(x);
  }

  double Integrator::uniform() const { return rnd_(rnd_gen_); }
}  // namespace cepgen
