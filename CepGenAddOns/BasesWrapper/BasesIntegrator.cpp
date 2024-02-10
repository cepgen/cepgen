/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/BasesWrapper/BasesCommonBlocks.h"

namespace cepgen {
  /// Bases integration algorithm
  class BasesIntegrator : public Integrator {
  public:
    explicit BasesIntegrator(const ParametersList& params) : Integrator(params) {
      bsinit_();
      bparm1_.ncall = steer<int>("numFunctionCalls");
      std::fill(bparm1_.ig.begin(), bparm1_.ig.end(), false);
      bscntl_.intv = steer<int>("intv");
      bscntl_.ipnt = steer<int>("verbose");
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Bases integration algorithm");
      desc.add<int>("numFunctionCalls", 50'000);
      desc.add<int>("intv", 1);
      desc.add<int>("verbose", 0);
      desc.add<std::vector<int> >("wildVars", {}).setDescription("list of 'wild' variables");
      return desc;
    }

    void setLimits(const std::vector<Limits>& lims) override {
      Integrator::setLimits(lims);
      for (size_t i = 0; i < limits_.size(); ++i) {
        bparm1_.xl[i] = limits_.at(i).min();
        bparm1_.xu[i] = limits_.at(i).max();
      }
    }

    Value integrate(Integrand& integr) override {
      checkLimits(integr);  // check the integration bounds
      bparm1_.ndim = integr.size();
      const auto wild_vars = steer<std::vector<int> >("wildVars");
      bparm1_.nwild = wild_vars.size();
      for (const auto& wc : wild_vars) {
        if (wc < 0 || wc >= bparm1_.ndim)
          throw CG_FATAL("BasesIntegrator:integrate") << "Invalid 'wild' variable coordinate set: " << wc << ".";
        bparm1_.ig[wc] = true;
      }
      double res, unc, ctime;
      int it1, it2;
      gIntegrand = &integr;
      bases_(integrand_bases, res, unc, ctime, it1, it2);
      CG_DEBUG("BasesIntegrator:integrate")
          << "Integration performed in " << ctime << " s. " << utils::s("iteration", it1, true)
          << " for the grid definition, " << utils::s("iteration", it2, true) << " for the integration.";
      return Value{res, unc};
    }

  private:
    static Integrand* gIntegrand;
    static double integrand_bases(double in[]) {
      if (!gIntegrand)
        throw CG_FATAL("BasesIntegrator") << "Integrand was not specified before integration.";
      return gIntegrand->eval(std::vector<double>(in, in + gIntegrand->size()));
    }
  };
  Integrand* BasesIntegrator::gIntegrand = nullptr;
}  // namespace cepgen
REGISTER_INTEGRATOR("bases", BasesIntegrator);
