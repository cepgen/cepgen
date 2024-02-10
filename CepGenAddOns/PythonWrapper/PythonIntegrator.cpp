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
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  class PythonIntegrator final : public Integrator {
  public:
    explicit PythonIntegrator(const ParametersList& params)
        : Integrator(params), env_(ParametersList().setName<std::string>("python_integrator")) {
      auto cfg = python::importModule(steer<std::string>("module"));
      if (!cfg)
        throw PY_ERROR << "Failed to import the Python module '" << steer<std::string>("module") << "'.";
      func_ = python::getAttribute(cfg.get(), "integrate");
      if (!func_ || !PyCallable_Check(func_.get()))
        throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
    }

    void setLimits(const std::vector<Limits>& lims) override { lims_ = python::set(lims); }

    Value integrate(Integrand& integrand) override {
      gIntegrand = &integrand;
      const auto iterations = steer<int>("iterations");
      const auto evals = steer<int>("evals");
      PyMethodDef py_integr = {"integrand", py_integrand, METH_VARARGS, "A python-wrapped integrand"};
      python::ObjectPtr function(PyCFunction_NewEx(&py_integr, nullptr, python::set<std::string>("integrand").get()));
      const auto value =
          lims_ ? python::call(func_, function.get(), (int)integrand.size(), iterations, 1000, evals, lims_.get())
                : python::call(func_, function.get(), (int)integrand.size(), iterations, 1000, evals);
      if (!value)
        throw PY_ERROR;
      const auto vals = python::getVector<double>(value);
      if (vals.size() < 2)
        throw CG_FATAL("PythonIntegrator")
            << "Wrong multiplicity of result returned from Python's integration algorithm: " << vals << ".";

      return Value{vals[0], vals[1]};
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Python integration algorithm");
      desc.add<std::string>("module", "IntegrationAlgos.Vegas")
          .setDescription("name of the Python module embedding the integrate() function");
      desc.add<int>("iterations", 10);
      desc.add<int>("evals", 1000);
      return desc;
    }
    static Integrand* gIntegrand;

  private:
    python::Environment env_;
    python::ObjectPtr func_, lims_{nullptr};
    static PyObject* py_integrand(PyObject* /*self*/, PyObject* args) {
      if (!gIntegrand)
        throw CG_FATAL("PythonIntegrator") << "Integrand was not initialised.";
      const auto c_args = python::getVector<double>(PyTuple_GetItem(args, 0));
      return python::set<double>(gIntegrand->eval(c_args)).release();
    }
  };
  Integrand* PythonIntegrator::gIntegrand = nullptr;
}  // namespace cepgen

REGISTER_INTEGRATOR("python", PythonIntegrator);
