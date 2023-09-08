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
  class IntegratorVegasPlus final : public Integrator {
  public:
    explicit IntegratorVegasPlus(const ParametersList& params) : Integrator(params), env_("vegas_plus") {
      auto cfg = python::importModule("VegasIntegration");
      if (!cfg)
        throw PY_ERROR << "Failed to import the Vegas python file.";
      func_ = python::getAttribute(cfg.get(), "integrate");
      if (!func_ || !PyCallable_Check(func_.get()))
        throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
    }

    void integrate(Integrand& integrand, double& result, double& abs_error) override {
      gIntegrand = &integrand;
      const auto iterations = steer<int>("iterations");
      const auto evals = steer<int>("evals");
      PyMethodDef py_integr = {"integrand", py_integrand, METH_VARARGS, "A python-wrapped integrand"};
      python::ObjectPtr function(PyCFunction_NewEx(&py_integr, nullptr, python::set<std::string>("integrand").get()));
      python::ObjectPtr value(PyObject_CallObject(
          func_.get(),
          python::newTuple(std::make_tuple(function.get(), (int)integrand.size(), iterations, 1000, evals))
              .release()));  // new
      if (!value)
        throw PY_ERROR;
      const auto vals = python::getVector<double>(value);
      if (vals.size() < 2)
        throw CG_FATAL("IntegratorVegasPlus")
            << "Wrong multiplicity of result returned from Python's Vegas: " << vals << ".";
      result_ = result = vals[0];
      err_result_ = abs_error = vals[1];
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Vegas+ MC integrator");
      desc.add<int>("iterations", 10);
      desc.add<int>("evals", 1000);
      return desc;
    }
    static Integrand* gIntegrand;

  private:
    python::Environment env_;
    python::ObjectPtr func_;
    static PyObject* py_integrand(PyObject* /*self*/, PyObject* args) {
      if (!gIntegrand)
        throw CG_FATAL("IntegratorVegasPlus") << "Integrand was not initialised.";
      const auto c_args = python::getVector<double>(PyTuple_GetItem(args, 0));
      return python::set<double>(gIntegrand->eval(c_args)).release();
    }
  };
  Integrand* IntegratorVegasPlus::gIntegrand = nullptr;
}  // namespace cepgen

REGISTER_INTEGRATOR("vegas_plus", IntegratorVegasPlus);
