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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/String.h"
#include "CepGenPython/Environment.h"
#include "CepGenPython/Error.h"

using namespace std::string_literals;

namespace cepgen::python {
  class Integrator final : public cepgen::Integrator {
  public:
    explicit Integrator(const ParametersList& params)
        : cepgen::Integrator(params), env_(ParametersList().setName("python_integrator")) {
      if (const auto cfg = ObjectPtr::importModule(steer<std::string>("module")); cfg) {
        if (func_ = cfg.attribute("integrate"); !func_ || !PyCallable_Check(func_.get()))
          throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
      } else
        throw PY_ERROR << "Failed to import the Python module '" << steer<std::string>("module") << "'.";
    }

    Value run(Integrand& integrand, const std::vector<Limits>& range) override {
      lims_ = ObjectPtr::make(range);
      gIntegrand = &integrand;
      const auto iterations = steer<int>("iterations");
      const auto evals = steer<int>("evals");
      PyMethodDef python_integrand = {"integrand", py_integrand, METH_VARARGS, "A python-wrapped integrand"};
      const ObjectPtr function(
          PyCFunction_NewEx(&python_integrand, nullptr, ObjectPtr::make<std::string>("integrand").get()));
      const auto value =
          lims_ ? func_(function.get(), static_cast<int>(integrand.size()), iterations, 1000, evals, lims_.get())
                : func_(function.get(), static_cast<int>(integrand.size()), iterations, 1000, evals);
      if (!value)
        throw PY_ERROR;
      const auto vals = value.vector<double>();
      if (vals.size() < 2)
        throw CG_FATAL("python:Integrator")
            << "Wrong multiplicity of result returned from Python's integration algorithm: " << vals << ".";

      return Value{vals[0], vals[1]};
    }

    static ParametersDescription description() {
      auto desc = cepgen::Integrator::description();
      desc.setDescription("Python integration algorithm");
      desc.add("module", "IntegrationAlgos.Vegas"s)
          .setDescription("name of the Python module embedding the integrate() function");
      desc.add("iterations", 10);
      desc.add("evals", 1000);
      return desc;
    }
    static Integrand* gIntegrand;

  private:
    Environment env_;
    ObjectPtr func_{nullptr}, lims_{nullptr};
    static PyObject* py_integrand(PyObject* /*self*/, PyObject* args) {
      if (!gIntegrand)
        throw CG_FATAL("python:Integrator") << "Integrand was not initialised.";
      const auto c_args = ObjectPtr::wrap(PyTuple_GetItem(args, 0)).vector<double>();
      return ObjectPtr::make<double>(gIntegrand->eval(c_args)).release();
    }
  };
  Integrand* Integrator::gIntegrand = nullptr;
}  // namespace cepgen::python
using PythonIntegrator = cepgen::python::Integrator;
REGISTER_INTEGRATOR("python", PythonIntegrator);
