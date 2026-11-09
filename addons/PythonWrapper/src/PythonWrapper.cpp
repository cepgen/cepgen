/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>
#include <CepGen/FormFactors/Parameterisation.h>
#include <CepGen/Generator.h>
#include <CepGen/Modules/FormFactorsFactory.h>
#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/Modules/StructureFunctionsFactory.h>
#include <CepGen/PartonFluxes/CollinearFlux.h>
#include <CepGen/PartonFluxes/KTFlux.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/StructureFunctions/Parameterisation.h>
#include <CepGen/StructureFunctions/SigmaRatio.h>
#include <pybind11/stl.h>

#include <sstream>

#include "CepGenPython/PythonUtils.h"

namespace py = pybind11;

namespace {
  /// Singletons are owned by the C++ side, and are only exposed through their static members
  template <typename T>
  using singleton_class = py::class_<T, std::unique_ptr<T, py::nodelete> >;

  template <typename Factory, typename Key>
  singleton_class<Factory> expose_factory(py::module_& module, const char* name, const char* description) {
    singleton_class<Factory> factory(module, name, description);
    factory.def_static("build", [](const Key& mod) { return Factory::get().build(mod); })
        .def_static("build",
                    [](const Key& mod, const py::dict& dict) {
                      return Factory::get().build(mod, cepgen::python::py_dict_to_plist(dict));
                    })
        .def_static("build",
                    [](const py::dict& dict) { return Factory::get().build(cepgen::python::py_dict_to_plist(dict)); })
        .def_static("describe",
                    [](const Key& mod) {
                      std::ostringstream os;
                      os << Factory::get().describeParameters(mod);
                      return os.str();
                    })
        .def_property_readonly_static("modules", [](const py::object&) { return Factory::get().modules(); });
    return factory;
  }
}  // namespace

/// Python interfacing module definition
PYBIND11_MODULE(pycepgen, module) {
  cepgen::initialise();

  // forward the complete message, rather than the fixed-size buffer returned by Exception::what()
  py::register_exception_translator([](std::exception_ptr exception) {
    try {
      if (exception)
        std::rethrow_exception(exception);
    } catch (const cepgen::Exception& error) {
      PyErr_SetString(PyExc_RuntimeError, error.message().c_str());
    }
  });

  py::class_<cepgen::Steerable>(module, "_Steerable", "base steerable object")
      .def_property(
          "parameters",
          [](const cepgen::Steerable& steerable) { return cepgen::python::plist_to_py_dict(steerable.parameters()); },
          [](cepgen::Steerable& steerable, const py::dict& dict) {
            steerable.setParameters(cepgen::python::py_dict_to_plist(dict));
          },
          "Operational parameters")
      .def_property_readonly(
          "name",
          [](const cepgen::Steerable& steerable) { return steerable.parameters().getString(cepgen::MODULE_NAME); },
          "Module name");

  py::class_<cepgen::strfun::Parameterisation, cepgen::Steerable>(
      module, "_StructureFunctions", "nucleon structure functions modelling")
      .def("F2", &cepgen::strfun::Parameterisation::F2)
      .def("FL", &cepgen::strfun::Parameterisation::FL)
      .def("F1", &cepgen::strfun::Parameterisation::F1);

  expose_factory<cepgen::StructureFunctionsFactory, std::string>(
      module, "StructureFunctionsFactory", "a structure functions evaluator objects factory")
      .def_static("build", [](int mod) { return cepgen::StructureFunctionsFactory::get().build(mod); });

  py::class_<cepgen::sigrat::Parameterisation, cepgen::Steerable>(
      module, "_SigmaRatio", "L/T cross section ratio modelling")
      .def("__call__", [](const cepgen::sigrat::Parameterisation& parameterisation, double xbj, double q2) {
        double uncertainty{0.};
        const auto sigma_ratio = parameterisation(xbj, q2, uncertainty);
        return py::make_tuple(sigma_ratio, uncertainty);
      });

  expose_factory<cepgen::SigmaRatiosFactory, int>(
      module, "SigmaRatiosFactory", "a longitudinal-to-transverse cross section ratio evaluator objects factory");

  py::class_<cepgen::formfac::Parameterisation, cepgen::Steerable>(
      module, "_FormFactors", "nucleon electromagnetic form factors modelling")
      .def("__call__", &cepgen::formfac::Parameterisation::operator(), py::return_value_policy::reference_internal);

  py::class_<cepgen::formfac::FormFactors>(module, "FormFactors", "nucleon electromagnetic form factors values")
      .def(py::init<>())
      .def_readonly("FE", &cepgen::formfac::FormFactors::FE, "Electric form factor")
      .def_readonly("FM", &cepgen::formfac::FormFactors::FM, "Magnetic form factor")
      .def_readonly("GE", &cepgen::formfac::FormFactors::GE, "Sachs electric form factor")
      .def_readonly("GM", &cepgen::formfac::FormFactors::GM, "Sachs magnetic form factor");

  expose_factory<cepgen::FormFactorsFactory, std::string>(
      module, "FormFactorsFactory", "an electromagnetic form factors evaluator objects factory");

  py::class_<cepgen::PartonFlux, cepgen::Steerable>(module, "_PartonFlux", "generic parton flux evaluator")
      .def_property_readonly("partonPdgId", &cepgen::PartonFlux::partonPdgId)
      .def_property_readonly("fragmenting", &cepgen::PartonFlux::fragmenting)
      .def_property_readonly("ktFactorised", &cepgen::PartonFlux::ktFactorised)
      .def(
          "__call__",
          [](const cepgen::PartonFlux& flux) -> py::object {
            if (flux.ktFactorised())
              return py::cast(&dynamic_cast<const cepgen::KTFlux&>(flux), py::return_value_policy::reference);
            return py::cast(&dynamic_cast<const cepgen::CollinearFlux&>(flux), py::return_value_policy::reference);
          },
          "Expose the flux evaluator object from its type");

  py::class_<cepgen::CollinearFlux, cepgen::PartonFlux>(
      module, "_CollinearFlux", "fractional momentum/virtuality-dependent parton flux evaluator")
      .def("fluxMX2", &cepgen::CollinearFlux::fluxMX2)
      .def("fluxQ2", &cepgen::CollinearFlux::fluxQ2);

  expose_factory<cepgen::CollinearFluxFactory, std::string>(
      module, "CollinearFluxFactory", "a collinear parton fluxes evaluator objects factory");

  py::class_<cepgen::KTFlux, cepgen::PartonFlux>(
      module, "_KTFlux", "fractional momentum/(transverse-longitudinal) virtuality-dependent parton flux evaluator")
      .def("fluxMX2", &cepgen::KTFlux::fluxMX2)
      .def("fluxQ2", &cepgen::KTFlux::fluxQ2);

  expose_factory<cepgen::KTFluxFactory, std::string>(
      module, "KTFluxFactory", "a kt-factorised parton fluxes evaluator objects factory");

  singleton_class<cepgen::PDG>(module, "PDG", "collection of particle definitions and properties")
      .def_static("colours", [](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().colours(pdgid); })
      .def_static("mass", [](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().mass(pdgid); })
      .def_static("width", [](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().width(pdgid); })
      .def_static("charge", [](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().charge(pdgid); });
}
