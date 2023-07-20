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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGenAddOns/BoostWrapper/PythonUtils.h"

namespace {
  BOOST_PYTHON_MODULE(pycepgen) {
    cepgen::initialise();

    py::class_<cepgen::Steerable>("_Steerable", "base steerable object", py::no_init)
        .add_property(
            "parameters",
            +[](const cepgen::Steerable& st) { return plist_to_py_dict(st.parameters()); },
            &cepgen::strfun::Parameterisation::setParameters,
            "Operational parameters")
        .add_property(
            "name",
            +[](const cepgen::Steerable& st) { return st.parameters().getString(cepgen::MODULE_NAME); },
            "Module name");

    py::class_<cepgen::strfun::Parameterisation, py::bases<cepgen::Steerable>, boost::noncopyable>(
        "_StructureFunctions", "nucleon structure functions modelling", py::no_init)
        .def("F2", &cepgen::strfun::Parameterisation::F2)
        .def("FL", &cepgen::strfun::Parameterisation::FL)
        .def("F1", &cepgen::strfun::Parameterisation::F1);

    EXPOSE_FACTORY(cepgen::StructureFunctionsFactory,
                   int,
                   "StructureFunctionsFactory",
                   "a structure functions evaluator objects factory");

    py::class_<cepgen::sigrat::Parameterisation, py::bases<cepgen::Steerable>, boost::noncopyable>(
        "_SigmaRatio", "L/T cross section ratio modelling", py::no_init)
        .def(
            "__call__", +[](const cepgen::sigrat::Parameterisation& par, double xbj, double q2) {
              double unc{0.};
              const auto sig_rat = par(xbj, q2, unc);
              return std_vector_to_py_tuple(std::vector<double>{sig_rat, unc});
            });

    EXPOSE_FACTORY(cepgen::SigmaRatiosFactory,
                   int,
                   "SigmaRatiosFactory",
                   "a longitudinal-to-transverse cross section ratio evaluator objects factory");

    py::class_<cepgen::formfac::Parameterisation, py::bases<cepgen::Steerable>, boost::noncopyable>(
        "_FormFactors", "nucleon electromagnetic form factors modelling", py::no_init)
        .def("__call__",
             py::make_function(&cepgen::formfac::Parameterisation::operator(), py::return_internal_reference()));

    py::class_<cepgen::formfac::FormFactors>("FormFactors", "nucleon electromagnetic form factors values")
        .add_property("FE", &cepgen::formfac::FormFactors::FE, "Electric form factor")
        .add_property("FM", &cepgen::formfac::FormFactors::FM, "Magnetic form factor")
        .add_property("GE", &cepgen::formfac::FormFactors::GE, "Sachs electric form factor")
        .add_property("GM", &cepgen::formfac::FormFactors::GM, "Sachs magnetic form factor");

    EXPOSE_FACTORY(cepgen::FormFactorsFactory,
                   std::string,
                   "FormFactorsFactory",
                   "an electromagnetic form factors evaluator objects factory");
  }
}  // namespace
