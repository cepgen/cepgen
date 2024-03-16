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
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGenAddOns/BoostWrapper/PythonObjectsWrappers.h"
#include "CepGenAddOns/BoostWrapper/PythonUtils.h"

namespace {
  /// Python interfacing module definition
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
        .add_static_property("name",
                             py::make_function(&cepgen::strfun::Parameterisation::name,
                                               py::return_value_policy<py::copy_const_reference>()))
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

    py::class_<PartonFluxWrap, py::bases<cepgen::Steerable>, boost::noncopyable>(
        "_PartonFlux", "generic parton flux evaluator", py::no_init)
        .add_property("partonPdgId", &cepgen::PartonFlux::partonPdgId)
        .add_property("fragmenting", &cepgen::PartonFlux::fragmenting)
        .add_property("ktFactorised", &cepgen::PartonFlux::ktFactorised)
        .def(
            "__call__",
            +[](const cepgen::PartonFlux& flux) {
              if (flux.ktFactorised())
                return adapt_reference(&dynamic_cast<const cepgen::KTFlux&>(flux));
              return adapt_reference(&dynamic_cast<const cepgen::CollinearFlux&>(flux));
            },
            "Expose the flux evaluator object from its type");

    py::class_<CollinearFluxWrap, py::bases<PartonFluxWrap>, boost::noncopyable>(
        "_CollinearFlux", "fractional momentum/virtuality-dependent parton flux evaluator", py::no_init)
        .def("fluxMX2", py::pure_virtual(&cepgen::CollinearFlux::fluxMX2))
        .def("fluxQ2", py::pure_virtual(&cepgen::CollinearFlux::fluxQ2));

    EXPOSE_FACTORY(cepgen::CollinearFluxFactory,
                   std::string,
                   "CollinearFluxFactory",
                   "a collinear parton fluxes evaluator objects factory");

    py::class_<KTFluxWrap, py::bases<PartonFluxWrap>, boost::noncopyable>(
        "_KTFlux",
        "fractional momentum/(transverse-longitudinal) virtuality-dependent parton flux evaluator",
        py::no_init)
        .def("fluxMX2", py::pure_virtual(&cepgen::KTFlux::fluxMX2))
        .def("fluxQ2", py::pure_virtual(&cepgen::KTFlux::fluxQ2));

    EXPOSE_FACTORY(
        cepgen::KTFluxFactory, std::string, "KTFluxFactory", "a kt-factorised parton fluxes evaluator objects factory");

    py::class_<cepgen::PDG, boost::noncopyable>("PDG", "collection of particle definitions and properties", py::no_init)
        .def(
            "colours", +[](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().colours(pdgid); })
        .def(
            "mass", +[](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().mass(pdgid); })
        .def(
            "width", +[](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().width(pdgid); })
        .def(
            "charge", +[](cepgen::pdgid_t pdgid) { return cepgen::PDG::get().charge(pdgid); });
  }
}  // namespace
