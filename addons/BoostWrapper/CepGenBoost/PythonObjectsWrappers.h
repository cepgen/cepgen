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

#ifndef CepGenBoost_PythonObjectsWrappers_h
#define CepGenBoost_PythonObjectsWrappers_h

#include <boost/python/wrapper.hpp>

#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/PartonFluxes/KTFlux.h"

namespace py = boost::python;

struct PartonFluxWrap : cepgen::PartonFlux, py::wrapper<cepgen::PartonFlux> {
  bool ktFactorised() const override {
    if (py::override ov = this->get_override("ktFactorised"))
      return ov();
    return PartonFlux::ktFactorised();
  }
  bool fragmenting() const override {
    if (py::override ov = this->get_override("fragmenting"))
      return ov();
    return PartonFlux::fragmenting();
  }
  cepgen::pdgid_t partonPdgId() const override {
    if (py::override ov = this->get_override("partonPdgId"))
      return ov();
    return PartonFlux::partonPdgId();
  }
  //double mass2() const override { return this->get_override("mass2")(); }
};

struct KTFluxWrap : cepgen::KTFlux, py::wrapper<cepgen::KTFlux> {
  double fluxMX2(double x, double kt2, double mx2) const override {
    if (py::override ov = this->get_override("fluxMX2"))
      return ov(x, kt2, mx2);
    return KTFlux::fluxMX2(x, kt2, mx2);
  }
  double fluxQ2(double x, double kt2, double q2) const override {
    if (py::override ov = this->get_override("fluxQ2"))
      return ov(x, kt2, q2);
    return KTFlux::fluxQ2(x, kt2, q2);
  }
};

struct CollinearFluxWrap : cepgen::CollinearFlux, py::wrapper<cepgen::CollinearFlux> {
  double fluxMX2(double x, double mx2) const override {
    if (py::override ov = this->get_override("fluxMX2"))
      return ov(x, mx2);
    return CollinearFlux::fluxMX2(x, mx2);
  }
  double fluxQ2(double x, double q2) const override {
    if (py::override ov = this->get_override("fluxQ2"))
      return ov(x, q2);
    return CollinearFlux::fluxQ2(x, q2);
  }
};

#endif
