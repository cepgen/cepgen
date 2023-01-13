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

#ifndef CepGen_Physics_PartonFlux_h
#define CepGen_Physics_PartonFlux_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class PartonFlux : public NamedModule<std::string> {
  public:
    explicit PartonFlux(const ParametersList&);

    static ParametersDescription description();

    virtual double operator()(double x, double kt2, double mf2 = 0.) = 0;

    virtual bool ktFactorised() const { return false; }
    virtual bool fragmenting() const { return true; }
    virtual int partonPdgId() const;

  protected:
    const double mp_, mp2_;
  };
}  // namespace cepgen

#endif
