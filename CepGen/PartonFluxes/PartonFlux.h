/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  class PartonFlux : public NamedModule<PartonFlux> {
  public:
    explicit PartonFlux(const ParametersList&);

    static ParametersDescription description();

    virtual bool ktFactorised() const { return false; }  ///< Is the flux parton kT-dependent?
    virtual bool fragmenting() const = 0;                ///< Is initiator particle fragmenting after parton emission?
    virtual pdgid_t partonPdgId() const = 0;             ///< Parton PDG identifier
    virtual double mass2() const = 0;                    ///< Initiator particle squared mass (in \f${\rm GeV}^2/c^4\f$)

  protected:
    const double alpha_over_pi_;
    const double mp_, mp2_;
    const Limits x_range_{0., 1.};
  };
}  // namespace cepgen

#endif
