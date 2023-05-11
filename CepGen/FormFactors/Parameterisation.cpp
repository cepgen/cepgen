/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace formfac {
    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule<std::string>(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("incomingParticle"))),
          mass2_(hi_.mass() * hi_.mass()),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_) {}

    double Parameterisation::tau(double q2) const {
      if (mp2_ <= 0.)
        throw CG_FATAL("FormFactors:tau") << "Invalid proton mass! check the form factors constructor!";
      return 0.25 * q2 / mp2_;
    }

    const FormFactors& Parameterisation::operator()(double q2) {
      last_value_.first = q2;
      auto& ff = last_value_.second;
      ff = compute(q2);
      const double GE2 = ff.GE * ff.GE, GM2 = ff.GM * ff.GM;
      ff.FE = (4. * mass2_ * GE2 + q2 * GM2) / (4. * mass2_ + q2);
      ff.FM = GM2;
      return ff;
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed form factors parameterisation");
      desc.addAs<pdgid_t, HeavyIon>("incomingParticle", HeavyIon::proton());
      return desc;
    }

    //------------------------------------------------------------------

    std::ostream& operator<<(std::ostream& os, const Parameterisation* ff) {
      if (!ff)
        return os << "[uninitialised form factors]";
      os << ff->name();
      if (ff->last_value_.first >= 0.)
        os << "(Q²=" << ff->last_value_.first << " GeV²): "
           << "FE=" << ff->last_value_.second.FE << ",FM=" << ff->last_value_.second.FM;
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& ff) { return os << &ff; }
  }  // namespace formfac
}  // namespace cepgen
