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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace formfac {
    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule<std::string>(params),
          pdg_id_(steer<pdgid_t>("pdgId")),
          mass2_(std::pow(HeavyIon::isHI(pdg_id_) ? HeavyIon::fromPdgId(pdg_id_).mass() : PDG::get().mass(pdg_id_), 2)),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_) {}

    double Parameterisation::tau(double q2) const {
      if (mp2_ <= 0.)
        throw CG_FATAL("FormFactors:tau") << "Invalid proton mass! check the form factors constructor!";
      return 0.25 * q2 / mp2_;
    }

    const FormFactors& Parameterisation::operator()(double q2) {
      q2_ = q2;
      compute();
      return last_ff_;
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed form factors parameterisation");
      desc.add<pdgid_t>("pdgId", PDG::invalid);
      return desc;
    }

    void Parameterisation::setFEFM(double fe, double fm) {
      last_ff_.FE = fe;
      last_ff_.FM = fm;
      const double tau = 0.25 * q2_ / mass2_;
      last_ff_.GM = std::sqrt(last_ff_.FM);
      last_ff_.GE = std::sqrt((1. + tau) * last_ff_.FE - tau * last_ff_.FM);
    }

    void Parameterisation::setGEGM(double ge, double gm) {
      last_ff_.GE = ge;
      last_ff_.GM = gm;
      last_ff_.FM = last_ff_.GM * last_ff_.GM;
      last_ff_.FE = (4. * mass2_ * last_ff_.GE * last_ff_.GE + q2_ * last_ff_.FM) / (4. * mass2_ + q2_);
    }

    //------------------------------------------------------------------

    std::ostream& operator<<(std::ostream& os, const Parameterisation* ff) {
      if (!ff)
        return os << "[uninitialised form factors]";
      os << ff->name();
      if (ff->q2_ >= 0.)
        os << "(Q²=" << ff->q2_ << " GeV²): " << ff->last_ff_;
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& ff) { return os << &ff; }
    std::ostream& operator<<(std::ostream& os, const FormFactors& ff) {
      return os << "FE=" << ff.FE << ", FM=" << ff.FM << ", GE=" << ff.GE << ", GM=" << ff.GM;
    }
  }  // namespace formfac
}  // namespace cepgen
