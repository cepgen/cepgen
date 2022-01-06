/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

namespace cepgen {
  namespace strfun {
    Parameterisation::Parameterisation(double f2, double fl)
        : NamedModule<int>(ParametersList()),
          F2(f2),
          FL(fl),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          mx_min_(mp_ + PDG::get().mass(PDG::piZero)),
          r_ratio_(sigrat::SigmaRatiosFactory::get().build((int)sigrat::Type::SibirtsevBlunden)) {}

    Parameterisation::Parameterisation(const Parameterisation& sf)
        : NamedModule<int>(sf.parameters()),
          F2(sf.F2),
          FL(sf.FL),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          mx_min_(mp_ + PDG::get().mass(PDG::piPlus)),
          old_vals_(sf.old_vals_),
          r_ratio_(sf.r_ratio_) {}

    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule<int>(params),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          mx_min_(mp_ + PDG::get().mass(PDG::piZero)) {
      CG_DEBUG("Parameterisation") << "Structure functions parameterisation to be built using following parameters:\n"
                                   << ParametersDescription(params_).describe(true);
      r_ratio_ = sigrat::SigmaRatiosFactory::get().build(steer<int>("sigmaRatio"));
    }

    Parameterisation& Parameterisation::operator=(const Parameterisation& sf) {
      F2 = sf.F2, FL = sf.FL;
      old_vals_ = sf.old_vals_;
      return *this;
    }

    Parameterisation& Parameterisation::operator()(double xbj, double q2) {
      if (old_vals_.first == xbj && old_vals_.second == q2)
        return *this;
      old_vals_ = {xbj, q2};
      return eval(xbj, q2);
    }

    Parameterisation& Parameterisation::eval(double, double) {
      CG_WARNING("StructureFunctions") << "Evaluation method called on base object!";
      return *this;
    }

    double Parameterisation::tau(double xbj, double q2) const { return 4. * xbj * xbj * mp2_ / q2; }

    double Parameterisation::F1(double xbj, double q2) const {
      if (xbj == 0. || q2 == 0.) {
        CG_WARNING("StructureFunctions:F1") << "Invalid range for Q² = " << q2 << " or xBj = " << xbj << ".";
        return 0.;
      }
      const double F1 = 0.5 * ((1 + tau(xbj, q2)) * F2 - FL) / xbj;
      CG_DEBUG_LOOP("StructureFunctions:F1") << "F1 for Q² = " << q2 << ", xBj = " << xbj << ": " << F1 << "\n\t"
                                             << "(F2 = " << F2 << ", FL = " << FL << ").";
      return F1;
    }

    Parameterisation& Parameterisation::computeFL(double xbj, double q2) {
      if (!r_ratio_)
        throw CG_FATAL("StructureFunctions:FL") << "Failed to retrieve a R-ratio calculator!";
      double r_error;
      return computeFL(xbj, q2, (*r_ratio_)(xbj, q2, r_error));
    }

    Parameterisation& Parameterisation::computeFL(double xbj, double q2, double r) {
      FL = F2 * (1. + tau(xbj, q2)) * (r / (1. + r));
      return *this;
    }

    std::string Parameterisation::describe() const {
      std::ostringstream os;
      os << (Type)name_;
      return os.str();
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation* sf) {
      os << sf->describe();
      if (sf->old_vals_ != std::pair<double, double>{0., 0.})
        os << " at (" << sf->old_vals_.first << ", " << sf->old_vals_.second << "): "
           << "F2 = " << sf->F2 << ", FL = " << sf->FL;
      return os;
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed structure functions parameterisation");
      desc.add<int>("sigmaRatio", (int)sigrat::Type::SibirtsevBlunden)
          .setDescription("Modelling for the sigma(L/T) ratio used in FL computation from F2");
      return desc;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& sf) { return os << &sf; }

    /// Human-readable format of a structure function type
    std::ostream& operator<<(std::ostream& os, const strfun::Type& sf) {
      switch (sf) {
        case strfun::Type::Invalid:
          return os << "<invalid>";
        case strfun::Type::Electron:
          return os << "Electron";
        case strfun::Type::ElasticProton:
          return os << "ElasticProton";
        case strfun::Type::SuriYennie:
          return os << "SuriYennie";
        case strfun::Type::SuriYennieAlt:
          return os << "SuriYennieAlt";
        case strfun::Type::SzczurekUleshchenko:
          return os << "SzczurekUleshchenko";
        case strfun::Type::FioreBrasse:
          return os << "FioreBrasse";
        case strfun::Type::FioreBrasseAlt:
          return os << "FioreBrasseAlt";
        case strfun::Type::ChristyBosted:
          return os << "ChristyBosted";
        case strfun::Type::CLAS:
          return os << "CLAS";
        case strfun::Type::BlockDurandHa:
          return os << "BlockDurandHa";
        case strfun::Type::ALLM91:
          return os << "ALLM91";
        case strfun::Type::ALLM97:
          return os << "ALLM97";
        case strfun::Type::HHT_ALLM:
          return os << "ALLM{HHT}";
        case strfun::Type::HHT_ALLM_FT:
          return os << "ALLM{HHT_FT}";
        case strfun::Type::GD07p:
          return os << "GD07p";
        case strfun::Type::GD11p:
          return os << "GD11p";
        case strfun::Type::Schaefer:
          return os << "LUXlike";
        case strfun::Type::Shamov:
          return os << "Shamov";
        case strfun::Type::MSTWgrid:
          return os << "MSTWgrid";
        case strfun::Type::Partonic:
          return os << "Partonic";
      }
      return os;
    }
  }  // namespace strfun
}  // namespace cepgen
