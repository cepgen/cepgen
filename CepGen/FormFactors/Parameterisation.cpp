/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace formfac {
    Parameterisation::Parameterisation()
        : NamedModule<std::string>(ParametersList()), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {}

    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule<std::string>(params), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {}

    Parameterisation::Parameterisation(const Parameterisation& param)
        : NamedModule<std::string>(param.parameters()),
          mp_(param.mp_),
          mp2_(param.mp2_),
          last_q2_(-1.),
          FE(param.FE),
          FM(param.FM),
          GE(param.GE),
          GM(param.GM) {}

    double Parameterisation::tau(double q2) const {
      if (mp2_ <= 0.)
        throw CG_FATAL("FormFactors:tau") << "Invalid proton mass! check the form factors constructor!";
      return 0.25 * q2 / mp2_;
    }

    Parameterisation& Parameterisation::operator()(const Beam::Mode& type,
                                                   double q2,
                                                   double mf2,
                                                   strfun::Parameterisation* sf) {
      last_q2_ = q2;
      switch (type) {
        case Beam::Mode::invalid:
        case Beam::Mode::CompositeScalar:
        case Beam::Mode::Other:
          throw CG_FATAL("FormFactors") << type << " mode is not yet supported!";
        case Beam::Mode::PointLikeScalar:
          FE = 1., FM = 0.;
          break;
        case Beam::Mode::PointLikeFermion:
          FE = FM = 1.;  // FE=U2, FM=U1 in LPAIR
          break;
        case Beam::Mode::ProtonElastic: {
          compute(q2);
          const double GE2 = GE * GE, GM2 = GM * GM;
          FE = (4. * mp2_ * GE2 + q2 * GM2) / (4. * mp2_ + q2);
          FM = GM2;
        } break;
        case Beam::Mode::ProtonInelastic: {
          if (!sf)
            throw CG_FATAL("FormFactors")
                << "Inelastic proton form factors computation requires a structure functions definition!";
          const double xbj = utils::xBj(q2, mp2_, mf2);
          switch ((strfun::Type)sf->name()) {
            case strfun::Type::ElasticProton:
              throw CG_FATAL("FormFactors") << "Elastic proton form factors requested!\n"
                                            << "Check your process definition!";
            case strfun::Type::SuriYennie: {  // this one requires its own object to deal with FM
              FE = sf->F2(xbj, q2) * xbj * mp_ / q2;
              FM = sf->FM(xbj, q2);
            } break;
            default: {
              FE = sf->F2(xbj, q2) * xbj / q2;
              FM = -2. * sf->F1(xbj, q2) / q2;
            } break;
          }
        } break;
      }
      return *this;
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed form factors parameterisation");
      return desc;
    }

    //------------------------------------------------------------------

    std::ostream& operator<<(std::ostream& os, const Parameterisation* ff) {
      if (!ff)
        return os << "[uninitialised form factors]";
      os << ff->name();
      if (ff->last_q2_ >= 0.)
        os << "(Q²=" << ff->last_q2_ << " GeV²): "
           << "FE=" << ff->FE << ",FM=" << ff->FM;
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& ff) { return os << &ff; }
  }  // namespace formfac
}  // namespace cepgen
