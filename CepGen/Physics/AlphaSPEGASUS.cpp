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
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"

namespace {
  extern "C" {
  void initalphas_(int& iord, double& fr2, double& mur, double& asmur, double& mc, double& mb, double& mt);
  double alphas_(double& mur);
  }
}  // namespace

namespace cepgen {
  class AlphaSPEGASUS : public Coupling {
  public:
    explicit AlphaSPEGASUS(const ParametersList& params)
        : Coupling(params),
          iord_(steer<int>("iord")),
          fr2_(steer<double>("fr2")),
          mur_(steer<double>("mur")),
          asmur_(steer<double>("asmur")) {
      double mc = PDG::get().mass(4), mb = PDG::get().mass(5), mt = PDG::get().mass(6);

      initalphas_(iord_, fr2_, mur_, asmur_, mc, mb, mt);
      CG_INFO("AlphaSPEGASUS:init") << "PEGASUS alpha(S) evolution algorithm initialised with parameters:\n\t"
                                    << "order: " << iord_ << ", fr2: " << fr2_ << ", "
                                    << "mur: " << mur_ << ", asmur: " << asmur_ << "\n\t"
                                    << "quark masses (GeV): charm: " << mc << ", bottom: " << mb << ", top: " << mt
                                    << ".";
    }

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("PEGASUS alpha(S) evolution algorithm");
      desc.add<int>("iord", 2).setDescription("Evolution order");
      desc.add<double>("fr2", 1.);
      desc.add<double>("mur", 1.);
      desc.add<double>("asmur", 0.68183);
      return desc;
    }

    double operator()(double q) const override { return alphas_(q); }

  private:
    int iord_;
    double fr2_;
    double mur_;
    double asmur_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("pegasus", AlphaSPEGASUS)
