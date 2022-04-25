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

#include <complex>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    ///\f${\cal W}_{1,2}\f$ structure functions parameterisation by Fiore et al \cite Fiore:2002re and Brasse et al \cite Brasse:1976bf
    class FioreBrasse : public Parameterisation {
    public:
      /// Fiore \cite Fiore:2002re and Brasse \cite Brasse:1976bf proton structure functions
      explicit FioreBrasse(const ParametersList& params)
          : Parameterisation(params), s0_(steer<double>("s0")), norm_(steer<double>("norm")) {
        for (const auto& res : steer<std::vector<ParametersList> >("resonances"))
          resonances_.emplace_back(res);
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Fiore-Brasse (low-mass resonances)");
        desc.add<double>("s0", 1.14);
        desc.add<double>("norm", 0.021).setDescription("absolute normalisation factor");
        // add the list of resonances
        desc.addParametersDescriptionVector("resonances",
                                            Resonance::description(),
                                            {ParametersList()  // N*(1520)
                                                 .set<double>("alpha0", -0.8377)
                                                 .set<double>("alpha1", 0.95)
                                                 .set<double>("alpha2", 0.1473)
                                                 .set<double>("a", 1.0)
                                                 .set<double>("q02", 2.4617)
                                                 .set<int>("spinTimesTwo", 3),
                                             ParametersList()  // N*(1680)
                                                 .set<double>("alpha0", -0.37)
                                                 .set<double>("alpha1", 0.95)
                                                 .set<double>("alpha2", 0.1471)
                                                 .set<double>("a", 0.5399)
                                                 .set<double>("q02", 2.4617)
                                                 .set<int>("spinTimesTwo", 5),
                                             ParametersList()  // Δ(1236)
                                                 .set<double>("alpha0", 0.0038)
                                                 .set<double>("alpha1", 0.85)
                                                 .set<double>("alpha2", 0.1969)
                                                 .set<double>("a", 4.2225)
                                                 .set<double>("q02", 1.5722)
                                                 .set<int>("spinTimesTwo", 3),
                                             ParametersList()  // exotic
                                                 .set<double>("alpha0", 0.5645)
                                                 .set<double>("alpha1", 0.1126)
                                                 .set<double>("alpha2", 1.3086)
                                                 .set<double>("a", 19.2694)
                                                 .set<double>("q02", 4.5259)
                                                 .set<int>("spinTimesTwo", 2)})
            .setDescription("collection of resonances parameters");
        return desc;
      }

      FioreBrasse& eval(double xbj, double q2) override;

    protected:
      /// Description of a single resonance in the modelling
      struct Resonance : SteeredObject<Resonance> {
        explicit Resonance(const ParametersList& params)
            : SteeredObject(params),
              alpha0(steer<double>("alpha0")),
              alpha1(steer<double>("alpha1")),
              alpha2(steer<double>("alpha2")),
              a(steer<double>("a")),
              q02(steer<double>("q02")),
              spinTimesTwo(steer<int>("spinTimesTwo")) {}

        static ParametersDescription description() {
          auto desc = ParametersDescription();
          desc.add<double>("alpha0", 0.);
          desc.add<double>("alpha1", 0.);
          desc.add<double>("alpha2", 0.);
          desc.add<double>("a", 0.).setDescription("resonance weight in total amplitude");
          desc.add<double>("q02", 0.);
          desc.add<int>("spinTimesTwo", 0).setDescription("spin of the resonance (x1/2)");
          return desc;
        }

        double alpha0, alpha1, alpha2, a, q02;
        int spinTimesTwo;
      };

    private:
      /// All resonances considered in this modelling
      std::vector<Resonance> resonances_;
      double s0_{0.}, norm_{0.};
    };

    FioreBrasse& FioreBrasse::eval(double xbj, double q2) {
      const double akin = 1. + 4. * mp2_ * xbj * xbj / q2;
      const double prefactor = q2 * (1. - xbj) / (4. * M_PI * constants::ALPHA_EM * akin);
      const double s = utils::mX2(xbj, q2, mp2_);

      double amplitude_res = 0.;
      const double sqrts0 = sqrt(s0_);
      for (unsigned short i = 0; i < 3; ++i) {  //FIXME 4??
        const auto& res = resonances_.at(i);

        std::complex<double> alpha;
        if (s > s0_)
          alpha = std::complex<double>(res.alpha0 + res.alpha2 * sqrts0 + res.alpha1 * s, res.alpha2 * sqrt(s - s0_));
        else
          alpha = std::complex<double>(res.alpha0 + res.alpha1 * s + res.alpha2 * (sqrts0 - sqrt(s0_ - s)), 0.);

        double formfactor = 1. / pow(1. + q2 / res.q02, 2);
        double denom = pow(res.spinTimesTwo * 0.5 - std::real(alpha), 2) + pow(std::imag(alpha), 2);
        double ampli_imag = res.a * formfactor * formfactor * std::imag(alpha) / denom;
        amplitude_res += ampli_imag;
      }
      double amplitude_bg = 0.;
      {
        const auto& res = resonances_.at(3);
        double sE = res.alpha2, sqrtsE = sqrt(sE);
        std::complex<double> alpha;
        if (s > sE)
          alpha = std::complex<double>(res.alpha0 + res.alpha1 * sqrtsE, res.alpha1 * sqrt(s - sE));
        else
          alpha = std::complex<double>(res.alpha0 + res.alpha1 * (sqrtsE - sqrt(sE - s)), 0.);
        double formfactor = 1. / pow(1. + q2 / res.q02, 2);
        double denom = pow(res.spinTimesTwo * 0.75 - std::real(alpha), 2) + pow(std::imag(alpha), 2);
        amplitude_bg = res.a * formfactor * formfactor * std::imag(alpha) / denom;
      }
      const double amplitude_tot = norm_ * (amplitude_res + amplitude_bg);

      CG_DEBUG_LOOP("FioreBrasse:amplitudes") << "Amplitudes:\n\t"
                                              << " resonance part:  " << amplitude_res << ",\n\t"
                                              << " background part: " << amplitude_bg << ",\n\t"
                                              << " total (with norm.): " << amplitude_tot << ".";

      setF2(prefactor * amplitude_tot);
      return *this;
    }

    class FioreBrasseAlt final : public FioreBrasse {
    public:
      explicit FioreBrasseAlt(const ParametersList& params) : FioreBrasse(params) {}

      static ParametersDescription description() {
        auto desc = FioreBrasse::description();
        desc.add<double>("s0", 1.2871);
        desc.add<double>("norm", 0.0207);
        // add the list of resonances
        desc.addParametersDescriptionVector("resonances",
                                            Resonance::description(),
                                            {ParametersList()  // N*(1520)
                                                 .set<double>("alpha0", -0.8070)
                                                 .set<double>("alpha1", 0.9632)
                                                 .set<double>("alpha2", 0.1387)
                                                 .set<double>("a", 1.0)
                                                 .set<double>("q02", 2.6066)
                                                 .set<int>("spinTimesTwo", 3),
                                             ParametersList()  // N*(1680)
                                                 .set<double>("alpha0", -0.3640)
                                                 .set<double>("alpha1", 0.9531)
                                                 .set<double>("alpha2", 0.1239)
                                                 .set<double>("a", 0.6086)
                                                 .set<double>("q02", 2.6066)
                                                 .set<int>("spinTimesTwo", 5),
                                             ParametersList()  // Δ(1236)
                                                 .set<double>("alpha0", -0.0065)
                                                 .set<double>("alpha1", 0.8355)
                                                 .set<double>("alpha2", 0.2320)
                                                 .set<double>("a", 4.7279)
                                                 .set<double>("q02", 1.4828)
                                                 .set<int>("spinTimesTwo", 3),
                                             ParametersList()  // exotic
                                                 .set<double>("alpha0", 0.5484)
                                                 .set<double>("alpha1", 0.1373)
                                                 .set<double>("alpha2", 1.3139)
                                                 .set<double>("a", 14.7267)
                                                 .set<double>("q02", 4.6041)
                                                 .set<int>("spinTimesTwo", 2)});
        return desc;
      }
    };
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::FioreBrasse, FioreBrasse, strfun::FioreBrasse)
REGISTER_STRFUN(strfun::Type::FioreBrasseAlt, FioreBrasseAlt, strfun::FioreBrasseAlt)
