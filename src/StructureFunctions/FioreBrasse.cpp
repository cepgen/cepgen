/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::strfun {
  ///\f${\cal W}_{1,2}\f$ structure functions parameterisation by Fiore et al \cite Fiore:2002re and Brasse et al. \cite Brasse:1976bf
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
      desc.add("s0", 1.14);
      desc.add("norm", 0.021).setDescription("absolute normalisation factor");
      // add the list of resonances
      desc.addParametersDescriptionVector("resonances",
                                          Resonance::description(),
                                          {ParametersList()
                                               .set("alpha0", -0.8377)
                                               .set("alpha1", 0.95)
                                               .set("alpha2", 0.1473)
                                               .set("a", 1.0)
                                               .set("q02", 2.4617)
                                               .set("spinTimesTwo", 3) /* N*(1520) */,
                                           ParametersList()
                                               .set("alpha0", -0.37)
                                               .set("alpha1", 0.95)
                                               .set("alpha2", 0.1471)
                                               .set("a", 0.5399)
                                               .set("q02", 2.4617)
                                               .set("spinTimesTwo", 5) /* N*(1680) */,
                                           ParametersList()
                                               .set("alpha0", 0.0038)
                                               .set("alpha1", 0.85)
                                               .set("alpha2", 0.1969)
                                               .set("a", 4.2225)
                                               .set("q02", 1.5722)
                                               .set("spinTimesTwo", 3) /* Δ(1236) */,
                                           ParametersList()
                                               .set("alpha0", 0.5645)
                                               .set("alpha1", 0.1126)
                                               .set("alpha2", 1.3086)
                                               .set("a", 19.2694)
                                               .set("q02", 4.5259)
                                               .set("spinTimesTwo", 2) /* exotic */})
          .setDescription("collection of resonances parameters");
      return desc;
    }

    void eval() override;

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
        desc.add("alpha0", 0.);
        desc.add("alpha1", 0.);
        desc.add("alpha2", 0.);
        desc.add("a", 0.).setDescription("resonance weight in total amplitude");
        desc.add("q02", 0.);
        desc.add("spinTimesTwo", 0).setDescription("spin of the resonance (x1/2)");
        return desc;
      }

      double alpha0, alpha1, alpha2, a, q02;
      int spinTimesTwo;
    };

  private:
    /// All resonances considered in this modelling
    std::vector<Resonance> resonances_;
    double s0_{0.}, norm_{0.};
    static constexpr double constant_prefactor_ = 0.25 * M_1_PI / constants::ALPHA_EM;
  };

  void FioreBrasse::eval() {
    const double prefactor = constant_prefactor_ * args_.q2 * (1. - args_.xbj) / gamma2(args_.xbj, args_.q2);
    const double s = utils::mX2(args_.xbj, args_.q2, mp2_);

    double amplitude_res = 0.;
    const double sqrts0 = std::sqrt(s0_);
    for (unsigned short i = 0; i < 3; ++i) {  //FIXME 4??
      const auto& res = resonances_.at(i);

      std::complex<double> alpha;
      if (s > s0_)
        alpha =
            std::complex<double>(res.alpha0 + res.alpha2 * sqrts0 + res.alpha1 * s, res.alpha2 * std::sqrt(s - s0_));
      else
        alpha = std::complex<double>(res.alpha0 + res.alpha1 * s + res.alpha2 * (sqrts0 - std::sqrt(s0_ - s)), 0.);

      double formfactor = std::pow(1. + args_.q2 / res.q02, -2);
      double denominator = std::pow(res.spinTimesTwo * 0.5 - std::real(alpha), 2) + std::pow(std::imag(alpha), 2);
      double imaginary_amplitude = res.a * formfactor * formfactor * std::imag(alpha) / denominator;
      amplitude_res += imaginary_amplitude;
    }
    double amplitude_bg = 0.;
    {
      const auto& res = resonances_.at(3);
      double sE = res.alpha2, sqrtsE = std::sqrt(sE);
      std::complex<double> alpha;
      if (s > sE)
        alpha = std::complex<double>(res.alpha0 + res.alpha1 * sqrtsE, res.alpha1 * std::sqrt(s - sE));
      else
        alpha = std::complex<double>(res.alpha0 + res.alpha1 * (sqrtsE - std::sqrt(sE - s)), 0.);
      double form_factor = std::pow(1. + args_.q2 / res.q02, -2);
      double denominator = std::pow(res.spinTimesTwo * 0.75 - std::real(alpha), 2) + std::pow(std::imag(alpha), 2);
      amplitude_bg = res.a * form_factor * form_factor * std::imag(alpha) / denominator;
    }
    const double amplitude_tot = norm_ * (amplitude_res + amplitude_bg);

    CG_DEBUG_LOOP("FioreBrasse:amplitudes") << "Amplitudes:\n\t"
                                            << " resonance part:  " << amplitude_res << ",\n\t"
                                            << " background part: " << amplitude_bg << ",\n\t"
                                            << " total (with norm.): " << amplitude_tot << ".";

    setF2(prefactor * amplitude_tot);
  }

  /// Alternative Fiore-Brasse parameterisation
  class FioreBrasseAlt final : public FioreBrasse {
  public:
    explicit FioreBrasseAlt(const ParametersList& params) : FioreBrasse(params) {}

    static ParametersDescription description() {
      auto desc = FioreBrasse::description();
      desc.add("s0", 1.2871);
      desc.add("norm", 0.0207);
      // add the list of resonances
      desc.addParametersDescriptionVector("resonances",
                                          Resonance::description(),
                                          {ParametersList()
                                               .set("alpha0", -0.8070)
                                               .set("alpha1", 0.9632)
                                               .set("alpha2", 0.1387)
                                               .set("a", 1.0)
                                               .set("q02", 2.6066)
                                               .set("spinTimesTwo", 3) /* N*(1520) */,
                                           ParametersList()
                                               .set("alpha0", -0.3640)
                                               .set("alpha1", 0.9531)
                                               .set("alpha2", 0.1239)
                                               .set("a", 0.6086)
                                               .set("q02", 2.6066)
                                               .set("spinTimesTwo", 5) /* N*(1680) */,
                                           ParametersList()
                                               .set("alpha0", -0.0065)
                                               .set("alpha1", 0.8355)
                                               .set("alpha2", 0.2320)
                                               .set("a", 4.7279)
                                               .set("q02", 1.4828)
                                               .set("spinTimesTwo", 3) /* Δ(1236) */,
                                           ParametersList()
                                               .set("alpha0", 0.5484)
                                               .set("alpha1", 0.1373)
                                               .set("alpha2", 1.3139)
                                               .set("a", 14.7267)
                                               .set("q02", 4.6041)
                                               .set("spinTimesTwo", 2) /* exotic */});
      return desc;
    }
  };
}  // namespace cepgen::strfun
using cepgen::strfun::FioreBrasse;
using cepgen::strfun::FioreBrasseAlt;
REGISTER_STRFUN("FioreBrasse", 101, FioreBrasse);
REGISTER_STRFUN("FioreBrasseAlt", 104, FioreBrasseAlt);
