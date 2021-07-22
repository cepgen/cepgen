#include <complex>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Physics.h"

namespace cepgen {
  namespace strfun {
    ///\f${\cal W}_{1,2}\f$ structure functions parameterisation by Fiore et al \cite Fiore:2002re and Brasse et al \cite Brasse:1976bf
    class FioreBrasse final : public Parameterisation {
    public:
      /// Fiore \cite Fiore:2002re and Brasse \cite Brasse:1976bf proton structure functions
      explicit FioreBrasse(const ParametersList& params = ParametersList());
      static std::string description() { return "Fiore-Brasse F2 parameterisation of low-mass resonances"; }

      FioreBrasse& eval(double xbj, double q2) override;

    private:
      /// General parameters for this modelling
      struct Parameters {
        static Parameters standard();
        static Parameters alternative();
        /// Description of a single resonance in the modelling
        struct Resonance {
          double alpha0, alpha1, alpha2, a, q02;
          float spin;
        };
        /// All resonances considered in this modelling
        std::vector<Resonance> resonances;
        double s0 = {0.}, norm = {0.};
      } fb_params_;
    };

    FioreBrasse::Parameters FioreBrasse::Parameters::standard() {
      Parameters p{};
      p.s0 = 1.14;
      p.norm = 0.021;
      p.resonances.emplace_back(Resonance{-0.8377, 0.95, 0.1473, 1.0, 2.4617, 3. / 2.});    // N*(1520)
      p.resonances.emplace_back(Resonance{-0.37, 0.95, 0.1471, 0.5399, 2.4617, 5. / 2.});   // N*(1680)
      p.resonances.emplace_back(Resonance{0.0038, 0.85, 0.1969, 4.2225, 1.5722, 3. / 2.});  // Δ(1236)
      p.resonances.emplace_back(Resonance{0.5645, 0.1126, 1.3086, 19.2694, 4.5259, 1.});    // exotic
      return p;
    }
    FioreBrasse::Parameters FioreBrasse::Parameters::alternative() {
      Parameters p{};
      p.s0 = 1.2871;
      p.norm = 0.0207;
      p.resonances.emplace_back(Resonance{-0.8070, 0.9632, 0.1387, 1.0, 2.6066, 3. / 2.});     // N*(1520)
      p.resonances.emplace_back(Resonance{-0.3640, 0.9531, 0.1239, 0.6086, 2.6066, 5. / 2.});  // N*(1680)
      p.resonances.emplace_back(Resonance{-0.0065, 0.8355, 0.2320, 4.7279, 1.4828, 3. / 2.});  // Δ(1236)
      p.resonances.emplace_back(Resonance{0.5484, 0.1373, 1.3139, 14.7267, 4.6041, 1.});       // exotic
      return p;
    }

    FioreBrasse::FioreBrasse(const ParametersList& params) : Parameterisation(params) {
      const auto& model = params.get<std::string>("model", "standard");
      if (model == "standard")
        fb_params_ = Parameters::standard();
      else if (model == "alternative")
        fb_params_ = Parameters::alternative();
      else
        throw CG_FATAL("FioreBrasse") << "Invalid modelling selected: " << model << "!";
    }

    FioreBrasse& FioreBrasse::eval(double xbj, double q2) {
      const double akin = 1. + 4. * mp2_ * xbj * xbj / q2;
      const double prefactor = q2 * (1. - xbj) / (4. * M_PI * constants::ALPHA_EM * akin);
      const double s = utils::mX2(xbj, q2, mp2_);

      double amplitude_res = 0.;
      for (unsigned short i = 0; i < 3; ++i) {  //FIXME 4??
        const Parameters::Resonance& res = fb_params_.resonances[i];
        const double sqrts0 = sqrt(fb_params_.s0);

        std::complex<double> alpha;
        if (s > fb_params_.s0)
          alpha = std::complex<double>(res.alpha0 + res.alpha2 * sqrts0 + res.alpha1 * s,
                                       res.alpha2 * sqrt(s - fb_params_.s0));
        else
          alpha =
              std::complex<double>(res.alpha0 + res.alpha1 * s + res.alpha2 * (sqrts0 - sqrt(fb_params_.s0 - s)), 0.);

        double formfactor = 1. / pow(1. + q2 / res.q02, 2);
        double denom = pow(res.spin - std::real(alpha), 2) + pow(std::imag(alpha), 2);
        double ampli_imag = res.a * formfactor * formfactor * std::imag(alpha) / denom;
        amplitude_res += ampli_imag;
      }
      double amplitude_bg = 0.;
      {
        const Parameters::Resonance& res = fb_params_.resonances[3];
        double sE = res.alpha2, sqrtsE = sqrt(sE);
        std::complex<double> alpha;
        if (s > sE)
          alpha = std::complex<double>(res.alpha0 + res.alpha1 * sqrtsE, res.alpha1 * sqrt(s - sE));
        else
          alpha = std::complex<double>(res.alpha0 + res.alpha1 * (sqrtsE - sqrt(sE - s)), 0.);
        double formfactor = 1. / pow(1. + q2 / res.q02, 2);
        double sp = 1.5 * res.spin;
        double denom = pow(sp - std::real(alpha), 2) + pow(std::imag(alpha), 2);
        amplitude_bg = res.a * formfactor * formfactor * std::imag(alpha) / denom;
      }
      const double amplitude_tot = fb_params_.norm * (amplitude_res + amplitude_bg);

      CG_DEBUG_LOOP("FioreBrasse:amplitudes") << "Amplitudes:\n\t"
                                              << " resonance part:  " << amplitude_res << ",\n\t"
                                              << " background part: " << amplitude_bg << ",\n\t"
                                              << " total (with norm.): " << amplitude_tot << ".";

      F2 = prefactor * amplitude_tot;
      return *this;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(FioreBrasse, strfun::FioreBrasse)
