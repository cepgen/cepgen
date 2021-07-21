#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include <cmath>

namespace cepgen {
  namespace formfac {
    /// \cite Mergell:1995bf
    class MergellEtAl : public Parameterisation {
    public:
      explicit MergellEtAl(const ParametersList&);
      static std::string description() { return "Mergell et al."; }

    private:
      void compute(double q2) override;

      const double a1rho_, a2rho_, b1rho_, b2rho_, c1rho_, c2rho_, d1rho_, d2rho_;
      const double inv_q20_;
      const double lambda_sq_;
      const double gamma_;
    };

    MergellEtAl::MergellEtAl(const ParametersList& params)
        : Parameterisation(params),
          a1rho_(params.get<double>("a1rho", 1.0317)),
          a2rho_(params.get<double>("a2rho", 5.7824)),
          b1rho_(params.get<double>("b1rho", 0.0875)),
          b2rho_(params.get<double>("b2rho", 0.3907)),
          c1rho_(params.get<double>("c1rho", 0.3176)),
          c2rho_(params.get<double>("c2rho", 0.1422)),
          d1rho_(params.get<double>("d1rho", 0.5496)),
          d2rho_(params.get<double>("d2rho", 0.5362)),
          inv_q20_(params.get<double>("q20inv", 1. / 0.35)),
          lambda_sq_(params.get<double>("lambdaSq", 9.733)),
          gamma_(params.get<double>("gamma", 2.148)) {}

    void MergellEtAl::compute(double q2) {
      const double log1 = std::pow(log((lambda_sq_ + q2) * inv_q20_), -gamma_);  // L(t=-q2) function in ref.

      // best fit parameterisation
      const double d1_1 = 0.611 + q2, d2_1 = 1.039 + q2, d3_1 = 2.560 + q2;
      const double Fs1 = (9.464 / d1_1 - 9.054 / d2_1 - 0.410 / d3_1) * log1;
      const double Fs2 = (-1.549 / d1_1 + 1.985 / d2_1 - 0.436 / d3_1) * log1;

      const double log2 = std::pow(log((lambda_sq_ - 0.500) * inv_q20_), +gamma_);
      const double log3 = std::pow(log((lambda_sq_ - 0.400) * inv_q20_), +gamma_);

      const double d1_2 = 2.103 + q2, d2_2 = 2.734 + q2, d3_2 = 2.835 + q2;
      const double Fv1 = (0.5 * (a1rho_ * log2 + b1rho_ * log3 * std::pow(1. + q2 / c1rho_, -2)) / (1. + q2 / d1rho_) -
                          38.885 / d1_2 + 425.007 / d2_2 - 389.742 / d3_2) *
                         log1;
      const double Fv2 = (0.5 * (a2rho_ * log2 + b2rho_ * log3 / (1. + q2 / c2rho_)) / (1. + q2 / d2rho_) -
                          73.535 / d1_2 + 83.211 / d2_2 - 29.467 / d3_2) *
                         log1;

      const double F1 = Fv1 + Fs1, F2 = Fv2 + Fs2;

      GE = F1 - tau(q2) * F2;
      GM = F1 + F2;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL("Mergell", MergellEtAl)
