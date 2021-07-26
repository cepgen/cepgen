#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

namespace cepgen {
  namespace formfac {
    /// \cite Arrington:2007ux
    class ArringtonEtAl final : public Parameterisation {
    public:
      explicit ArringtonEtAl(const ParametersList&);

      static std::string description() { return "Arrington et al."; }

    private:
      void compute(double q2) override;

      const int mode_;
      std::vector<double> a_e_, b_e_;
      std::vector<double> a_m_, b_m_;
    };

    ArringtonEtAl::ArringtonEtAl(const ParametersList& params)
        : Parameterisation(params), mode_(params.get<int>("mode")) {
      switch (mode_) {
        case 0:  // original
          a_e_ = {3.439, -1.602, 0.068};
          b_e_ = {15.055, 48.061, 99.304, 0.012, 8.650};
          a_m_ = {-1.465, 1.260, 0.262};
          b_m_ = {9.627, 0., 0., 11.179, 13.245};
          break;
        case 1:  // fit of quoted Ge+dGe values
          a_e_ = {4.309, -1.108, -0.324};
          b_e_ = {15.340, 58.321, 124.11, 3.927, 0.589};
          a_m_ = {-1.472, 1.210, 0.334};
          b_m_ = {9.486, 0., 0., 9.440, 15.416};
          break;
        case 2:  // fit of quoted Ge-dGe values
          a_e_ = {4.286, -1.281, -0.486};
          b_e_ = {16.308, 54.535, 138.03, 7.005, 0.014};
          a_m_ = {-1.374, 1.080, 0.124};
          b_m_ = {10.003, 0., 0., 7.680, 9.009};
          break;
        case 3:  // fit of quoted Ge values
          a_e_ = {4.109, -1.052, -0.375};
          b_e_ = {15.602, 55.519, 123.96, 11.403, 1.931};
          a_m_ = {-1.436, 1.196, 0.210};
          b_m_ = {9.721, 0., 0., 9.623, 11.817};
          break;
      }
    }

    void ArringtonEtAl::compute(double q2) {
      const double tau_val = tau(q2);

      double num_e = 1., den_e = 1.;
      for (size_t i = 0; i < a_e_.size(); ++i)
        num_e += a_e_.at(i) * pow(tau_val, 1. + i);
      for (size_t i = 0; i < b_e_.size(); ++i)
        den_e += b_e_.at(i) * pow(tau_val, 1. + i);
      GE = num_e / den_e;

      double num_m = 1., den_m = 1.;
      for (size_t i = 0; i < a_m_.size(); ++i)
        num_m += a_m_.at(i) * pow(tau_val, 1. + i);
      for (size_t i = 0; i < b_m_.size(); ++i)
        den_m += b_m_.at(i) * pow(tau_val, 1. + i);
      GM = MU * num_m / den_m;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL("Arrington", ArringtonEtAl)
