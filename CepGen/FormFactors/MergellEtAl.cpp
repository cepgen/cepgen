#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include <cassert>
#include <cmath>

namespace cepgen {
  namespace formfac {
    /// \cite Mergell:1995bf
    class MergellEtAl : public Parameterisation {
    public:
      MergellEtAl(const ParametersList&);
      static std::string description() { return "Mergell et al."; }

    private:
      void compute(double q2) override;
      static constexpr double Q2_RESCL = 9.733, INV_DENUM = 1. / 0.350;
      static constexpr double EXPO = 2.148;
      const std::vector<double> par1_, par2_;
    };

    MergellEtAl::MergellEtAl(const ParametersList& params)
        : Parameterisation(params),
          par1_(params.get<std::vector<double> >("par1", {1.0317, 0.0875, 0.3176, 0.5496})),
          par2_(params.get<std::vector<double> >("par2", {5.7824, 0.3907, 0.1422, 0.5362})) {
      assert(par1_.size() == 4);
      assert(par2_.size() == 4);
    }

    void MergellEtAl::compute(double q2) {
      const double log1 = std::pow(log((Q2_RESCL + q2) * INV_DENUM), -EXPO);
      const double d1_1 = 0.611 + q2, d2_1 = 1.039 + q2, d3_1 = 2.560 + q2;

      const double Fs1 = (9.464 / d1_1 - 9.054 / d2_1 - 0.410 / d3_1) * log1;
      const double Fs2 = (-1.549 / d1_1 + 1.985 / d2_1 - 0.436 / d3_1) * log1;

      const double log2 = std::pow(log((Q2_RESCL - 0.500) * INV_DENUM), +EXPO);
      const double log3 = std::pow(log((Q2_RESCL - 0.400) * INV_DENUM), +EXPO);
      const double d1_2 = 2.103 + q2, d2_2 = 2.734 + q2, d3_2 = 2.835 + q2;

      const double Fv1 = (0.5 * (par1_.at(0) * log2 + par1_.at(1) * log3 * std::pow(1. + q2 / par1_.at(2), -2)) /
                              (1. + q2 / par1_.at(3)) -
                          38.885 / d1_2 + 425.007 / d2_2 - 389.742 / d3_2) *
                         log1;
      const double Fv2 =
          (0.5 * (par2_.at(0) * log2 + par2_.at(1) * log3 / (1. + q2 / par2_.at(2))) / (1. + q2 / par2_.at(3)) -
           73.535 / d1_2 + 83.211 / d2_2 - 29.467 / d3_2) *
          log1;

      const double F1 = Fv1 + Fs1, F2 = Fv2 + Fs2;

      GE = F1 - tau(q2) * F2;
      GM = F1 + F2;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL("Mergell", MergellEtAl)
