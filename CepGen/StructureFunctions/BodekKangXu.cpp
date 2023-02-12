/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{1,2}\f$ modelling by Bodek, Kang, and Xu \cite Bodek:2021bde
    class BodekKangXu : public Parameterisation {
    public:
      explicit BodekKangXu(const ParametersList& params)
          : Parameterisation(params),
            constants_(steer<std::vector<double> >("constants")),
            pi_em_sq_(constants_.empty() ? 0. : std::pow(constants_.at(0) - mp_, 2.)),
            spins_(steer<std::vector<int> >("spins")),
            r_(steer<double>("r")) {
        if (constants_.size() != 24)
          throw CG_FATAL("BodekKangXu") << "Invalid parameters multiplicity given. Should have size 24, has size "
                                        << constants_.size() << ".";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Bodek, Kang, and Xu");
        desc.add<std::vector<double> >(
            "constants",
            {1.0741163,  0.75531124,  3.3506491,    1.7447015,    3.5102405, 1.040004,    1.2299128, 0.10625394,
             0.48132786, 1.5101467,   0.081661975,  0.65587179,   1.7176216, 0.12551987,  0.7473379, 1.953819,
             0.19891522, -0.17498537, 0.0096701919, -0.035256748, 3.5185207, -0.59993696, 4.7615828, 0.41167589});
        desc.add<std::vector<int> >("spins", {1, 2, 3, 2});
        desc.add<double>("r", 0.18);
        desc.add<double>("q0", 1.);
        return desc;
      }

      BodekKangXu& eval(double xbj, double q2) override {
        const auto mx2 = utils::mX2(xbj, q2, mp2_);
        if (mx2 < mp2_) {
          setF1F2(0., 0.);
          return *this;
        }
        const auto q0 = 0.5 * q2 / mp_ / xbj;
        const auto w2h = gp_h(q0, q2) * bodek(std::sqrt(mx2), q2) / q0, w1h = (1. + q0 * q0 / q2) / (1. + r_) * w2h;
        setF1F2(mp_ * w1h, q0 * w2h);
        return *this;
      }

    private:
      float gp_h(float q0, float q2) const {
        const auto gi = 2. * mp_ * q0;
        const auto ww = (gi + 1.642) / (q2 + 0.376), t = 1. - 1. / ww;
        const auto wp = 0.256 * std::pow(t, 3.) + 2.178 * std::pow(t, 4.) + 0.898 * std::pow(t, 5.) -
                        6.716 * std::pow(t, 6.) + 3.756 * std::pow(t, 7.);
        return wp * ww * q2 / gi;
      }
      float bodek(double w, double q2) const {
        const size_t NRES = 4, NBKG = 5;

        if (w <= mp_)
          return 0.;
        const auto w2 = w * w;
        const float omega = 1. + (w2 - mp2_) / q2, x = 1. / omega;
        const float xpx = constants_.at(21) + constants_.at(22) * std::pow(x - constants_.at(23), 2.);

        auto b1 = 0., b2 = 0.;
        if (w != constants_.at(0))
          b1 = std::max(0., w - constants_.at(0)) / (w - constants_.at(0)) * constants_.at(1);
        const auto eb1 = constants_.at(2) * (w - constants_.at(0));

        if (eb1 <= 25.) {
          b1 *= (1. - exp(-eb1));
          b2 = 0.;
        }
        if (w != constants_.at(3))
          b2 = std::max(0., w - constants_.at(3)) / (w - constants_.at(3)) * (1. - constants_.at(1));

        const auto eb2 = constants_.at(4) * (w2 - std::pow(constants_.at(3), 2.));

        if (eb2 <= 25.)
          b2 *= (1. - exp(-eb2));

        const auto BBKG = b1 + b2, BRES = constants_.at(1) + b2;

        auto res = 0., ressum = 0.;
        for (size_t i = 0; i < NRES; ++i) {
          const size_t index = i * 3 + 1 + NBKG;
          auto ram = constants_.at(index), rma = constants_.at(index + 1), rwd = constants_.at(index + 2);
          if (i == 0)
            ram += constants_.at(17) * q2 + constants_.at(18) * q2 * q2;
          if (i == 2)
            rma *= (1. + constants_.at(19) / (1. + constants_.at(20) * q2));
          const auto qstarn = std::sqrt(std::max(0., std::pow((w2 + mp2_ - pi_em_sq_) / (2. * w), 2.) - mp2_));
          const auto qstar0 =
              std::sqrt(std::max(0., std::pow((rma * rma - mp2_ + pi_em_sq_) / (2. * rma), 2.) - pi_em_sq_));
          if (qstar0 <= 1.e-10) {
            res = 0.;
            ressum += res;
            continue;
          }

          const auto term = prefactor_ * qstarn, term0 = prefactor_ * qstar0;
          const size_t j = 2 * spins_.at(i);
          const auto gamres =
              0.5 * (rwd * std::pow(term / term0, j + 1) * (1. + std::pow(term0, j)) / (1. + std::pow(term, j)));
          const auto brwig = M_1_PI * gamres / (std::pow(w - rma, 2.) + std::pow(gamres, 2.));
          res = ram * brwig / 2. / mp_;
          ressum += res;
        }

        return BBKG * (1. + (1. - BBKG) * xpx) + ressum * (1. - BRES);
      }
      const std::vector<double> constants_;
      const double pi_em_sq_;
      const std::vector<int> spins_;
      const double r_;
      static constexpr double prefactor_ = 6.08974;
    };
  }  // namespace strfun
}  // namespace cepgen
typedef cepgen::strfun::BodekKangXu BodekKangXu;
REGISTER_STRFUN(304, BodekKangXu);
