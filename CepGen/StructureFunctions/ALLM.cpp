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

#include <cassert>
#include <cmath>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Physics.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{2,L}\f$ parameterisation by Abramowicz, Levin, Levy, and Maor \cite Abramowicz:1991xz\cite Abramowicz:1997ms
    class ALLM final : public Parameterisation {
    public:
      class Parameters {
      private:
        struct Trajectory {
          explicit Trajectory(const ParametersList& params = ParametersList());
          std::vector<double> a, b, c;
        };

      public:
        explicit Parameters(const ParametersList&);
        /// Pre-HERA data fit (694 data points)
        static Parameters allm91();
        /// Fixed target and HERA photoproduction total cross sections (1356 points)
        static Parameters allm97();
        static Parameters hht_allm();
        static Parameters hht_allm_ft();
        static Parameters gd07p();
        static Parameters gd11p();

        Trajectory pomeron, reggeon;
        /// Effective photon squared mass
        double m02;
        /// Effective pomeron squared mass
        double mp2;
        /// Effective reggeon squared mass
        double mr2;
        double q02;
        /// Squared QCD scale
        double lambda2;
        Type type;
      };

      explicit ALLM(const ParametersList& params = ParametersList());
      static std::string description() { return "Abramowicz, Levin, Levy, and Maor parametrisation of F2/FL"; }

      ALLM& eval(double xbj, double q2) override;
      std::string describe() const override { return descr_; }

    private:
      Parameters mod_params_;
      std::string descr_;
    };

    ALLM::ALLM(const ParametersList& params)
        : Parameterisation(params), mod_params_(params.get<ParametersList>("parameterisation")) {
      const auto& model = params.get<std::string>("model");
      if (model == "GD07p") {
        mod_params_ = Parameters::gd07p();
        descr_ = "ALLM{GD07p}";
      } else if (model == "GD11p") {
        mod_params_ = Parameters::gd11p();
        descr_ = "ALLM{GD11p}";
      } else if (model == "ALLM91") {
        mod_params_ = Parameters::allm91();
        descr_ = "ALLM{91}";
      } else if (model == "ALLM97") {
        mod_params_ = Parameters::allm97();
        descr_ = "ALLM{97}";
      } else if (model == "HHT_ALLM") {
        mod_params_ = Parameters::hht_allm();
        descr_ = "ALLM{HHT}";
      } else if (model == "HHT_ALLM_FT") {
        mod_params_ = Parameters::hht_allm_ft();
        descr_ = "ALLM{HHT_FT}";
      }
      CG_DEBUG("ALLM") << "ALLM structure functions builder initialised.\n"
                       << "Parameterisation (" << mod_params_.type << "):\n"
                       << " *) Pomeron trajectory:\n"
                       << "   a = {" << mod_params_.pomeron.a.at(0) << ", " << mod_params_.pomeron.a.at(1) << ", "
                       << mod_params_.pomeron.a.at(2) << "}\n"
                       << "   b = {" << mod_params_.pomeron.b.at(0) << ", " << mod_params_.pomeron.b.at(1) << ", "
                       << mod_params_.pomeron.b.at(2) << "}\n"
                       << "   c = {" << mod_params_.pomeron.c.at(0) << ", " << mod_params_.pomeron.c.at(1) << ", "
                       << mod_params_.pomeron.c.at(2) << "}\n"
                       << " *) Reggeon trajectory:\n"
                       << "   a = {" << mod_params_.reggeon.a.at(0) << ", " << mod_params_.reggeon.a.at(1) << ", "
                       << mod_params_.reggeon.a.at(2) << "}\n"
                       << "   b = {" << mod_params_.reggeon.b.at(0) << ", " << mod_params_.reggeon.b.at(1) << ", "
                       << mod_params_.reggeon.b.at(2) << "}\n"
                       << "   c = {" << mod_params_.reggeon.c.at(0) << ", " << mod_params_.reggeon.c.at(1) << ", "
                       << mod_params_.reggeon.c.at(2) << "}\n"
                       << " masses: m_0^2=" << mod_params_.m02 << ", m_p^2=" << mod_params_.mp2
                       << ", m_r^2=" << mod_params_.mr2 << " GeV^2\n"
                       << " q_0^2=" << mod_params_.q02 << ", Lambda^2=" << mod_params_.lambda2 << " GeV^2.";
    }

    ALLM& ALLM::eval(double xbj, double q2) {
      const double w2_eff = utils::mX2(xbj, q2, mp2_) - mp2_;
      const double xp = (q2 + mod_params_.mp2) / (q2 + w2_eff + mod_params_.mp2),
                   xr = (q2 + mod_params_.mr2) / (q2 + w2_eff + mod_params_.mr2);

      const double xlog1 = log((q2 + mod_params_.q02) / mod_params_.lambda2),
                   xlog2 = log(mod_params_.q02 / mod_params_.lambda2);
      const double t = log(xlog1 / xlog2);

      const double apom = mod_params_.pomeron.a.at(0) + (mod_params_.pomeron.a.at(0) - mod_params_.pomeron.a.at(1)) *
                                                            (1. / (1. + pow(t, mod_params_.pomeron.a.at(2))) - 1.);
      const double bpom =
          mod_params_.pomeron.b.at(0) + mod_params_.pomeron.b.at(1) * pow(t, mod_params_.pomeron.b.at(2));
      const double cpom = mod_params_.pomeron.c.at(0) + (mod_params_.pomeron.c.at(0) - mod_params_.pomeron.c.at(1)) *
                                                            (1. / (1. + pow(t, mod_params_.pomeron.c.at(2))) - 1.);

      const double areg =
          mod_params_.reggeon.a.at(0) + mod_params_.reggeon.a.at(1) * pow(t, mod_params_.reggeon.a.at(2));
      const double breg =
          mod_params_.reggeon.b.at(0) + mod_params_.reggeon.b.at(1) * pow(t, mod_params_.reggeon.b.at(2));
      const double creg =
          mod_params_.reggeon.c.at(0) + mod_params_.reggeon.c.at(1) * pow(t, mod_params_.reggeon.c.at(2));

      const double F2_Pom = cpom * pow(xp, apom) * pow(1. - xbj, bpom),
                   F2_Reg = creg * pow(xr, areg) * pow(1. - xbj, breg);

      F2 = q2 / (q2 + mod_params_.m02) * (F2_Pom + F2_Reg);

      return *this;
    }

    //---------------------------------------------------------------------------------------------
    // parameterisation object
    //---------------------------------------------------------------------------------------------

    ALLM::Parameters::Parameters(const ParametersList& params)
        : pomeron(params.get<ParametersList>("pomeronTrajectory")),
          reggeon(params.get<ParametersList>("reggeonTrajectory")),
          m02(params.get<double>("m02")),
          mp2(params.get<double>("mp2")),
          mr2(params.get<double>("mr2")),
          q02(params.get<double>("q02")),
          lambda2(params.get<double>("lambda2")),
          type(params.getAs<int, Type>("type", Type::Invalid)) {}

    ALLM::Parameters ALLM::Parameters::allm91() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.04503, -0.36407, 8.17091})
                                                       .set<std::vector<double> >("b", {0.49222, 0.52116, 3.5515})
                                                       .set<std::vector<double> >("c", {0.26550, 0.04856, 1.04682}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.60408, 0.17353, 1.61812})
                                                       .set<std::vector<double> >("b", {1.26066, 1.83624, 0.81141})
                                                       .set<std::vector<double> >("c", {0.67639, 0.49027, 2.66275}))
                              .set<double>("m02", 0.30508)
                              .set<double>("mp2", 10.676)
                              .set<double>("mr2", 0.20623)
                              .set<double>("q02", 0.27799)
                              .set<double>("lambda2", 0.06527)
                              .set<int>("type", (int)Type::ALLM91));
      return p;
    }

    ALLM::Parameters ALLM::Parameters::allm97() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.0808, -0.44812, 1.1709})
                                                       .set<std::vector<double> >("b", {0.36292, 1.8917, 1.8439})
                                                       .set<std::vector<double> >("c", {0.28067, 0.22291, 2.1979}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.58400, 0.37888, 2.6063})
                                                       .set<std::vector<double> >("b", {0.01147, 3.7582, 0.49338})
                                                       .set<std::vector<double> >("c", {0.80107, 0.97307, 3.4924}))
                              .set<double>("m02", 0.31985)
                              .set<double>("mp2", 49.457)
                              .set<double>("mr2", 0.15052)
                              .set<double>("q02", 0.52544)
                              .set<double>("lambda2", 0.06526)
                              .set<int>("type", (int)Type::ALLM97));
      return p;
    }

    ALLM::Parameters ALLM::Parameters::hht_allm() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.835, -0.446, 10.6})
                                                       .set<std::vector<double> >("b", {-45.8, 55.7, -0.031})
                                                       .set<std::vector<double> >("c", {0.412, 0.164, 17.7}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.706, 0.185, -16.4})
                                                       .set<std::vector<double> >("b", {-1.29, 4.51, 1.16})
                                                       .set<std::vector<double> >("c", {-1.04, 2.97, 0.163}))
                              .set<double>("m02", 0.446)
                              .set<double>("mp2", 74.2)
                              .set<double>("mr2", 29.3)
                              .set<double>("q02", 4.74e-5)
                              .set<double>("lambda2", 2.2e-8));
      return p;
    }

    ALLM::Parameters ALLM::Parameters::hht_allm_ft() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.075, -0.470, 9.2})
                                                       .set<std::vector<double> >("b", {-0.477, 54.0, 0.073})
                                                       .set<std::vector<double> >("c", {0.356, 0.171, 18.6}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.882, 0.082, -8.5})
                                                       .set<std::vector<double> >("b", {0.339, 3.38, 1.07})
                                                       .set<std::vector<double> >("c", {-0.636, 3.37, -0.660}))
                              .set<double>("m02", 0.388)
                              .set<double>("mp2", 50.8)
                              .set<double>("mr2", 0.838)
                              .set<double>("q02", 1.87e-5)
                              .set<double>("lambda2", 4.4e-9));
      return p;
    }

    ALLM::Parameters ALLM::Parameters::gd07p() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.105, -0.495, 1.29})
                                                       .set<std::vector<double> >("b", {-1.42, 4.51, 0.551})
                                                       .set<std::vector<double> >("c", {0.339, 0.127, 1.16}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.374, 0.998, 0.775})
                                                       .set<std::vector<double> >("b", {2.71, 1.83, 1.26})
                                                       .set<std::vector<double> >("c", {0.838, 2.36, 1.77}))
                              .set<double>("m02", 0.454)
                              .set<double>("mp2", 30.7)
                              .set<double>("mr2", 0.117)
                              .set<double>("q02", 1.15)
                              .set<double>("lambda2", 0.06527)
                              .set<int>("type", (int)Type::GD07p));
      return p;
    }

    ALLM::Parameters ALLM::Parameters::gd11p() {
      static Parameters p(ParametersList()
                              .set<ParametersList>("pomeronTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {-0.11895, -0.4783, 1.353})
                                                       .set<std::vector<double> >("b", {1.0833, 2.656, 1.771})
                                                       .set<std::vector<double> >("c", {0.3638, 0.1211, 1.166}))
                              .set<ParametersList>("reggeonTrajectory",
                                                   ParametersList()
                                                       .set<std::vector<double> >("a", {0.3425, 1.0603, 0.5164})
                                                       .set<std::vector<double> >("b", {-10.408, 14.857, 0.07739})
                                                       .set<std::vector<double> >("c", {1.3633, 2.256, 2.209}))
                              .set<double>("m02", 0.5063)
                              .set<double>("mp2", 34.75)
                              .set<double>("mr2", 0.03190)
                              .set<double>("q02", 1.374)
                              .set<double>("lambda2", 0.06527)
                              .set<int>("type", (int)Type::GD11p));
      return p;
    }

    ALLM::Parameters::Trajectory::Trajectory(const ParametersList& params)
        : a(params.get<std::vector<double> >("a", {0., 0., 0.})),
          b(params.get<std::vector<double> >("b", {0., 0., 0.})),
          c(params.get<std::vector<double> >("c", {0., 0., 0.})) {
      assert(a.size() == 3);
      assert(b.size() == 3);
      assert(c.size() == 3);
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN_PARAMS(ALLM91, strfun::ALLM, ParametersList().set<std::string>("model", "ALLM91"))
REGISTER_STRFUN_PARAMS(ALLM97, strfun::ALLM, ParametersList().set<std::string>("model", "ALLM97"))
REGISTER_STRFUN_PARAMS(GD07p, strfun::ALLM, ParametersList().set<std::string>("model", "GD07p"))
REGISTER_STRFUN_PARAMS(GD11p, strfun::ALLM, ParametersList().set<std::string>("model", "GD11p"))
