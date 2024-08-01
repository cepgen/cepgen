/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// \f$F_{2,L}\f$ parameterisation by Abramowicz, Levin, Levy, and Maor \cite Abramowicz:1991xz\cite Abramowicz:1997ms
  class ALLM : public Parameterisation {
  public:
    explicit ALLM(const ParametersList&);

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Abramowicz, Levin, Levy, and Maor (continuum)");
      desc.add<ParametersDescription>("pomeronTrajectory", Trajectory::description());
      desc.add<ParametersDescription>("reggeonTrajectory", Trajectory::description());
      desc.add<double>("m02", 0.).setDescription("effective photon squared mass (in GeV^2)");
      desc.add<double>("mp2", 0.).setDescription("effective pomeron squared mass (in GeV^2)");
      desc.add<double>("mr2", 0.).setDescription("effective reggeon squared mass (in GeV^2)");
      desc.add<double>("q02", 0.);
      desc.add<double>("lambda2", 0.).setDescription("squared QCD scale (in GeV^2");
      return desc;
    }

    void eval() override;

  private:
    class Trajectory : public SteeredObject<Trajectory> {
    public:
      explicit Trajectory(const ParametersList& params)
          : SteeredObject(params),
            a_(steer<std::vector<double> >("a")),
            b_(steer<std::vector<double> >("b")),
            c_(steer<std::vector<double> >("c")) {
        assert(a_.size() == 3);
        assert(b_.size() == 3);
        assert(c_.size() == 3);
      }

      static ParametersDescription description() {
        auto desc = ParametersDescription();
        desc.add<std::vector<double> >("a", {0., 0., 0.});
        desc.add<std::vector<double> >("b", {0., 0., 0.});
        desc.add<std::vector<double> >("c", {0., 0., 0.});
        return desc;
      }

      friend std::ostream& operator<<(std::ostream& os, const Trajectory& tr) {
        auto print_vec = [](const auto& vec) -> std::string {
          std::ostringstream os;
          std::string sep;
          for (const auto& it : vec)
            os << sep << it, sep = ", ";
          return os.str();
        };
        return os << "[a = " << print_vec(tr.a_) << ", b = " << print_vec(tr.b_) << ", c = " << print_vec(tr.c_) << "]";
      }

      double eval1(char mod, double t) const {
        const auto& par = get(mod);
        return par.at(0) + (par.at(0) - par.at(1)) * (1. / (1. + pow(t, par.at(2))) - 1.);
      }
      double eval2(char mod, double t) const {
        const auto& par = get(mod);
        return par.at(0) + par.at(1) * pow(t, par.at(2));
      }

    private:
      const std::vector<double>& get(char mod) const {
        switch (mod) {
          case 'a':
            return a_;
          case 'b':
            return b_;
          case 'c':
            return c_;
          default:
            throw CG_FATAL("ALLM:Trajectory") << "Failed to retrieve parameterisation for char{'" << mod << "'}!";
        }
      }
      std::vector<double> a_, b_, c_;
    };
    Trajectory pomeron_, reggeon_;
    /// Effective photon squared mass
    double m02_;
    /// Effective pomeron squared mass
    double mpom2_;
    /// Effective reggeon squared mass
    double mreg2_;
    double q02_;
    /// Squared QCD scale
    double lambda2_;
  };

  ALLM::ALLM(const ParametersList& params)
      : Parameterisation(params),
        pomeron_(steer<ParametersList>("pomeronTrajectory")),
        reggeon_(steer<ParametersList>("reggeonTrajectory")),
        m02_(steer<double>("m02")),
        mpom2_(steer<double>("mp2")),
        mreg2_(steer<double>("mr2")),
        q02_(steer<double>("q02")),
        lambda2_(steer<double>("lambda2")) {
    CG_DEBUG("ALLM") << "ALLM structure functions builder initialised.\n"
                     << " *) Pomeron trajectory: " << pomeron_ << "\n"
                     << " *) Reggeon trajectory: " << reggeon_ << "\n"
                     << " masses: m_0^2=" << m02_ << ", m_p^2=" << mpom2_ << ", m_r^2=" << mreg2_ << " GeV^2\n"
                     << " q_0^2=" << q02_ << ", Lambda^2=" << lambda2_ << " GeV^2.";
  }

  void ALLM::eval() {
    const double w2_eff = utils::mX2(args_.xbj, args_.q2, mp2_) - mp2_;
    const double xp = (args_.q2 + mpom2_) / (args_.q2 + w2_eff + mpom2_),
                 xr = (args_.q2 + mreg2_) / (args_.q2 + w2_eff + mreg2_);

    const double xlog1 = log((args_.q2 + q02_) / lambda2_), xlog2 = log(q02_ / lambda2_);
    const double t = log(xlog1 / xlog2);

    const double apom = pomeron_.eval1('a', t), bpom = pomeron_.eval2('b', t), cpom = pomeron_.eval1('c', t);
    const double areg = reggeon_.eval2('a', t), breg = reggeon_.eval2('b', t), creg = reggeon_.eval2('c', t);

    const double F2_Pom = cpom * pow(xp, apom) * pow(1. - args_.xbj, bpom),
                 F2_Reg = creg * pow(xr, areg) * pow(1. - args_.xbj, breg);

    setF2(args_.q2 / (args_.q2 + m02_) * (F2_Pom + F2_Reg));
  }

  //---------------------------------------------------------------------------------------------
  // ALLM parameterisations
  //---------------------------------------------------------------------------------------------

  /// Pre-HERA data fit (694 data points)
  struct ALLM91 final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("ALLM91 (continuum, pre-HERA 694 points fit)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.04503, -0.36407, 8.17091});
      pom_params.add<std::vector<double> >("b", {0.49222, 0.52116, 3.5515});
      pom_params.add<std::vector<double> >("c", {0.26550, 0.04856, 1.04682});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.60408, 0.17353, 1.61812});
      reg_params.add<std::vector<double> >("b", {1.26066, 1.83624, 0.81141});
      reg_params.add<std::vector<double> >("c", {0.67639, 0.49027, 2.66275});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.30508);
      desc.add<double>("mp2", 10.676);
      desc.add<double>("mr2", 0.20623);
      desc.add<double>("q02", 0.27799);
      desc.add<double>("lambda2", 0.06527);
      return desc;
    }
  };

  /// Fixed target and HERA photoproduction total cross sections (1356 points)
  struct ALLM97 final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("ALLM97 (continuum, FT/HERA photoprod. tot.x-s 1356 points fit)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.0808, -0.44812, 1.1709});
      pom_params.add<std::vector<double> >("b", {0.36292, 1.8917, 1.8439});
      pom_params.add<std::vector<double> >("c", {0.28067, 0.22291, 2.1979});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.58400, 0.37888, 2.6063});
      reg_params.add<std::vector<double> >("b", {0.01147, 3.7582, 0.49338});
      reg_params.add<std::vector<double> >("c", {0.80107, 0.97307, 3.4924});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.31985);
      desc.add<double>("mp2", 49.457);
      desc.add<double>("mr2", 0.15052);
      desc.add<double>("q02", 0.52544);
      desc.add<double>("lambda2", 0.06526);
      return desc;
    }
  };

  struct HHTALLM final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("HHT (continuum)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.835, -0.446, 10.6});
      pom_params.add<std::vector<double> >("b", {-45.8, 55.7, -0.031});
      pom_params.add<std::vector<double> >("c", {0.412, 0.164, 17.7});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.706, 0.185, -16.4});
      reg_params.add<std::vector<double> >("b", {-1.29, 4.51, 1.16});
      reg_params.add<std::vector<double> >("c", {-1.04, 2.97, 0.163});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.446);
      desc.add<double>("mp2", 74.2);
      desc.add<double>("mr2", 29.3);
      desc.add<double>("q02", 4.74e-5);
      desc.add<double>("lambda2", 2.2e-8);
      return desc;
    }
  };

  struct HHTALLMFT final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("HHT_FT (continuum)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.075, -0.470, 9.2});
      pom_params.add<std::vector<double> >("b", {-0.477, 54.0, 0.073});
      pom_params.add<std::vector<double> >("c", {0.356, 0.171, 18.6});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.882, 0.082, -8.5});
      reg_params.add<std::vector<double> >("b", {0.339, 3.38, 1.07});
      reg_params.add<std::vector<double> >("c", {-0.636, 3.37, -0.660});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.388);
      desc.add<double>("mp2", 50.8);
      desc.add<double>("mr2", 0.838);
      desc.add<double>("q02", 1.87e-5);
      desc.add<double>("lambda2", 4.4e-9);
      return desc;
    }
  };

  struct GD07p final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("GD07p (continuum)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.105, -0.495, 1.29});
      pom_params.add<std::vector<double> >("b", {-1.42, 4.51, 0.551});
      pom_params.add<std::vector<double> >("c", {0.339, 0.127, 1.16});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.374, 0.998, 0.775});
      reg_params.add<std::vector<double> >("b", {2.71, 1.83, 1.26});
      reg_params.add<std::vector<double> >("c", {0.838, 2.36, 1.77});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.454);
      desc.add<double>("mp2", 30.7);
      desc.add<double>("mr2", 0.117);
      desc.add<double>("q02", 1.15);
      desc.add<double>("lambda2", 0.06527);
      return desc;
    }
  };

  struct GD11p final : public ALLM {
    using ALLM::ALLM;
    static ParametersDescription description() {
      auto desc = ALLM::description();
      desc.setDescription("GD11p (continuum)");
      ParametersDescription pom_params;
      pom_params.add<std::vector<double> >("a", {-0.11895, -0.4783, 1.353});
      pom_params.add<std::vector<double> >("b", {1.0833, 2.656, 1.771});
      pom_params.add<std::vector<double> >("c", {0.3638, 0.1211, 1.166});
      desc.add("pomeronTrajectory", pom_params);
      ParametersDescription reg_params;
      reg_params.add<std::vector<double> >("a", {0.3425, 1.0603, 0.5164});
      reg_params.add<std::vector<double> >("b", {-10.408, 14.857, 0.07739});
      reg_params.add<std::vector<double> >("c", {1.3633, 2.256, 2.209});
      desc.add("reggeonTrajectory", reg_params);
      desc.add<double>("m02", 0.5063);
      desc.add<double>("mp2", 34.75);
      desc.add<double>("mr2", 0.03190);
      desc.add<double>("q02", 1.374);
      desc.add<double>("lambda2", 0.06527);
      return desc;
    }
  };
}  // namespace cepgen::strfun
using cepgen::strfun::ALLM91;
using cepgen::strfun::ALLM97;
using cepgen::strfun::GD07p;
using cepgen::strfun::GD11p;
using cepgen::strfun::HHTALLM;
using cepgen::strfun::HHTALLMFT;
REGISTER_STRFUN("ALLM91", 201, ALLM91);
REGISTER_STRFUN("ALLM97", 202, ALLM97);
REGISTER_STRFUN("GD07p", 203, GD07p);
REGISTER_STRFUN("GD11p", 204, GD11p);
REGISTER_STRFUN("htALLM", 206, HHTALLM);
REGISTER_STRFUN("htALLMft", 207, HHTALLMFT);
