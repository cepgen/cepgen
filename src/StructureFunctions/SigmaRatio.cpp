/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;
using namespace cepgen::sigrat;

Parameterisation::Parameterisation(const ParametersList& params)
    : NamedModule(params), mp_(PDG::get().mass(PDG::proton)), mp2_(mp_ * mp_) {}

double Parameterisation::theta(double xbj, double q2) {
  return 1. + 12. * (q2 / (q2 + 1.)) * (0.125 * 0.125 / (0.125 * 0.125 + xbj * xbj));
}

ParametersDescription Parameterisation::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed sigma ratio parameterisation");
  return desc;
}

//---------------------------------------------------------------------------------------------

namespace cepgen::sigrat {
  /// E143 experimental R measurement \cite Abe:1998ym
  class E143 final : public Parameterisation {
  public:
    explicit E143(const ParametersList& params)
        : Parameterisation(params),
          q2_b_(steer<double>("q2_b")),
          lambda2_(steer<double>("lambda2")),
          a_(steer<std::vector<double> >("a")),
          b_(steer<std::vector<double> >("b")),
          c_(steer<std::vector<double> >("c")) {
      if (a_.size() != 6)
        throw CG_FATAL("E143") << "Parameter 'a' should have 6 components! Parsed " << a_ << ".";
      if (b_.size() != 6)
        throw CG_FATAL("E143") << "Parameter 'b' should have 6 components! Parsed " << b_ << ".";
      if (c_.size() != 6)
        throw CG_FATAL("E143") << "Parameter 'c' should have 6 components! Parsed " << c_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("E143 (experimental)");
      desc.add("q2_b", 0.34);
      desc.add("lambda2", 0.2 * 0.2);
      desc.add("a", std::vector{0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285});
      desc.add("b", std::vector{0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317});
      desc.add("c", std::vector{0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415});
      return desc;
    }

    double operator()(double xbj, double q2, double& err) const override {
      const double u = q2 / q2_b_;
      const double inv_xl = 1. / std::log(q2 / lambda2_);
      const double pa = (1. + a_.at(3) * xbj + a_.at(4) * xbj * xbj) * std::pow(xbj, a_.at(5));
      const double pb = (1. + b_.at(3) * xbj + b_.at(4) * xbj * xbj) * std::pow(xbj, b_.at(5));
      const double q2_thr = c_.at(3) * xbj + c_.at(4) * xbj * xbj + c_.at(5) * xbj * xbj * xbj;
      const double th = theta(xbj, q2);
      // here come the three fits
      const double ra = a_.at(0) * inv_xl * th + a_.at(1) / std::pow(pow(q2, 4) + std::pow(a_.at(2), 4), 0.25) * pa,
                   rb = b_.at(0) * inv_xl * th + (b_.at(1) / q2 + b_.at(2) / (q2 * q2 + 0.3 * 0.3)) * pb,
                   rc = c_.at(0) * inv_xl * th + c_.at(1) / utils::fastHypot(q2 - q2_thr, c_.at(2));

      const double r = (ra + rb + rc) / 3.;  // R is the average of the three fits
      // numerical safety for low-Q²
      err = 0.0078 - 0.013 * xbj + (0.070 - 0.39 * xbj + 0.70 * xbj * xbj) / (1.7 + q2);
      if (q2 > q2_b_)
        return r;
      return r * 0.5 * (3. * u - u * u * u);
    }

  private:
    double q2_b_, lambda2_;
    std::vector<double> a_, b_, c_;
  };

  //---------------------------------------------------------------------------------------------

  /// SLAC experimental R measurement \cite Whitlow:1990gk
  /// \warning valid for \f$Q^2\f$ > 0.3 GeV\f$^2\f$
  class R1990 final : public Parameterisation {
  public:
    explicit R1990(const ParametersList& params)
        : Parameterisation(params), lambda2_(steer<double>("lambda2")), b_(steer<std::vector<double> >("b")) {
      if (b_.size() != 3)
        throw CG_FATAL("R1990") << "Parameter 'b' should have 3 components! Parsed " << b_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("SLAC (experimental)");
      desc.add("lambda2", 0.04);
      desc.add("b", std::vector{0.0635, 0.5747, -0.3534});
      return desc;
    }

    double operator()(double xbj, double q2, double& err) const override {
      err = 0.;
      return b_.at(0) + theta(xbj, q2) / std::log(q2 / lambda2_) + b_.at(1) / q2 + b_.at(2) / (q2 * q2 + 0.09);
    }

  private:
    double lambda2_;
    std::vector<double> b_;
  };

  //---------------------------------------------------------------------------------------------

  /// CLAS experimental R measurement
  class CLAS final : public Parameterisation {
  public:
    explicit CLAS(const ParametersList& params)
        : Parameterisation(params),
          p_(steer<std::vector<double> >("p")),
          wth_(steer<double>("wth")),
          q20_(steer<double>("q20")) {
      if (p_.size() != 3)
        throw CG_FATAL("CLAS") << "Parameter 'p' should have 3 components! Parsed " << p_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("CLAS (experimental)");
      desc.add("p", std::vector{0.041, 0.592, 0.331});
      desc.add("wth", 2.5);
      desc.add("q20", 0.3);
      return desc;
    }

    double operator()(double xbj, double q2, double& err) const override {
      err = 0.;
      //--- 2 kinematic regions: resonances ( w < wth ), and DIS ( w > wth )
      const double w2 = utils::mX2(xbj, q2, mp2_), w = sqrt(w2);
      const double xth = q2 / (q2 + wth_ * wth_ - mp2_);  // xth = x( W = wth )
      const double zeta = std::log(25. * q2);
      const double xi_tmp = (w < wth_) ? theta(xth, q2) : theta(xbj, q2);
      const double tmp = p_.at(0) * xi_tmp / zeta + p_.at(1) / q2 - p_.at(2) / (q20_ * q20_ + q2 * q2);
      if (w >= wth_)
        return tmp;
      return tmp * std::pow((1. - xbj) / (1. - xth), 3);
    }

  private:
    std::vector<double> p_;
    double wth_, q20_;
  };

  //---------------------------------------------------------------------------------------------

  /// Sibirtsev & Blunden parameterisation of the R ratio \cite Sibirtsev:2013cga
  class SibirtsevBlunden final : public Parameterisation {
  public:
    explicit SibirtsevBlunden(const ParametersList& params)
        : Parameterisation(params),
          a_(steer<double>("a")),
          b1_(steer<double>("b1")),
          b2_(steer<double>("b2")),
          c_(steer<double>("c")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Sibirtsev-Blunden (theoretical)");
      desc.add("a", 0.014);
      desc.add("b1", -0.07);
      desc.add("b2", -0.8);
      desc.add("c", 41.);
      return desc;
    }

    double operator()(double /* xbj */, double q2, double& err) const override {
      err = 0.;
      return a_ * q2 * (exp(b1_ * q2) + c_ * exp(b2_ * q2));  // equation (10) of reference paper
    }

  private:
    double a_, b1_, b2_, c_;
  };
}  // namespace cepgen::sigrat
using sigrat::CLAS;
using sigrat::E143;
using sigrat::R1990;
using sigrat::SibirtsevBlunden;
REGISTER_SIGMA_RATIO("E143", 1, E143);
REGISTER_SIGMA_RATIO("R1990", 2, R1990);
REGISTER_SIGMA_RATIO("CLAS", 3, CLAS);
REGISTER_SIGMA_RATIO("SibirtsevBlunden", 4, SibirtsevBlunden);
