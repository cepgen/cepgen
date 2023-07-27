/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule<int>(params),
          r_ratio_(SigmaRatiosFactory::get().build(steer<int>("sigmaRatio"))),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          mx_min_(mp_ + PDG::get().mass(PDG::piZero)) {
      CG_DEBUG("Parameterisation") << "Structure functions parameterisation to be built using following parameters:\n"
                                   << ParametersDescription(params_).describe(true);
    }

    Parameterisation& Parameterisation::operator()(double xbj, double q2) {
      const auto args = Arguments{xbj, q2};
      if (args == args_)
        return *this;
      clear();
      if (!args.valid()) {
        CG_WARNING("StructureFunctions") << "Invalid range for Q² = " << q2 << " or xBj = " << xbj << ".";
        return *this;
      }
      args_ = args;
      eval();
      return *this;
    }

    Parameterisation& Parameterisation::clear() {
      f2_ = fl_ = w1_ = w2_ = fe_ = fm_ = 0.;
      fl_computed_ = false;
      return *this;
    }

    double Parameterisation::F2(double xbj, double q2) { return operator()(xbj, q2).f2_; }

    double Parameterisation::FL(double xbj, double q2) {
      if (!fl_computed_)
        computeFL(xbj, q2);
      return operator()(xbj, q2).fl_;
    }

    double Parameterisation::W1(double xbj, double q2) { return operator()(xbj, q2).w1_; }

    double Parameterisation::W2(double xbj, double q2) { return operator()(xbj, q2).w2_; }

    double Parameterisation::FE(double xbj, double q2) { return operator()(xbj, q2).fe_; }

    double Parameterisation::FM(double xbj, double q2) { return operator()(xbj, q2).fm_; }

    double Parameterisation::F1(double xbj, double q2) {
      return 0.5 * (gamma2(xbj, q2) * F2(xbj, q2) - FL(xbj, q2)) / xbj;
    }

    Parameterisation& Parameterisation::setF1F2(double f1, double f2) {
      return (*this)
          .setF2(f2)  // trivial
          .setFL(gamma2(args_.xbj, args_.q2) * f2_ - 2. * f1 * args_.xbj);
    }

    Parameterisation& Parameterisation::setF2(double f2) {
      f2_ = f2;
      return *this;
    }

    Parameterisation& Parameterisation::setFL(double fl) {
      fl_ = fl;
      fl_computed_ = true;
      return *this;
    }

    Parameterisation& Parameterisation::setW1(double w1) {
      w1_ = w1;
      return *this;
    }

    Parameterisation& Parameterisation::setW2(double w2) {
      w2_ = w2;
      return *this;
    }

    Parameterisation& Parameterisation::setFE(double fe) {
      fe_ = fe;
      return *this;
    }

    Parameterisation& Parameterisation::setFM(double fm) {
      fm_ = fm;
      return *this;
    }

    double Parameterisation::tau(double xbj, double q2) const { return 4. * xbj * xbj * mp2_ / q2; }

    double Parameterisation::gamma2(double xbj, double q2) const { return 1. + tau(xbj, q2); }

    Parameterisation& Parameterisation::computeFL(double xbj, double q2) {
      if (!fl_computed_ && !r_ratio_)
        throw CG_FATAL("StructureFunctions:FL") << "Failed to retrieve a R-ratio calculator!";
      double r_error;
      return computeFL(xbj, q2, (*r_ratio_)(xbj, q2, r_error));
    }

    Parameterisation& Parameterisation::computeFL(double xbj, double q2, double r) {
      if (!fl_computed_)
        setFL(f2_ * gamma2(xbj, q2) * (r / (1. + r)));
      return *this;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& sf) {
      os << sf.description().description();
      if (sf.args_.valid())
        os << " at " << sf.args_ << ": F2 = " << sf.f2_ << ", FL = " << sf.fl_;
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation::Arguments& args) {
      return os << "(" << args.xbj << ", " << args.q2 << ")";
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed structure functions parameterisation");
      desc.add<int>("sigmaRatio", 4 /* SibirtsevBlunden */)
          .setDescription("Modelling for the sigma(L/T) ratio used in FL computation from F2");
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen
