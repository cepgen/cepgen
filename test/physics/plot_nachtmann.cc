/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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
#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Physics/NachtmannAmplitudes.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"

using namespace std::complex_literals;

int main(int argc, char* argv[]) {
  std::string integr_name, plotter;
  bool logy, draw_grid;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("integrator,i", "type of integrator used", &integr_name, "Vegas")
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();
  auto integrator = cepgen::IntegratorFactory::get().build(integr_name);

  const std::vector<short> lambdas = {-1, 0, 1};
  auto dsig_dcosth = [&integrator, &lambdas](
                         const cepgen::NachtmannAmplitudes& ampl, double shat, double costh) -> double {
    const auto kin = cepgen::NachtmannAmplitudes::Kinematics::fromScosTheta(shat, costh, 80.379);
    auto dsig = [&ampl, &kin, &lambdas](double theta1, double phi1, double theta2, double phi2) -> double {
      auto l = [](short lam, double theta, double phi) -> std::complex<double> {
        auto dp = [](double theta) -> double { return 0.5 * (1. + cos(theta)) * M_SQRT2; };
        auto d0 = [](double theta) -> double { return sin(theta); };
        auto dm = [](double theta) -> double { return 0.5 * (1. - cos(theta)) * M_SQRT2; };
        if (lam < 0)
          return dp(theta) * exp(std::complex<double>(0., -phi));
        if (lam > 0)
          return dm(theta) * exp(std::complex<double>(0., +phi));
        return -d0(theta);
      };
      auto lbar = [&l](short lam, double theta, double phi) -> std::complex<double> {
        return std::conj(l(lam, theta, phi));
      };

      auto p = [&kin, &ampl](short lam3, short lam4, short lam3p, short lam4p) -> std::complex<double> {
        std::complex<double> p{0.};
        for (const auto& lam1 : {-1, 1})
          for (const auto& lam2 : {-1, 1}) {
            CG_DEBUG("dsig_dcosth") << kin << "?" << ampl(kin, lam1, lam2, lam3, lam4) << "?"
                                    << std::conj(ampl(kin, lam1, lam2, lam3p, lam4p));
            p += ampl(kin, lam1, lam2, lam3, lam4) * std::conj(ampl(kin, lam1, lam2, lam3p, lam4p));
          }
        return p;
      };

      auto d = [&l](short lam, short lamp, double theta, double phi) -> std::complex<double> {
        return l(lam, theta, phi) * std::conj(l(lamp, theta, phi));
      };
      auto dbar = [&lbar](short lam, short lamp, double theta, double phi) -> std::complex<double> {
        return lbar(lam, theta, phi) * std::conj(lbar(lamp, theta, phi));
      };

      double dsig{0.};
      // indices contraction
      for (const auto& lam3p : lambdas)
        for (const auto& lam4p : lambdas)
          for (const auto& lam3 : lambdas)
            for (const auto& lam4 : lambdas) {
              CG_DEBUG("dsig_dcosth") << p(lam3, lam4, lam3p, lam4p) << ":" << d(lam3, lam3p, theta1, phi1) << ":"
                                      << dbar(lam4, lam4p, theta2, phi2);
              dsig += std::norm(p(lam3, lam4, lam3p, lam4p) * d(lam3, lam3p, theta1, phi1) *
                                dbar(lam4, lam4p, theta2, phi2));
            }
      CG_DEBUG("dsig_dcosth") << theta1 << ", " << phi1 << ", " << theta2 << ", " << phi2 << ", " << dsig;
      return dsig;
    };

    cepgen::FunctionIntegrand integr(4, [&dsig](const std::vector<double>& params) {
      double theta1 = acos(params[0] * 2. - 1.), theta2 = acos(params[2] * 2. - 1.);
      double phi1 = params[1] * 2. * M_PI, phi2 = params[3] * 2. * M_PI;
      return dsig(theta1, phi1, theta2, phi2);
    });
    double val, unc;
    integrator->setIntegrand(integr);
    integrator->integrate(val, unc);
    return 3. * kin.beta * std::pow(2., -13) * std::pow(M_1_PI, -3) / shat * val;
  };

  for (auto i = (int)cepgen::NachtmannAmplitudes::Mode::SM; i <= (int)cepgen::NachtmannAmplitudes::Mode::WbarB; ++i) {
    // first get mode name
    if (i > 0)
      break;
    std::ostringstream os;
    os << (cepgen::NachtmannAmplitudes::Mode)i;
    const auto name = os.str();
    CG_LOG << "Computing " << name << ".";
    // then fill the plot
    auto ampl = cepgen::NachtmannAmplitudes(cepgen::ParametersList().set<int>("mode", i));
    if (ampl.mode() == cepgen::NachtmannAmplitudes::Mode::SM) {
      cepgen::utils::Graph1D gr_sm_cth_400gev("sm_cth_400gev", "400 GeV;cos#theta;d#sigma/dcos#theta (pb)"),
          gr_sm_cth_2400gev("sm_cth_2400gev", "2400 GeV;cos#theta;d#sigma/dcos#theta (pb)");
      for (double x = -1.; x <= 1.; x += 0.1) {
        CG_LOG << x;
        gr_sm_cth_400gev.addPoint(x, dsig_dcosth(ampl, 400., x));
        gr_sm_cth_2400gev.addPoint(x, dsig_dcosth(ampl, 2400., x));
      }
      if (!plotter.empty()) {
        auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
        cepgen::utils::Drawer::Mode dm;
        if (logy)
          dm |= cepgen::utils::Drawer::Mode::logy;
        if (draw_grid)
          dm |= cepgen::utils::Drawer::Mode::grid;
        plt->draw(gr_sm_cth_400gev, dm);
        plt->draw(gr_sm_cth_2400gev, dm);
      }
    }
  }

  return 0;
}
