/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
 *                2021       Sergey Kulagin
 *                           Vladislav Barinov
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
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/GSLDerivator.h"
#include "CepGen/Utils/GridHandler.h"

#define DEBUG_GRID

#ifdef DEBUG_GRID
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#endif

namespace cepgen {
  namespace strfun {
    /// Kulagin and Barinov hybrid parameterisation
    /// \cite Kulagin:2021mee
    class KulaginBarinov : public Parameterisation {
    public:
      explicit KulaginBarinov(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Kulagin-Barinov (hybrid)");
        desc.addParametersDescriptionVector(
            "resonances",
            ResonanceObject::description(),
            {ParametersList()  // Delta(1232)
                 .set<double>("mass", 1.2270)
                 .set<double>("width", 0.11279)
                 .set<int>("angularMomentum", 1)
                 .set<double>("x0", 0.055384)
                 .set<std::vector<double> >("a", {0.31115, 2.0294, 1.6713, 2.76})
                 .set<std::vector<double> >("c", {0.05029, 0., 0.42522})
                 .set<ParametersList>(
                     "branchingRatios",
                     ParametersList().set<double>("singlePi", 1.).set<double>("doublePi", 0.).set<double>("eta", 0.)),
             ParametersList()  // N(1440)
                 .set<double>("mass", 1.4497)
                 .set<double>("width", 0.40223)
                 .set<int>("angularMomentum", 1)
                 .set<double>("x0", 0.1125)
                 .set<std::vector<double> >("a", {0.089547, 0.18087, 0.23431, 4.1173})
                 .set<std::vector<double> >("c", {0., 0.23847, 1.4982})
                 .set<ParametersList>(
                     "branchingRatios",
                     ParametersList().set<double>("singlePi", 0.65).set<double>("doublePi", 0.35).set<double>("eta", 0.)),
             ParametersList()  // R1
                 .set<double>("mass", 1.5123)
                 .set<double>("width", 0.094542)
                 .set<int>("angularMomentum", 2)
                 .set<double>("x0", 0.4959)
                 .set<std::vector<double> >("a", {0.10677, 0.24897, 0.55621, 3.0798})
                 .set<std::vector<double> >("c", {0.091979, -0.10652, 1.0758})
                 .set<ParametersList>(
                     "branchingRatios",
                     ParametersList().set<double>("singlePi", 0.75).set<double>("doublePi", 0.25).set<double>("eta", 0.)),
             ParametersList()  // R2
                 .set<double>("mass", 1.5764)
                 .set<double>("width", 0.50046)
                 .set<int>("angularMomentum", 0)
                 .set<double>("x0", 0.30969)
                 .set<std::vector<double> >("a", {0.38953, -0.17962, 0.37638, 2.9622})
                 .set<std::vector<double> >("c", {0., 0., 0.})
                 .set<ParametersList>(
                     "branchingRatios",
                     ParametersList().set<double>("singlePi", 0.15).set<double>("doublePi", 0.85).set<double>("eta", 0.)),
             ParametersList()  // R3
                 .set<double>("mass", 1.7002)
                 .set<double>("width", 0.11768)
                 .set<int>("angularMomentum", 2)
                 .set<double>("x0", 0.25831)
                 .set<std::vector<double> >("a", {0.067075, 0.097330, 0.27891, 3.5372})
                 .set<std::vector<double> >("c", {0.12027, 0., 0.89367})
                 .set<ParametersList>("branchingRatios",
                                      ParametersList()
                                          .set<double>("singlePi", 0.15)
                                          .set<double>("doublePi", 0.6)
                                          .set<double>("eta", 0.25))});
        // DIS block
        auto dis_desc = ParametersDescription();
        dis_desc.add<double>("bg1l", 3.4742);
        dis_desc.add<double>("bg2l", 0.54193);
        dis_desc.add<double>("pml", 1.1).setDescription("exponent of t dependence for FL");
        dis_desc.add<double>("bg1t", 0.14453);
        dis_desc.add<double>("bg2t", 3.1297);
        dis_desc.add<double>("pmt", 1.6302).setDescription("exponent of t dependence for FT");
        desc.add<ParametersDescription>("disParameters", dis_desc);

        desc.add<double>("eps1", 1.e-6);
        desc.add<double>("t0", 2.);
        desc.add<std::string>("gridFile", "a08tmc.dat");
        return desc;
      }

      KulaginBarinov& eval(double xbj, double q2) override;

    private:
      static constexpr double prefactor_ = M_1_PI * M_1_PI / constants::ALPHA_EM;
      const double eps1_, t0_;
      const std::string sfs_grid_file_;
      enum Polarisation { L, T };
      class Resonance : public ResonanceObject, public SteeredObject<Resonance> {
      public:
        explicit Resonance(const ParametersList& params)
            : ResonanceObject(params),
              SteeredObject<Resonance>(params),
              a_(SteeredObject<Resonance>::steer<std::vector<double> >("a")),
              c_(SteeredObject<Resonance>::steer<std::vector<double> >("c")) {}

        static ParametersDescription description() {
          auto desc = ResonanceObject::description();
          desc.add<std::vector<double> >("a", std::vector<double>(4, 0.));
          desc.add<std::vector<double> >("c", std::vector<double>(3, 0.));
          return desc;
        }

        bool computeStrFuns(const KinematicsBlock& kin, double& fl, double& ft) const {
          // compute contributions to the total resonance width
          double width_t = partialWidth(kin);
          if (width_t <= 0.)
            return false;
          // off-shell effect on electrocouplings
          double f_gamma = photonWidth(kin);
          const double mass2 = mass_ * mass_;

          // Breit-Wigner factor together with off-shell factor
          double f_bw = f_gamma * kcmr() * mass2 * (width_t / width_) /
                        (std::pow(kin.w2 - mass2, 2) + mass2 * std::pow(width_t, 2));

          // compute structure functions using model of resonance helicity amplitudes
          fl = f_bw * std::pow((c_.at(0) + c_.at(1) * kin.q2) * exp(-c_.at(2) * kin.q2), 2);
          ft = f_bw * std::pow((a_.at(0) + a_.at(1) * kin.q2) * std::pow(1. + a_.at(2) * kin.q2, -a_.at(3)), 2);
          return true;
        }

      private:
        const std::vector<double> a_;
        const std::vector<double> c_;
      };
      std::vector<Resonance> resonances_;
      const struct DISParameters : SteeredObject<DISParameters> {
        explicit DISParameters(const ParametersList& params)
            : SteeredObject(params),
              bg1l(steer<double>("bg1l")),
              bg2l(steer<double>("bg2l")),
              pml(steer<double>("pml")),
              bg1t(steer<double>("bg1t")),
              bg2t(steer<double>("bg2t")),
              pmt(steer<double>("pmt")) {}
        double bg1l, bg2l, pml;
        double bg1t, bg2t, pmt;
      } dis_params_;
      GridHandler<2, 2> sfs_grid_{GridType::linear};
      utils::GSLDerivator deriv_;
      const double mpi2_, meta2_;
    };

    KulaginBarinov::KulaginBarinov(const ParametersList& params)
        : Parameterisation(params),
          eps1_(steer<double>("eps1")),
          t0_(steer<double>("t0")),
          sfs_grid_file_(steer<std::string>("gridFile")),
          dis_params_(steer<ParametersList>("disParameters")),
          deriv_(steer<ParametersList>("derivator")),
          mpi2_(std::pow(PDG::get().mass(PDG::piZero), 2)),
          meta2_(std::pow(PDG::get().mass(PDG::eta), 2)) {
      for (const auto& res : steer<std::vector<ParametersList> >("resonances"))
        resonances_.emplace_back(res);
      {  // build the FT and F2 grid
        if (!utils::fileExists(sfs_grid_file_))
          throw CG_FATAL("KulaginBarinov")
              << "Failed to load the DIS structure functions interpolation grid from '" << sfs_grid_file_ << "'!";
        CG_INFO("KulaginBarinov") << "Loading A08 structure function values from '" << sfs_grid_file_ << "' file.";
        std::ifstream grid_file(sfs_grid_file_);
        static const size_t num_xbj = 99, num_q2 = 70, num_sf = 2;
        static const double min_xbj = 1.01e-5, min_q2 = 0.8, max_q2 = 1.e3;
        //--- xbj & Q2 binning
        const size_t nxbb = num_xbj / 2;
        const double x1 = 0.3, xlog1 = log(x1), delx = (xlog1 - log(min_xbj)) / (nxbb - 1),
                     delx1 = std::pow(1. - x1, 2) / (nxbb + 1);
        const double dels = (log(log(max_q2 / 0.04)) - log(log(min_q2 / 0.04))) / (num_q2 - 1);
        for (size_t idx_xbj = 0; idx_xbj < num_xbj; ++idx_xbj) {  // xbj grid
          const double xbj = idx_xbj < nxbb ? exp(log(min_xbj) + delx * idx_xbj)
                                            : 1. - std::sqrt(fabs(std::pow(1. - x1, 2) - delx1 * (idx_xbj - nxbb + 1)));
          for (size_t idx_q2 = 0; idx_q2 < num_q2; ++idx_q2) {  // Q^2 grid
            const double q2 = 0.04 * exp(exp(log(log(min_q2 / 0.04)) + dels * idx_q2));
            std::array<double, num_sf> sfs;
            for (size_t idx_sf = 0; idx_sf < num_sf; ++idx_sf)
              grid_file >> sfs[idx_sf];
            CG_DEBUG("KulaginBarinov:grid") << "Inserting new values into grid: " << std::vector<double>{xbj, q2} << "("
                                            << std::vector<size_t>{idx_xbj, idx_q2} << "): " << sfs;
            sfs_grid_.insert({xbj, q2}, sfs);
          }
        }
        sfs_grid_.init();
        CG_DEBUG("KulaginBarinov:grid") << "Grid boundaries: " << sfs_grid_.boundaries();
#ifdef DEBUG_GRID
        utils::DrawerFactory::get().build("root")->draw(sfs_grid_, utils::Drawer::Mode::logz);
#endif
      }
    }

    KulaginBarinov& KulaginBarinov::eval(double xbj, double q2) {
      const double w2 = utils::mX2(xbj, q2, mp2_), w = std::sqrt(w2);
      double fl{0.}, f2{0.};
      const double gamma2 = 1. + tau(xbj, q2);
      {  //--- resonances region
        const auto kin = Resonance::KinematicsBlock(w2, q2, mp2_, mpi2_, meta2_);
        // proton c.m. energy and momentum (needed to compute additional kinematics factor for FL)
        const double ecm = utils::energyFromW(w, -q2, mp2_), q2cm = ecm * ecm - mp2_;
        double fl_res{0.}, ft_res{0.};
        for (const auto& res : resonances_) {  // sum over the resonant contributions
          double fl_sng_res, ft_sng_res;
          if (!res.computeStrFuns(kin, fl_sng_res, ft_sng_res)) {
            setFL(0.);
            setF2(0.);
            return *this;
          }
          fl_res += fl_sng_res;
          ft_res += ft_sng_res;
        }  // end of resonance loop

        // finalize calculation of sfn and xsec taking into account normalization factors
        ft_res *= prefactor_ * xbj * mp_;
        fl_res *= 2. * prefactor_ * xbj * mp_ * q2 / q2cm;
        const double f2_res = (fl_res + ft_res) / gamma2;
        fl += fl_res;
        f2 += f2_res;
      }  //--- end of resonances region
      {  //--- perturbative region
        const auto sfs = sfs_grid_.eval({xbj, q2});
        double f2_dis = sfs.at(0), fl_dis = sfs.at(1);
        if (q2 > t0_) {
        } else if (q2 > 1.e-12) {
          const double t = std::max(q2, t0_);
          // make extrapolation in q2 from q2=t down to q2=0
          // compute derivative of SF with respect to q2 at q2=t

          double ddt{0.}, ddl{0.};
          if (xbj >= 1.e-6 && q2 > 0.8) {  // check the range of validity
            // parameterisation of Twist-4 correction from A08 analysis arXiv:0710.0124 [hep-ph] (assuming F2ht=FTht)
            auto sfnht = [](double xbj, double q2) -> double {
              return (std::pow(xbj, 0.9) * std::pow(1. - xbj, 3.63) * (xbj - 0.356) *
                      (1.0974 + 47.7352 * std::pow(xbj, 4))) /
                     q2;
            };
            // DIS structure function model using the results of A08 analysis arXiv:0710.0124 [hep-ph]
            ddt = deriv_.eval(
                [&](double q2) -> double {
                  return sfs_grid_.eval({xbj, q2}).at(0) + sfnht(xbj, q2);
                },  // TM-corrected FT with twist-4 correction
                q2,
                q2 * 1.e-2);
            ddl = deriv_.eval(
                [&](double q2) -> double {
                  const auto vals = sfs_grid_.eval({xbj, q2});  // FT, F2
                  return gamma2 * (vals.at(1) + sfnht(xbj, q2)) - (vals.at(0) + sfnht(xbj, q2));
                },  // TM-corrected FL with twist-4 correction
                q2,
                q2 * 1.e-2);
          }
          const double fl_der = q2 * (std::pow(q2, dis_params_.pml - 1.) / std::pow(t, dis_params_.pml) *
                                      (fl_dis + log(t / q2) * (dis_params_.pml * fl_dis - t * ddl)));

          /// Regge fit to total photoproduction cross section from hep-ph/9908218
          /// \param[in] w2 invariant mass squared (in GeV^2)
          /// \return cross section value in mb
          auto photot = [](double w2) -> double {
            return 0.0598 * std::pow(w2, 0.0933) + 0.1164 * std::pow(w2, -0.357);
          };
          const double ft_dis = gamma2 * f2_dis - fl_dis;
          // extrapolation in q2 ; note cross section in mb units, and converted to GeV units
          const double f0t = photot(w2) / constants::G_EM_SQ * M_1_PI / (constants::GEVM2_TO_PB * 1.e-9);
          const double ft_der =
              q2 * (f0t + std::pow(q2, dis_params_.pmt - 1.) / std::pow(t, dis_params_.pmt) *
                              (ft_dis - f0t * t +
                               log(t / q2) * (dis_params_.pmt * ft_dis - t * ddt - (dis_params_.pmt - 1.) * f0t * t)));
          const auto f2_der = (fl_der + ft_der) / gamma2;

          fl_dis = fl_der;
          f2_dis = f2_der;
        } else  // Q^2 -> 0 limit
          f2_dis = fl_dis = 0.;
        const double bl = 1. - exp(-dis_params_.bg1l * std::pow(w2 - mx_min_, dis_params_.bg2l)),
                     bt = 1. - exp(-dis_params_.bg1t * std::pow(w2 - mx_min_, dis_params_.bg2t));
        fl += fl_dis * bl;
        f2 += f2_dis * bt;
      }  //--- end of perturbative region
      setFL(fl);
      setF2(f2);
      return *this;
    }
  }  // namespace strfun
}  // namespace cepgen
REGISTER_STRFUN(303, KulaginBarinov, strfun::KulaginBarinov)
