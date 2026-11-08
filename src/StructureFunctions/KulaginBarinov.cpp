/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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
#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/GridHandler.h"

using namespace std::string_literals;

namespace cepgen::strfun {
  /// Kulagin and Barinov hybrid parameterisation
  /// \cite Kulagin:2021mee
  class KulaginBarinov final : public Parameterisation {
  public:
    explicit KulaginBarinov(const ParametersList& params)
        : Parameterisation(params),
          t0_(steer<double>("t0")),
          q2_range_(steer<Limits>("Q2range")),
          q2_grid_range_(steer<Limits>("Q2gridRange")),
          sfs_grid_file_(steerPath("gridFile")),
          dis_params_(steer<ParametersList>("disParameters")),
          derivator_(DerivatorFactory::get().build(steer<ParametersList>("derivator"))),
          mpi2_(std::pow(PDG::get().mass(PDG::piZero), 2)),
          meta2_(std::pow(PDG::get().mass(PDG::eta), 2)) {
      for (const auto& resonance : steer<std::vector<ParametersList> >("resonances"))
        resonances_.emplace_back(resonance);
      {  // build the FT and F2 grid
        if (!utils::fileExists(sfs_grid_file_))
          throw CG_FATAL("KulaginBarinov")
              << "Failed to load the DIS structure functions interpolation grid from '" << sfs_grid_file_ << "'!";
        CG_INFO("KulaginBarinov") << "Loading A08 structure function values from '" << sfs_grid_file_ << "' file.";
        std::ifstream grid_file(sfs_grid_file_);
        static constexpr size_t num_xbj = 99, num_q2 = 70, num_sf = 2;
        static constexpr double min_xbj = 1.01e-5;
        //--- xbj & Q2 binning
        constexpr size_t num_bins_xbj = num_xbj / 2;
        const double x1 = 0.3, log_x1 = std::log(x1), delta_x = (log_x1 - std::log(min_xbj)) / (num_bins_xbj - 1),
                     delta_x1 = std::pow(1. - x1, 2) / (num_bins_xbj + 1);
        const double deltas =
            (std::log(std::log(q2_grid_range_.max() / 0.04)) - std::log(std::log(q2_grid_range_.min() / 0.04))) /
            (num_q2 - 1);
        // parameterisation of Twist-4 correction from A08 analysis arXiv:0710.0124 [hep-ph] (assuming F2ht=FTht)
        auto sf_higher_twist = [](double xbj, double q2) -> double {
          return (std::pow(xbj, 0.9) * std::pow(1. - xbj, 3.63) * (xbj - 0.356) *
                  (1.0974 + 47.7352 * std::pow(xbj, 4))) /
                 q2;
        };

        for (size_t idx_xbj = 0; idx_xbj < num_xbj; ++idx_xbj) {  // xbj grid
          const double xbj =
              idx_xbj < num_bins_xbj
                  ? std::exp(log(min_xbj) + delta_x * idx_xbj)
                  : 1. - std::sqrt(std::fabs(std::pow(1. - x1, 2) - delta_x1 * (idx_xbj - num_bins_xbj + 1)));
          for (size_t idx_q2 = 0; idx_q2 < num_q2; ++idx_q2) {  // Q^2 grid
            const double q2 =
                0.04 * std::exp(std::exp(std::log(std::log(q2_grid_range_.min() / 0.04)) + deltas * idx_q2));
            std::array<double, num_sf> sfs{};
            for (size_t idx_sf = 0; idx_sf < num_sf; ++idx_sf) {
              grid_file >> sfs[idx_sf];  // FT, F2
              sfs[idx_sf] += sf_higher_twist(xbj, q2);
            }
            CG_DEBUG("KulaginBarinov:grid") << "Inserting new values into grid: " << std::vector{xbj, q2} << "("
                                            << std::vector{idx_xbj, idx_q2} << "): " << sfs;
            sfs_grid_.insert({xbj, q2}, sfs);
          }
        }
        sfs_grid_.initialise();
        CG_DEBUG("KulaginBarinov:grid") << "Grid boundaries: " << sfs_grid_.boundaries();
      }
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Kulagin-Barinov (hybrid)");
      desc.add("derivator", DerivatorFactory::get().describeParameters("gsl"));
      desc.addParametersDescriptionVector(
          "resonances",
          ResonanceObject::description(),
          {ParametersList()  // Delta(1232)
               .set("mass", 1.2270)
               .set("width", 0.11279)
               .set("angularMomentum", 1)
               .set("x0", 0.055384)
               .set("a", std::vector{0.31115, 2.0294, 1.6713, 2.76})
               .set("c", std::vector{0.05029, 0., 0.42522})
               .set("branchingRatios", ParametersList().set("singlePi", 1.).set("doublePi", 0.).set("eta", 0.)),
           ParametersList()  // N(1440)
               .set("mass", 1.4497)
               .set("width", 0.40223)
               .set("angularMomentum", 1)
               .set("x0", 0.1125)
               .set("a", std::vector{0.089547, 0.18087, 0.23431, 4.1173})
               .set("c", std::vector{0., 0.23847, 1.4982})
               .set("branchingRatios", ParametersList().set("singlePi", 0.65).set("doublePi", 0.35).set("eta", 0.)),
           ParametersList()  // R1
               .set("mass", 1.5123)
               .set("width", 0.094542)
               .set("angularMomentum", 2)
               .set("x0", 0.4959)
               .set("a", std::vector{0.10677, 0.24897, 0.55621, 3.0798})
               .set("c", std::vector{0.091979, -0.10652, 1.0758})
               .set("branchingRatios", ParametersList().set("singlePi", 0.75).set("doublePi", 0.25).set("eta", 0.)),
           ParametersList()  // R2
               .set("mass", 1.5764)
               .set("width", 0.50046)
               .set("angularMomentum", 0)
               .set("x0", 0.30969)
               .set("a", std::vector{0.38953, -0.17962, 0.37638, 2.9622})
               .set("c", std::vector{0., 0., 0.})
               .set("branchingRatios", ParametersList().set("singlePi", 0.15).set("doublePi", 0.85).set("eta", 0.)),
           ParametersList()  // R3
               .set("mass", 1.7002)
               .set("width", 0.11768)
               .set("angularMomentum", 2)
               .set("x0", 0.25831)
               .set("a", std::vector{0.067075, 0.097330, 0.27891, 3.5372})
               .set("c", std::vector{0.12027, 0., 0.89367})
               .set("branchingRatios", ParametersList().set("singlePi", 0.15).set("doublePi", 0.6).set("eta", 0.25))});
      // DIS block
      auto dis_desc = ParametersDescription();
      dis_desc.add("bg1l", 3.4742);
      dis_desc.add("bg2l", 0.54193);
      dis_desc.add("pml", 1.1).setDescription("exponent of t dependence for FL");
      dis_desc.add("bg1t", 0.14453);
      dis_desc.add("bg2t", 3.1297);
      dis_desc.add("pmt", 1.6302).setDescription("exponent of t dependence for FT");
      desc.add("disParameters", dis_desc);

      desc.add("t0", 2.);
      desc.add("Q2range", Limits{1.e-12, 1.e3});
      desc.add("Q2gridRange", Limits{0.8, 1.e3}).setDescription("Q^2 range covered by the grid");
      desc.add("gridFile", "a08tmc.dat"s).setDescription("path to the DIS grid");
      return desc;
    }

    void eval() override;

  private:
    static constexpr double prefactor_ = M_1_PI * M_1_PI / constants::ALPHA_EM;
    const double t0_;
    const Limits q2_range_, q2_grid_range_;
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
        desc.add("a", std::vector(4, 0.));
        desc.add("c", std::vector(3, 0.));
        return desc;
      }

      bool computeStructureFunctions(const KinematicsBlock& kin, double& fl, double& ft) const {
        // compute contributions to the total resonance width
        const double width_t = partialWidth(kin);
        if (width_t <= 0.)
          return false;
        // off-shell effect on electro-couplings
        const double f_gamma = photonWidth(kin) / width_;
        const double mass2 = mass_ * mass_;

        // Breit-Wigner factor together with off-shell factor
        const double f_bw =
            f_gamma * kcmr() * mass2 * width_t / (std::pow(kin.w2 - mass2, 2) + mass2 * std::pow(width_t, 2));

        // compute structure functions using model of resonance helicity amplitudes
        fl = f_bw * std::pow((c_.at(0) + c_.at(1) * kin.q2) * std::exp(-c_.at(2) * kin.q2), 2);
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
    std::unique_ptr<utils::Derivator> derivator_;
    const double mpi2_, meta2_;
  };

  void KulaginBarinov::eval() {
    const double w2 = utils::mX2(args_.xbj, args_.q2, mp2_), w = std::sqrt(w2);
    double fl{0.}, f2{0.};
    {  //--- resonances region
      const auto kin = Resonance::KinematicsBlock(w2, args_.q2, mp2_, mpi2_, meta2_);
      // proton c.m. energy and momentum (needed to compute additional kinematics factor for FL)
      const double ecm = utils::energyFromW(w, -args_.q2, mp2_), q2cm = ecm * ecm - mp2_;
      double fl_res{0.}, ft_res{0.};
      for (const auto& res : resonances_) {  // sum over the resonant contributions
        double fl_sng_res, ft_sng_res;
        if (!res.computeStructureFunctions(kin, fl_sng_res, ft_sng_res)) {
          setFL(0.);
          setF2(0.);
          return;
        }
        fl_res += fl_sng_res;
        ft_res += ft_sng_res;
      }  // end of resonance loop

      // finalize calculation of structure functions and x-section taking into account normalization factors
      ft_res *= prefactor_ * args_.xbj * mp_;
      fl_res *= 2. * prefactor_ * args_.xbj * mp_ * args_.q2 / q2cm;
      const double f2_res = (fl_res + ft_res) / gamma2(args_.xbj, args_.q2);
      fl += fl_res;
      f2 += f2_res;
    }  //--- end of resonances region
    {  //--- perturbative region
      double f2_dis = 0., fl_dis = 0.;
      if (q2_range_.contains(args_.q2)) {
        double ft_dis = 0.;
        const double t = std::max(args_.q2, t0_), xbj_t = utils::xBj(t, mp2_, w2), gam2 = gamma2(xbj_t, t);
        if (t > q2_grid_range_.min()) {
          const auto sfs = sfs_grid_.eval({xbj_t, t});  // FT, F2
          ft_dis = sfs.at(0);
          f2_dis = sfs.at(1);
          fl_dis = gam2 * f2_dis - ft_dis;
        }
        if (args_.q2 < t0_) {
          // make extrapolation in q2 from q2=t down to q2=0
          // compute derivative of SF with respect to q2 at q2=t

          double ddt{0.}, ddl{0.};
          if (xbj_t >= 1.e-6) {  // check the range of validity
            // DIS structure function model using the results of A08 analysis arXiv:0710.0124 [hep-ph]
            ddt = derivator_->derivate(
                [this, &xbj_t](double qsq) -> double {
                  return sfs_grid_.eval({xbj_t, qsq}).at(0);
                },  // TM-corrected FT with twist-4 correction
                t,
                t * 1.e-2);
            ddl = derivator_->derivate(
                [this, &xbj_t](double qsq) -> double {
                  const auto vals = sfs_grid_.eval({xbj_t, qsq});  // FT, F2
                  const auto &ft_l = vals.at(0), &f2_l = vals.at(1);
                  return gamma2(xbj_t, qsq) * f2_l - ft_l;
                },  // TM-corrected FL with twist-4 correction
                t,
                t * 1.e-2);
          }
          const double fl_der = args_.q2 * (std::pow(args_.q2, dis_params_.pml - 1.) / std::pow(t, dis_params_.pml) *
                                            (fl_dis + std::log(t / args_.q2) * (dis_params_.pml * fl_dis - t * ddl)));

          /// Regge fit to total photoproduction cross section from hep-ph/9908218
          /// \param[in] w2 invariant mass squared (in GeV^2)
          /// \return cross section value in mb
          auto photoproduction_t = [](double w2) -> double {
            return 0.0598 * std::pow(w2, 0.0933) + 0.1164 * std::pow(w2, -0.357);
          };
          // extrapolation in q2 ; note cross-section in mb units, and converted to GeV units
          const double f0t = photoproduction_t(w2) / constants::G_EM_SQ * M_1_PI / (constants::GEVM2_TO_PB * 1.e-9);
          const double ft_der =
              args_.q2 * (f0t + std::pow(args_.q2, dis_params_.pmt - 1.) / std::pow(t, dis_params_.pmt) *
                                    (ft_dis - f0t * t +
                                     std::log(t / args_.q2) *
                                         (dis_params_.pmt * ft_dis - t * ddt - (dis_params_.pmt - 1.) * f0t * t)));
          fl_dis = fl_der;
          ft_dis = ft_der;
        }
        const double bl = 1. - std::exp(-dis_params_.bg1l * std::pow(w2 - mx_min_, dis_params_.bg2l)),
                     bt = 1. - std::exp(-dis_params_.bg1t * std::pow(w2 - mx_min_, dis_params_.bg2t));
        fl_dis *= bl;
        ft_dis *= bt;
        f2_dis = (fl_dis + ft_dis) / gamma2(args_.xbj, args_.q2);
      }
      fl += fl_dis;
      f2 += f2_dis;
    }  //--- end of perturbative region
    setFL(fl);
    setF2(f2);
  }
}  // namespace cepgen::strfun
using cepgen::strfun::KulaginBarinov;
REGISTER_STRFUN("KulaginBarinov", 303, KulaginBarinov);
