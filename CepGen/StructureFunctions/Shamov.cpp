/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/GridHandler.h"

namespace cepgen {
  namespace strfun {
    class Shamov final : public Parameterisation {
    public:
      explicit Shamov(const ParametersList&);

      static ParametersDescription description();

      Shamov& eval(double xbj, double q2) override;

    private:
      static constexpr double prefac_ =
          2. * M_PI * M_PI * constants::ALPHA_EM * constants::GEVM2_TO_PB * 1e-9;  // pb/GeV -> mb/GeV
      static constexpr std::array<float, 19> gmq_ = {{0.000,
                                                      0.200,
                                                      0.300,
                                                      0.400,
                                                      0.470,
                                                      0.480,
                                                      0.500,
                                                      0.600,
                                                      0.630,
                                                      0.631,
                                                      0.770,
                                                      0.780,
                                                      0.790,
                                                      0.970,
                                                      0.980,
                                                      1.150,
                                                      1.340,
                                                      1.570,
                                                      2.340}};
      static constexpr std::array<float, 19> gmv_ = {{3.000,
                                                      1.770,
                                                      1.380,
                                                      1.170,
                                                      0.978,
                                                      0.961,
                                                      0.964,
                                                      0.766,
                                                      0.735,
                                                      0.719,
                                                      0.570,
                                                      0.572,
                                                      0.553,
                                                      0.460,
                                                      0.446,
                                                      0.326,
                                                      0.269,
                                                      0.209,
                                                      0.102}};
      static constexpr std::array<float, 291> gp_en_ = {
          {1.070, 1.100, 1.130, 1.150, 1.174, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.306, 1.324, 1.342, 1.359,
           1.376, 1.393, 1.410, 1.426, 1.443, 1.456, 1.459, 1.475, 1.491, 1.506, 1.522, 1.537, 1.543, 1.552, 1.567,
           1.582, 1.597, 1.612, 1.617, 1.626, 1.640, 1.655, 1.669, 1.683, 1.697, 1.710, 1.724, 1.738, 1.743, 1.751,
           1.764, 1.778, 1.791, 1.804, 1.817, 1.830, 1.831, 1.843, 1.855, 1.868, 1.880, 1.893, 1.898, 1.905, 1.917,
           1.922, 1.930, 1.942, 1.954, 1.966, 1.978, 1.989, 2.001, 2.013, 2.025, 2.036, 2.041, 2.048, 2.059, 2.070,
           2.082, 2.093, 2.104, 2.115, 2.117, 2.126, 2.137, 2.148, 2.153, 2.154, 2.159, 2.170, 2.174, 2.181, 2.191,
           2.202, 2.213, 2.223, 2.234, 2.238, 2.244, 2.255, 2.265, 2.275, 2.286, 2.296, 2.300, 2.306, 2.316, 2.320,
           2.326, 2.336, 2.346, 2.356, 2.360, 2.366, 2.376, 2.386, 2.396, 2.400, 2.406, 2.415, 2.419, 2.425, 2.435,
           2.444, 2.454, 2.464, 2.473, 2.477, 2.483, 2.492, 2.501, 2.511, 2.513, 2.520, 2.529, 2.533, 2.539, 2.548,
           2.551, 2.557, 2.566, 2.575, 2.584, 2.593, 2.602, 2.611, 2.620, 2.624, 2.629, 2.638, 2.642, 2.647, 2.649,
           2.656, 2.665, 2.674, 2.682, 2.691, 2.695, 2.700, 2.708, 2.717, 2.726, 2.734, 2.743, 2.746, 2.751, 2.760,
           2.763, 2.768, 2.777, 2.785, 2.794, 2.797, 2.802, 2.810, 2.819, 2.827, 2.830, 2.835, 2.844, 2.847, 2.852,
           2.860, 2.868, 2.876, 2.877, 2.885, 2.893, 2.896, 2.897, 2.901, 2.909, 2.917, 2.919, 2.925, 2.933, 2.941,
           2.944, 2.949, 2.957, 2.958, 2.965, 3.032, 3.035, 3.038, 3.114, 3.130, 3.147, 3.207, 3.218, 3.253, 3.296,
           3.305, 3.383, 3.389, 3.416, 3.471, 3.479, 3.540, 3.551, 3.582, 3.634, 3.683, 3.784, 3.867, 3.868, 3.911,
           3.937, 4.029, 4.087, 4.144, 4.204, 4.257, 4.282, 4.327, 4.334, 4.379, 4.396, 4.458, 4.514, 4.580, 4.581,
           4.645, 4.715, 4.773, 4.860, 4.912, 4.931, 5.001, 5.065, 5.140, 5.141, 5.213, 5.356, 5.374, 5.541, 5.623,
           5.704, 5.862, 5.935, 6.665, 7.271, 7.672, 7.733, 8.066, 8.329, 8.485, 9.125, 9.186, 9.576, 9.979, 15.10,
           20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0, 120.0, 150.0, 170.0, 200.0, 240.0, 270.0,
           310.0, 350.0, 400.0, 500.0, 7000., 14000.}};
      static constexpr std::array<float, 291> gp_cs_ = {
          {0.000, 0.052, 0.108, 0.175, 0.424, 0.487, 0.527, 0.478, 0.407, 0.334, 0.244, 0.225, 0.200, 0.178, 0.177,
           0.187, 0.194, 0.212, 0.223, 0.233, 0.211, 0.240, 0.265, 0.279, 0.276, 0.261, 0.245, 0.201, 0.221, 0.206,
           0.214, 0.209, 0.202, 0.193, 0.205, 0.201, 0.212, 0.218, 0.215, 0.192, 0.191, 0.175, 0.165, 0.182, 0.159,
           0.162, 0.150, 0.149, 0.144, 0.156, 0.150, 0.147, 0.154, 0.154, 0.154, 0.147, 0.154, 0.154, 0.144, 0.152,
           0.151, 0.156, 0.154, 0.146, 0.139, 0.156, 0.150, 0.150, 0.145, 0.139, 0.145, 0.146, 0.142, 0.141, 0.142,
           0.143, 0.149, 0.154, 0.135, 0.115, 0.148, 0.144, 0.144, 0.143, 0.142, 0.149, 0.144, 0.148, 0.138, 0.132,
           0.145, 0.138, 0.145, 0.136, 0.138, 0.138, 0.139, 0.136, 0.129, 0.136, 0.140, 0.144, 0.133, 0.139, 0.139,
           0.143, 0.140, 0.140, 0.139, 0.134, 0.141, 0.130, 0.136, 0.124, 0.132, 0.128, 0.130, 0.142, 0.132, 0.134,
           0.139, 0.133, 0.144, 0.133, 0.135, 0.136, 0.130, 0.134, 0.134, 0.111, 0.130, 0.131, 0.130, 0.129, 0.140,
           0.134, 0.138, 0.129, 0.144, 0.128, 0.133, 0.132, 0.127, 0.128, 0.127, 0.124, 0.124, 0.133, 0.127, 0.127,
           0.121, 0.134, 0.129, 0.134, 0.123, 0.129, 0.132, 0.121, 0.137, 0.123, 0.130, 0.135, 0.127, 0.129, 0.128,
           0.126, 0.123, 0.122, 0.120, 0.119, 0.128, 0.134, 0.132, 0.125, 0.122, 0.126, 0.114, 0.127, 0.135, 0.122,
           0.140, 0.125, 0.136, 0.116, 0.136, 0.124, 0.131, 0.127, 0.146, 0.136, 0.120, 0.131, 0.142, 0.132, 0.129,
           0.132, 0.134, 0.133, 0.133, 0.127, 0.128, 0.116, 0.125, 0.127, 0.124, 0.122, 0.134, 0.122, 0.128, 0.122,
           0.122, 0.130, 0.118, 0.125, 0.124, 0.121, 0.116, 0.122, 0.129, 0.122, 0.120, 0.118, 0.126, 0.122, 0.124,
           0.124, 0.122, 0.120, 0.119, 0.123, 0.115, 0.124, 0.114, 0.123, 0.117, 0.120, 0.122, 0.128, 0.124, 0.112,
           0.128, 0.121, 0.122, 0.115, 0.115, 0.118, 0.116, 0.111, 0.115, 0.118, 0.119, 0.113, 0.113, 0.115, 0.114,
           0.113, 0.115, 0.117, 0.115, 0.114, 0.114, 0.114, 0.115, 0.112, 0.114, 0.115, 0.115, 0.115, 0.114, 0.113,
           0.114, 0.115, 0.115, 0.114, 0.114, 0.117, 0.116, 0.116, 0.116, 0.114, 0.117, 0.116, 0.114, 0.116, 0.120,
           0.116, 0.118, 0.143, 0.143, 0.143, 0.143}};
      static constexpr std::array<float, 291> gp_nr_ = {
          {0.683, 0.539, 0.374, 0.267, 0.378, 0.143, 0.085, 0.137, 0.233, 0.297, 0.249, 0.333, 0.364, 0.365, 0.408,
           0.460, 0.476, 0.494, 0.472, 0.431, 0.303, 0.374, 0.376, 0.384, 0.394, 0.412, 0.444, 0.357, 0.458, 0.486,
           0.553, 0.576, 0.577, 0.559, 0.580, 0.545, 0.530, 0.531, 0.568, 0.590, 0.655, 0.677, 0.699, 0.740, 0.720,
           0.750, 0.749, 0.765, 0.771, 0.800, 0.803, 0.799, 0.817, 0.825, 0.831, 0.830, 0.844, 0.846, 0.839, 0.852,
           0.853, 0.860, 0.863, 0.859, 0.856, 0.876, 0.874, 0.877, 0.876, 0.873, 0.881, 0.883, 0.881, 0.883, 0.886,
           0.889, 0.895, 0.901, 0.889, 0.870, 0.900, 0.900, 0.901, 0.901, 0.900, 0.906, 0.904, 0.908, 0.902, 0.899,
           0.909, 0.906, 0.912, 0.908, 0.910, 0.910, 0.912, 0.911, 0.908, 0.914, 0.917, 0.920, 0.914, 0.919, 0.919,
           0.922, 0.921, 0.922, 0.923, 0.920, 0.925, 0.919, 0.924, 0.917, 0.923, 0.921, 0.923, 0.929, 0.925, 0.927,
           0.930, 0.927, 0.934, 0.929, 0.930, 0.931, 0.929, 0.931, 0.932, 0.918, 0.931, 0.932, 0.931, 0.931, 0.937,
           0.935, 0.937, 0.933, 0.940, 0.934, 0.937, 0.937, 0.935, 0.936, 0.936, 0.935, 0.935, 0.939, 0.937, 0.937,
           0.935, 0.941, 0.939, 0.942, 0.937, 0.941, 0.942, 0.937, 0.945, 0.939, 0.943, 0.945, 0.942, 0.943, 0.943,
           0.942, 0.941, 0.941, 0.941, 0.941, 0.945, 0.948, 0.947, 0.944, 0.943, 0.946, 0.940, 0.946, 0.950, 0.945,
           0.952, 0.947, 0.951, 0.943, 0.952, 0.947, 0.950, 0.948, 0.955, 0.952, 0.946, 0.951, 0.955, 0.952, 0.951,
           0.952, 0.953, 0.953, 0.953, 0.951, 0.954, 0.949, 0.953, 0.956, 0.955, 0.955, 0.960, 0.957, 0.960, 0.959,
           0.959, 0.963, 0.960, 0.962, 0.963, 0.962, 0.962, 0.964, 0.966, 0.966, 0.966, 0.967, 0.970, 0.969, 0.970,
           0.971, 0.971, 0.971, 0.972, 0.973, 0.972, 0.974, 0.973, 0.975, 0.974, 0.975, 0.976, 0.977, 0.977, 0.975,
           0.978, 0.978, 0.978, 0.978, 0.978, 0.979, 0.979, 0.978, 0.979, 0.980, 0.981, 0.980, 0.980, 0.982, 0.982,
           0.982, 0.983, 0.984, 0.986, 0.987, 0.988, 0.988, 0.989, 0.989, 0.990, 0.991, 0.991, 0.991, 0.992, 0.992,
           0.992, 0.993, 0.993, 0.993, 0.993, 0.994, 0.994, 0.994, 0.994, 0.995, 0.995, 0.995, 0.995, 0.995, 0.996,
           0.996, 0.996, 1.000, 1.000, 1.000, 1.000}};

      const enum class Mode {
        SuriYennie = 0,  ///< Suri & Yennie
        /// real photon cross section for q^2=0,
        /// q^2 dependence as for Delta(1232)
        RealRes = 1,
        /// real photon cross section for q^2=0,
        /// q^2 dependence: resonant contribution as for Delta(1232),
        /// nonresonant contribution according to S&Y
        RealResAndNonRes = 2,
        /// real photon cross section for q^2=0,
        /// q^2 dependence according to S&Y
        RealAndSuriYennieNonRes = 3,
        /// real photon cross section for q^2=0,
        /// q^2 dependence: resonant contribution as for Delta(1232),
        /// some fit for nonresonant contribution
        RealAndFitNonRes = 4
      } mode_;
      size_t fit_model_;
      const double gm0_;
      const double gmb_;
      const double q20_;
      const double r_power_;
      const double lowq2_;

      GridHandler<1, 2> sigma_grid_{GridType::linear};
      GridHandler<1, 1> gm_grid_{GridType::linear};
      std::unique_ptr<Parameterisation> sy_sf_;
      bool non_resonant_{true};
    };

    constexpr std::array<float, 19> Shamov::gmq_, Shamov::gmv_;
    constexpr std::array<float, 291> Shamov::gp_en_, Shamov::gp_cs_, Shamov::gp_nr_;

    Shamov::Shamov(const ParametersList& params)
        : Parameterisation(params),
          mode_(steerAs<int, Mode>("mode")),
          fit_model_(steer<int>("fitModel")),
          gm0_(steer<double>("gm0")),
          gmb_(steer<double>("gmb")),
          q20_(steer<double>("q20")),
          r_power_(steer<double>("rPower")),
          lowq2_(steer<double>("lowQ2")),
          sy_sf_(StructureFunctionsFactory::get().build((int)strfun::Type::SuriYennie,
                                                        steer<ParametersList>("syParams"))) {
      //----- initialise the interpolation grids

      //--- grid E -> (cross section, norm)
      for (size_t i = 0; i < gp_en_.size(); ++i)
        sigma_grid_.insert({gp_en_.at(i)}, {gp_cs_.at(i), gp_nr_.at(i)});
      sigma_grid_.init();

      //--- grid Q -> gamma_v
      for (size_t i = 0; i < gmv_.size(); ++i)
        gm_grid_.insert({gmq_.at(i)}, {gmv_.at(i)});
      gm_grid_.init();
      if (mode_ == Mode::SuriYennie || mode_ == Mode::RealAndSuriYennieNonRes || mode_ == Mode::RealResAndNonRes ||
          mode_ == Mode::RealAndFitNonRes)
        non_resonant_ = true;
    }

    Shamov& Shamov::eval(double xbj, double q2) {
      //--- Suri & Yennie structure functions
      if (mode_ == Mode::SuriYennie) {  // trivial composite case
        Parameterisation::operator=((*sy_sf_)(xbj, q2));
        return *this;
      }

      const double mx2 = utils::mX2(xbj, q2, mp2_), mx = sqrt(mx2);

      const auto sigma = sigma_grid_.eval({mx});
      double sgp = sigma[0];  // cross section value at MX

      if (mode_ == Mode::RealAndFitNonRes && mx > 1.5)
        non_resonant_ = true;
      //if (mx > sigma_grid_.max()[0] || rand() * 1. / RAND_MAX < sigma[1])
      //  non_resonant_ = true;

      //--- q^2 dependence of sigma
      double Gm;
      if (non_resonant_) {
        if (mode_ == Mode::RealAndFitNonRes)
          /// q^2 dependence of the non-resonant gamma-p cross section
          /// Some fit of the know e-p data. Good only for Q^2 < 6 GeV^2
          Gm = std::pow(1. + std::pow(q2 / q20_, 2), -r_power_);
        else
          Gm = sy_sf_->W1(xbj, q2) / sy_sf_->W1(xbj, lowq2_);
      } else {                        // resonant
        if (q2 >= gm_grid_.max()[0])  // above grid range
          Gm = gm0_ * exp(-gmb_ * q2);
        else
          Gm = gm_grid_.eval({q2})[0];
        Gm /= 3.;  // due to data normalization
      }
      sgp *= Gm;  // cross section with some q^2 dependence

      //--- for W -> cross section
      const double s1 = prefac_ * 4. * mp_ / (mx2 - mp2_);  // mb/GeV
      const double s2 =
          prefac_ * (std::pow(mx2 - mp2_, 2) + 2. * (mx2 + mp2_) * q2 + q2 * q2) / mp_ / q2 / (mx2 - mp2_);  // mb/GeV
      //---  ratio (\sigma_T+\sigma_L)/\sigma_T according to S & Y
      const double ratio = (s2 * sy_sf_->W2(xbj, q2)) / (s1 * sy_sf_->W1(xbj, q2));

      if (mode_ == Mode::RealAndFitNonRes) {  // transverse sigma only
        setW1(sgp / s1);
        setW2(sgp * ratio / s2);
      } else {  // longitudinal + transverse component for sigma
        setW1(sgp / ratio / s1);
        setW2(sgp / s2);
      }
      const double nu = 0.5 * (q2 + mx2 - mp2_) / mp_;
      setF2(W2(xbj, q2) * nu / mp_);
      return *this;
    }

    ParametersDescription Shamov::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Shamov (hybrid, soft)");
      desc.addAs<int, Mode>("mode", Mode::RealResAndNonRes).setDescription("sub-structure functions choice");
      desc.add<int>("fitModel", 2);
      desc.add<double>("q20", 0.65 /* 0.36 */)
          .setDescription("first parameter for non-resonant gamma-p cross section q^2 dependence");
      desc.add<double>("rPower", 0.71 /* 0.52 */)
          .setDescription("second parameter for non-resonant gamma-p cross section q^2 dependence");
      desc.add<double>("gm0", 1.).setDescription("scaling factor for the magnetic form factor template fit");
      desc.add<double>("gmb", 0.984).setDescription("exponential parameter for the magnetic form factor template fit");
      desc.add<double>("lowQ2", 1.e-7);
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::Shamov, Shamov, strfun::Shamov);
