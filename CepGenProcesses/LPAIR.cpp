/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Algebra.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;

/**
 * Analytic matrix element \cite Vermaseren:1982cz for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process
 * The integrand variables are mapped as:
 * - 0 = \f$t_1\f$, first incoming photon's virtuality
 * - 1 = \f$t_2\f$, second incoming photon's virtuality
 * - 2 = \f$s_2\f$ mapping
 * - 3 = yy4 = \f$\cos\left(\pi x_3\right)\f$ definition
 * - 4 = \f$w_4\f$, the two-photon system's invariant mass
 * - 5 = xx6 = \f$\frac{1}{2}\left(1-\cos\theta^{\rm CM}_6\right)\f$ definition (3D rotation of the first outgoing lepton with respect to the two-photon centre-of-mass system). If the \a nm_ optimisation flag is set this angle coefficient value becomes
 *   \f[\frac{1}{2}\left(\frac{a_{\rm map}}{b_{\rm map}}\frac{\beta-1}{\beta+1}+1\right)\f] with
 *   \f$a_{\rm map}=\frac{1}{2}\left(w_4-t_1-t_2\right)\f$, \f$b_{\rm map}=\frac{1}{2}\sqrt{\left(\left(w_4-t_1-t_2\right)^2-4t_1t_2\right)\left(1-4\frac{w_6}{w_4}\right)}\f$, and
 *   \f$\beta=\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)^{2x_5-1}\f$
 *   and the Jacobian element is scaled by a factor
 *   \f$\frac{1}{2}\frac{\left(a_{\rm map}^2-b_{\rm map}^2\cos^2\theta^{\rm CM}_6\right)}{a_{\rm map}b_{\rm map}}\log\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)\f$
 * - 6 = _phicm6_, or \f$\phi_6^{\rm CM}\f$ the rotation angle of the dilepton system in the centre-of-mass system
 * - 7 = \f$x_q\f$, \f$w_X\f$ mappings, as used in the single- and double-dissociative cases only
 */
class LPAIR final : public cepgen::proc::Process {
public:
  explicit LPAIR(const ParametersList& params)
      : proc::Process(params), pair_(steer<ParticleProperties>("pair")), symmetrise_(steer<bool>("symmetrise")) {}
  LPAIR(const LPAIR& oth) : proc::Process(oth), pair_(oth.pair_), symmetrise_(oth.symmetrise_) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new LPAIR(*this)); }

  void addEventContent() override {
    proc::Process::setEventContent({{Particle::IncomingBeam1, {PDG::proton}},
                                    {Particle::IncomingBeam2, {PDG::proton}},
                                    {Particle::Parton1, {PDG::photon}},
                                    {Particle::Parton2, {PDG::photon}},
                                    {Particle::OutgoingBeam1, {PDG::proton}},
                                    {Particle::OutgoingBeam2, {PDG::proton}},
                                    {Particle::CentralSystem, {pair_.pdgid, pair_.pdgid}}});
  }

  double computeWeight() override;

  void prepareKinematics() override {
    ml2_ = pair_.mass * pair_.mass;
    charge_factor_ = std::pow(pair_.charge / 3., 2);
    beams_mode_ = kinematics().incomingBeams().mode();
    ep1_ = pA().energy();
    ep2_ = pB().energy();
    w12_ = mA2() - mB2();  // mass difference between the two incoming particles
    ss_ = s() + w12_;
    if (const auto rl1 = ss_ * ss_ - 4. * mA2() * s(); utils::positive(rl1))
      sl1_ = std::sqrt(rl1);
    else
      throw CG_FATAL("LPAIR:prepareKinematics") << "Invalid rl1 = " << rl1 << ".";
    p_cm_ = 0.5 * sl1_ * inverseSqrtS();
    mom_prefactor_ = 2. / sl1_;
    p12_ = 0.5 * (s() - mA2() - mB2());
    e1mp1_ = mA2() / (ep1_ + p_cm_);

    CG_DEBUG_LOOP("LPAIR:repareKinematics") << std::scientific << "w12 = " << w12_ << std::fixed
                                            << " -> incoming particles' energy = " << ep1_ << ", " << ep2_ << ".";

    formfac_ = FormFactorsFactory::get().build(kinematics().incomingBeams().formFactors());
    strfun_ = StructureFunctionsFactory::get().build(kinematics().incomingBeams().structureFunctions());

    //--- first define the squared mass range for the diphoton/dilepton system
    const auto w_limits = kinematics()
                              .cuts()
                              .central.mass_sum.compute([](double ext) { return std::pow(ext, 2); })
                              .truncate(Limits{4. * ml2_, s()});
    CG_DEBUG_LOOP("LPAIR:prepareKinematics") << "w limits = " << w_limits << "\n\t"
                                             << "wmax/wmin = " << w_limits.max() / w_limits.min();

    //--- variables mapping
    defineVariable(m_u_t1_, Mapping::linear, {0., 1.}, "u_t1");
    defineVariable(m_u_t2_, Mapping::linear, {0., 1.}, "u_t2");
    defineVariable(m_u_s2_, Mapping::linear, {0., 1.}, "u_s2");
    defineVariable(m_w4_, Mapping::power_law, w_limits, "w4");
    defineVariable(m_theta4_, Mapping::linear, {0., M_PI}, "theta4");
    defineVariable(m_phi6_cm_, Mapping::linear, {0., 2. * M_PI}, "phi6cm");
    defineVariable(m_x6_, Mapping::linear, {-1., 1.}, "x6");

    const auto mx_range = [&](double m_in) {
      return kinematics()
          .cuts()
          .remnants.mx.truncate(Limits{mp_ + PDG::get().mass(PDG::piPlus), sqrtS() - m_in - 2. * pair_.mass})
          .compute([](double m) { return m * m; });
    };
    if (!kinematics().incomingBeams().positive().elastic())  // first outgoing beam particle or remnant mass
      defineVariable(mX2(), Mapping::power_law, mx_range(mA()), "MX2");
    if (!kinematics().incomingBeams().negative().elastic())  // second outgoing beam particle or remnant mass
      defineVariable(mY2(), Mapping::power_law, mx_range(mB()), "MY2");
  }

  void fillKinematics() override {
    // boost of the incoming beams
    pA() = Momentum(0., 0., +p_cm_, ep1_).betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
    pB() = Momentum(0., 0., -p_cm_, ep2_).betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
    // boost of the outgoing beams
    pX().betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
    pY().betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
    // incoming partons
    q1() = pA() - pX();
    q2() = pB() - pY();

    // randomly rotate all particles
    const short rany = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;
    const double ranphi = rnd_gen_->uniform(0., 2. * M_PI);
    const bool mirror = rnd_gen_->uniformInt(0, 1) == 1;
    for (auto* mom : {&q1(), &q2(), &pc(0), &pc(1), &pX(), &pY()}) {
      mom->rotatePhi(ranphi, rany);
      if (symmetrise_ && mirror)
        mom->mirrorZ();
    }
    CG_DEBUG_LOOP("LPAIR:gmufil") << "boosted+rotated PX=" << pX() << "\n\t"
                                  << "boosted+rotated PY=" << pY() << "\n\t"
                                  << "boosted+rotated P(l1)=" << pc(0) << "\n\t"
                                  << "boosted+rotated P(l2)=" << pc(1);

    // first outgoing beam
    auto& op1 = event().oneWithRole(Particle::OutgoingBeam1);
    if (kinematics().incomingBeams().positive().elastic())
      op1.setStatus(Particle::Status::FinalState);  // stable proton
    else {
      op1.setStatus(Particle::Status::Unfragmented);  // fragmenting remnants
      pX().setMass(mX());
    }

    // second outgoing beam
    auto& op2 = event().oneWithRole(Particle::OutgoingBeam2);
    if (kinematics().incomingBeams().negative().elastic())
      op2.setStatus(Particle::Status::FinalState);  // stable proton
    else {
      op2.setStatus(Particle::Status::Unfragmented);  // fragmenting remnants
      pY().setMass(mY());
    }

    // central system
    const short ransign = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;
    event()[Particle::CentralSystem][0].get().setChargeSign(+ransign).setStatus(Particle::Status::FinalState);
    event()[Particle::CentralSystem][1].get().setChargeSign(-ransign).setStatus(Particle::Status::FinalState);
  }

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("γγ → l⁺l¯ (LPAIR)");
    desc.add<int>("nopt", 0).setDescription("Optimised mode? (inherited from LPAIR, by default disabled = 0)");
    desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("Lepton pair considered");
    desc.add<bool>("symmetrise", false).setDescription("Symmetrise along z the central system?");
    return desc;
  }

private:
  static constexpr double constb_ = 0.5 * M_1_PI * M_1_PI * M_1_PI;
  /**
   * Calculate energies and momenta of the
   *  1st, 2nd (resp. the "proton-like" and the "electron-like" incoming particles),
   *  3rd (the "proton-like" outgoing particle),
   *  4th (the two-photons central system), and
   *  5th (the "electron-like" outgoing particle) particles in the overall centre-of-mass frame.
   * \brief Energies/momenta computation for the various particles, in the CM system
   * \return Success state of the operation
   */
  bool orient();
  /**
   * Compute the expression of the matrix element squared for the \f$\gamma\gamma\rightarrow\ell^{+}\ell^{-}\f$ process.
   * It returns the value of the convolution of the form factor or structure functions with the central two-photons matrix element squared.
   * \brief Computes the matrix element squared for the requested process
   * \return Full matrix element for the two-photon production of a pair of spin\f$-\frac{1}{2}-\f$point particles.
   *  It is noted as \f[
   *  M = \frac{1}{4bt_1 t_2}\sum_{i=1}^2\sum_{j=1}^2 u_i v_j t_{ij} = \frac{1}{4}\frac{u_1 v_1 t_{11}+u_2 v_1 t_{21}+u_1 v_2 t_{12}+u_2 v_2 t_{22}}{t_1 t_2 b}
   * \f] where \f$b\f$ = \a bb_ is defined in \a ComputeWeight as : \f[
   *  b = t_1 t_2+\left(w_{\gamma\gamma}\sin^2{\theta^{\rm CM}_6}+4m_\ell\cos^2{\theta^{\rm CM}_6}\right) p_g^2
   * \f]
   */
  double periPP() const;
  /**
   * Describe the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in terms of Lorentz-invariant variables.
   * These variables (along with others) will then be fed into the \a PeriPP method (thus are essential for the evaluation of the full matrix element).
   * \return Value of the Jacobian after the operation
   */
  double pickin();

  const ParticleProperties pair_;
  const bool symmetrise_;

  double ml2_{0.};  ///< squared mass of the outgoing leptons
  double charge_factor_{0.};
  mode::Kinematics beams_mode_;
  double ep1_{0.};  ///< energy of the first proton-like incoming particle
  double ep2_{0.};  ///< energy of the second proton-like incoming particle
  double w12_{0.};  ///< \f$\delta_2=m_1^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
  double ss_{0.};
  double p12_{0.};  ///< \f$p_{12} = \frac{1}{2}\left(s-m_{p_1}^2-m_{p_2}^2\right)\f$
  double e1mp1_{0.};

  // mapped variables
  double m_u_t1_{0.};
  double m_u_t2_{0.};
  double m_u_s2_{0.};
  double m_w4_{0.};       ///< squared mass of the two-photon system
  double m_theta4_{0.};   ///< polar angle of the two-photon system
  double m_phi6_cm_{0.};  ///< azimutal angle of the first outgoing lepton
  double m_x6_{0.};

  double w31_{0.};  ///< \f$\delta_1=m_3^2-m_1^2\f$ as defined in \cite Vermaseren:1982cz
  double w52_{0.};  ///< \f$\delta_4=m_5^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
  double p_cm_{0.}, mom_prefactor_{0.};

  //-- two-photon system
  double ec4_{0.};         ///< energy of the two-photon system
  double pc4_{0.};         ///< 3-momentum norm of the two-photon system
  double pt4_{0.};         ///< transverse momentum of the two-photon system
  double mc4_{0.};         ///< mass of the two-photon system
  double cos_theta4_{0.};  ///< cosine of the polar angle for the two-photon system
  double sin_theta4_{0.};  ///< sine of the polar angle for the two-photon system

  double p1k2_{0.}, p2k1_{0.};

  double q1dq_{0.}, q1dq2_{0.};

  double s1_{0.}, s2_{0.};
  double sa1_{0.}, sa2_{0.}, sl1_{0.};

  double epsilon_{0.};
  double alpha4_{0.}, beta4_{0.}, gamma4_{0.};
  double alpha5_{0.}, gamma5_{0.}, alpha6_{0.}, gamma6_{0.};
  double bb_{0.};

  double gram_{0.};
  /// Deltas such as \f$\delta_5=m_4^2-t_1\f$ as defined in Vermaseren's paper
  /// \cite Vermaseren:1982cz for the full definition of these quantities
  std::array<double, 5> deltas_;
  /**
   * Invariant used to tame divergences in the matrix element computation. It is defined as
   * \f[\Delta = \left(p_1\cdot p_2\right)\left(q_1\cdot q_2\right)-\left(p_1\cdot q_2\right)\left(p_2\cdot q_1\right)\f]
   * with \f$p_i, q_i\f$ the 4-momenta associated to the incoming proton-like particle and to the photon emitted from it.
   */
  double delta_{0.}, delta3_{0.}, delta5_{0.};
  struct {
    double gamma{0.}, beta_gamma{0.};
  } boost_props_;

  std::unique_ptr<formfac::Parameterisation> formfac_;
  std::unique_ptr<strfun::Parameterisation> strfun_;
};

//---------------------------------------------------------------------------------------------

double LPAIR::pickin() {
  const auto map_expo = [](double expo, const Limits& lim, const std::string& var_name = "") {
    /**
     * Define modified variables of integration to avoid peaks integrations (see \cite Vermaseren:1982cz for details)
     * Return a set of two modified variables of integration to maintain the stability of the integrand. These two new variables are :
     * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
     * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
     * \param[in] expo Exponent
     * \param[in] lim Min/maximal value of the variable
     * \param[in] var_name The variable name
     * \return A pair containing the value and the bin width the new variable definition
     */
    const double y = lim.max() / lim.min(), out = lim.min() * std::pow(y, expo), dout = out * std::log(y);
    CG_DEBUG_LOOP("LPAIR:map") << "Mapping variable \"" << var_name << "\" in range (" << lim << ")"
                               << " (max/min = " << y << ")\n\t"
                               << "exponent = " << expo << " => "
                               << "x = " << out << ", dx = " << dout;
    return std::make_pair(out, dout);
  };
  const auto [s2_val, s2_width] =
      map_expo(m_u_s2_, Limits{mc4_ + mY(), sqrtS() - mX()}.compute([](double lim) { return lim * lim; }), "s2");
  s2_ = s2_val;
  if (s2_width <= 0.)
    return 0.;

  const double sp = s() + mX2() - s2_, d3 = s2_ - mB2();
  const double rl2 = sp * sp - 4. * s() * mX2();  // lambda(s, m3**2, sigma)
  if (!utils::positive(rl2)) {
    CG_DEBUG_LOOP("LPAIR:pickin") << "Invalid rl2 = " << rl2 << ".";
    return 0.;
  }
  // definition from eq. (A.4) and (A.5) in [1]
  const auto t1_max = mA2() + mX2() - 0.5 * (ss_ * sp + sl1_ * std::sqrt(rl2)) / s(),
             t1_min = (w31_ * d3 + (d3 - w31_) * (d3 * mA2() - w31_ * mB2()) / s()) / t1_max;
  const auto [t1_val, t1_width] =
      map_expo(m_u_t1_, Limits{t1_min, t1_max}, "t1");  // definition of the first photon propagator (t1 < 0)
  t1() = t1_val;
  if (t1_width >= 0.)
    return 0.;

  const auto r1 = s2_ - t1() + mB2(), r2 = s2_ - m_w4_ + mY2(),
             rl4 = (r1 * r1 - 4. * s2_ * mB2()) * (r2 * r2 - 4. * s2_ * mY2());
  if (!utils::positive(rl4)) {
    CG_DEBUG_LOOP("LPAIR:pickin") << "Invalid rl4 = " << rl4 << ".";
    return 0.;
  }

  const auto d4 = m_w4_ - t1();
  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  const auto t2_max = mB2() + mY2() - 0.5 * (r1 * r2 + std::sqrt(rl4)) / s2_,
             t2_min = (w52_ * d4 + (d4 - w52_) * (d4 * mB2() - w52_ * t1()) / s2_) / t2_max;
  const auto [t2_val, t2_width] =
      map_expo(m_u_t2_, Limits{t2_min, t2_max}, "t2");  // definition of the second photon propagator (t2 < 0)
  t2() = t2_val;
  if (t2_width >= 0.)
    return 0.;

  const auto r3 = m_w4_ - t1() - t2();
  if (gamma4_ = t1() * t2() - 0.25 * r3 * r3; gamma4_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "gamma4 = " << gamma4_ << " >= 0";
    return 0.;
  }

  const auto compute_deltas =
      [this](double var, short sign, double t_1, double mi2_1, double mf2_1, double t_2, double mi2_2, double mf2_2)
      -> std::tuple<double, double, double, double> {
    const auto del1 = t_1 - mi2_2, del2 = t_1 - mi2_1 - mf2_1, del3 = m_w4_ - mf2_2;
    const auto m2diff = mf2_1 - mi2_1;
    const auto compute_sa = [](double t, double mi2, double mf2) {
      return mi2 * t - 0.25 * std::pow(mf2 - mi2 - t, 2);
    };
    const auto sa_1 = compute_sa(t_1, mi2_1, mf2_1), sa_2 = compute_sa(t_2, mi2_2, mf2_2);
    const auto compute_boundaries = [](double sb, double sd, double se) {
      std::pair<double, double> out;
      if (std::fabs((sb - sd) / sd) >= 1.) {
        out.first = sb - sd;
        out.second = se / out.first;
      } else {
        out.second = sb + sd;
        out.first = se / out.second;
      }
      return out;
    };
    double var_pm = 0., var_mp = 0., var_min = 0., var_max = 0.;
    if (mi2_1 == 0.)
      var_max =
          (s() * (t_1 * (s() + del1 - mf2_1) - mi2_2 * mf2_1) + mi2_2 * mf2_1 * (mf2_1 - del1)) / ((s() + w12_) * del2);
    else {
      const auto inv_w1 = 1. / mi2_1;
      const auto sb = mf2_1 + 0.5 * (s() * (t_1 - m2diff) + w12_ * del2) * inv_w1,
                 sd = sl1_ * std::sqrt(-sa_1) * inv_w1,
                 se = (s() * (t_1 * (s() + del2 - mi2_2) - mi2_2 * m2diff) + mf2_1 * (mi2_2 * mf2_1 + w12_ * del1)) *
                      inv_w1;
      std::tie(var_pm, var_max) = compute_boundaries(sb, sd, se);
    }
    {
      const auto inv_t = 1. / t_2;
      const auto sb = mi2_2 + t_1 - 0.5 * (m_w4_ - t_1 - t_2) * (mf2_2 - mi2_2 - t_2) * inv_t,
                 sd = 2. * sign * std::sqrt(sa_2 * gamma4_) * inv_t,
                 se = del3 * del1 + (del3 - del1) * (del3 * mi2_2 - del1 * mf2_2) * inv_t;
      std::tie(var_mp, var_min) = compute_boundaries(sb, sd, se);
    }
    return std::make_tuple(-0.25 * (var_max - var) * (mi2_1 != 0. ? (var_pm - var) * mi2_1 : ss_ * del2),
                           -0.25 * (var_min - var) * (var_mp - var) * t_2,
                           sa_1,
                           0.5 * (var - t_1 - mi2_2));
  };

  if (std::tie(deltas_[0], deltas_[1], sa1_, p2k1_) = compute_deltas(s2_, -1, t1(), mA2(), mX2(), t2(), mB2(), mY2());
      sa1_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "sa1_ = " << sa1_ << " >= 0";
    return 0.;
  }

  const auto dd = deltas_[0] * deltas_[1];
  if (!utils::positive(dd)) {
    CG_WARNING("LPAIR:pickin") << "Invalid dd = " << dd << ".";
    return 0.;
  }

  const auto ap = s2_ * t1() - 0.25 * std::pow(s2_ + t1() - mB2(), 2);
  if (utils::positive(ap)) {
    CG_WARNING("LPAIR:pickin") << "ap = " << ap << " should be strictly negative.";
    return 0.;
  }
  const auto inv_ap = 1. / ap;
  const auto st = s2_ - t1() - mB2();
  delta_ = 0.5 *
           ((mB2() * r3 + 0.5 * (w52_ - t2()) * st) * (p12_ * t1() - 0.25 * (t1() - w31_) * st) -
            std::cos(m_theta4_) * st * std::sqrt(dd)) *
           inv_ap;

  s1_ = t2() + mA2() + 2. * (p12_ * r3 - 2. * delta_) / st;

  const auto jacobian = s2_width * t1_width * t2_width * 0.125 * 0.5 / (sl1_ * std::sqrt(-ap));
  if (!utils::positive(jacobian)) {
    CG_WARNING("LPAIR:pickin") << "Null Jacobian.\n\t"
                               << "ds2=" << s2_width << ", dt1=" << t1_width << ", dt2=" << t2_width << ".";
    return 0.;
  }
  CG_DEBUG_LOOP("LPAIR:pickin") << "ds2=" << s2_width << ", dt1=" << t1_width << ", dt2=" << t2_width << "\n\t"
                                << "Jacobian=" << std::scientific << jacobian << std::fixed;

  gram_ = std::pow(std::sin(m_theta4_), 2) * dd * inv_ap;

  if (std::tie(deltas_[2], deltas_[3], sa2_, p1k2_) = compute_deltas(s1_, +1, t2(), mB2(), mY2(), t1(), mA2(), mX2());
      sa2_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "sa2_ = " << sa2_ << " >= 0";
    return 0.;
  }
  CG_DEBUG_LOOP("LPAIR:pickin") << std::scientific << "deltas = " << deltas_ << std::fixed;

  if (deltas_[4] = deltas_[0] + deltas_[2] +
                   ((p12_ * (t1() - w31_) * 0.5 - mA2() * p2k1_) * (p2k1_ * (t2() - w52_) - mB2() * r3) -
                    delta_ * (2. * p12_ * p2k1_ - mB2() * (t1() - w31_))) /
                       p2k1_;
      !utils::positive(deltas_[4])) {
    CG_WARNING("LPAIR:pickin") << "Invalid dd5=" << deltas_[4] << ", with all deltas=" << deltas_ << ".";
    return 0.;
  }

  return jacobian;
}

//---------------------------------------------------------------------------------------------

bool LPAIR::orient() {
  const auto re = 0.5 * inverseSqrtS();
  delta3_ = re * (s2_ - mX2() + w12_);
  delta5_ = re * (s1_ - mY2() - w12_);

  //----- central two-photon/lepton system
  if (ec4_ = delta3_ + delta5_; ec4_ < mc4_) {
    CG_WARNING("LPAIR:orient") << "ec4_ = " << ec4_ << " < mc4_ = " << mc4_ << "\n\t"
                               << "==> delta3 = " << delta3_ << ", delta5 = " << delta5_;
    return false;
  }
  if (pc4_ = utils::fastSqrtSqDiff(ec4_, mc4_); pc4_ == 0.) {  // protons' momenta are not along the z-axis
    CG_WARNING("LPAIR:orient") << "pzc4 is null and should not be...";
    return false;
  }

  CG_DEBUG_LOOP("LPAIR:orient") << "Central system's energy: E4 = " << ec4_ << "\n\t"
                                << "               momentum: p4 = " << pc4_ << "\n\t"
                                << "         invariant mass: m4 = " << mc4_ << ".";

  pt4_ = mom_prefactor_ * std::sqrt(deltas_[4]);
  if (sin_theta4_ = pt4_ / pc4_; !Limits{-1., 1.}.contains(sin_theta4_)) {
    CG_WARNING("LPAIR:orient") << "Invalid sin(theta4): " << sin_theta4_ << ".";
    return false;
  }
  const auto p14 = +0.5 * (s1_ + t1() - t2() - mX2());
  cos_theta4_ = std::sqrt(1. - sin_theta4_ * sin_theta4_) * (ep1_ * ec4_ < p14 ? -1. : 1.);
  alpha4_ = 1. - cos_theta4_;
  beta4_ = 1. + cos_theta4_;
  if (cos_theta4_ < 0.)
    beta4_ = sin_theta4_ * sin_theta4_ / alpha4_;
  else
    alpha4_ = sin_theta4_ * sin_theta4_ / beta4_;

  CG_DEBUG_LOOP("LPAIR:orient") << "cos(theta4) = " << cos_theta4_ << "\t"
                                << "sin(theta4) = " << sin_theta4_ << "\n\t"
                                << "alpha4 = " << alpha4_ << ", beta4 = " << beta4_;

  //----- outgoing beam states
  const auto rr = mom_prefactor_ * std::sqrt(-gram_) / pt4_;

  //--- beam 1 -> 3
  const auto ep3 = ep1_ - delta3_, pp3 = std::sqrt(ep3 * ep3 - mX2()), pt3 = mom_prefactor_ * std::sqrt(deltas_[0]);
  if (pt3 > pp3) {
    CG_WARNING("LPAIR:orient") << "Invalid momentum for outgoing beam 1.";
    return false;
  }
  if (pt3 < rr) {
    CG_WARNING("LPAIR:orient") << "Invalid momentum balance for outgoing beam 1.";
    return false;
  }
  pX() = Momentum::fromPThetaPhiE(pp3, -std::asin(pt3 / pp3), std::asin(-rr / pt3), ep3);
  CG_DEBUG_LOOP("LPAIR:orient") << "Positive-z beam state:\n\t" << std::scientific << "energy: E3 = " << ep3
                                << ", pt3 = " << pt3 << "\n\t"
                                << "momentum = " << pX() << ".";

  //--- beam 2 -> 5
  const auto ep5 = ep2_ - delta5_, pp5 = std::sqrt(ep5 * ep5 - mY2()), pt5 = mom_prefactor_ * std::sqrt(deltas_[2]);
  if (pt5 > pp5) {
    CG_WARNING("LPAIR:orient") << "Invalid momentum for outgoing beam 2.";
    return false;
  }
  if (pt5 < rr) {
    CG_WARNING("LPAIR:orient") << "Invalid momentum balance for outgoing beam 2.";
    return false;
  }
  pY() = Momentum::fromPThetaPhiE(pp5, M_PI + std::asin(pt5 / pp5), std::asin(+rr / pt5), ep5);
  CG_DEBUG_LOOP("LPAIR:orient") << "Negative-z beam state:\n\t" << std::scientific << "energy: E5 = " << ep5
                                << ", pt5 = " << pt5 << "\n\t"
                                << "momentum = " << pY() << ".";

  // x-axis mirroring
  if (const double a1 = pX().px() - pY().px();
      std::fabs(pt4_ + pX().px() + pY().px()) >= std::fabs(std::fabs(a1) - pt4_)) {
    CG_DEBUG_LOOP("LPAIR:orient") << "|pt4+pt3*cos(phi3)+pt5*cos(phi5)| < | |a1|-pt4 |\n\t"
                                  << "pt4 = " << pt4_ << ".";
    if (a1 < 0.)
      pY().mirrorX();
    else
      pX().mirrorX();
  }
  return true;
}

//---------------------------------------------------------------------------------------------

double LPAIR::computeWeight() {
  w31_ = mX2() - mA2();     // mass difference between the first outgoing particle and the first incoming particle
  w52_ = mY2() - mB2();     // mass difference between the second outgoing particle and the second incoming particle
  mc4_ = std::sqrt(m_w4_);  // compute the two-photon energy for this point

  CG_DEBUG_LOOP("LPAIR:weight") << "Masses dump:\n\t"
                                << "m1 = " << mA() << ", m2 = " << mB() << ", m3 = " << mX() << ", m4 = " << mc4_
                                << ", m5 = " << mY() << ".\n\t"
                                << "w1 = " << mA2() << ", w2 = " << mB2() << ", w3 = " << mX2() << ", w4 = " << m_w4_
                                << ", w5 = " << mY2() << ".\n\t"
                                << "w31 = " << w31_ << ", w52 = " << w52_ << ", w12 = " << w12_ << ".";

  auto jacobian = pickin();
  if (!utils::positive(jacobian)) {
    CG_DEBUG_LOOP("LPAIR:weight") << "Pickin failed.";
    return 0.;
  }
  if (!orient()) {
    CG_DEBUG_LOOP("LPAIR:weight") << "Orient failed.";
    return 0.;
  }

  const double ecm6 = m_w4_ / (2. * mc4_), pp6cm = std::sqrt(ecm6 * ecm6 - ml2_);

  jacobian *= pp6cm / mc4_;

  // Let the most obscure part of this code begin...

  const double e3mp3 = mX2() / (pX().energy() + pX().p());
  const double theta_x = pX().theta(), al3 = std::pow(std::sin(theta_x), 2) / (1. + theta_x);

  // 2-photon system kinematics ?!
  const double eg = (m_w4_ + t1() - t2()) / (2. * mc4_);

  const double gamma4 = ec4_ / mc4_;
  const Momentum pg(-pX().px() * cos_theta4_ - (pX().p() * al3 + e3mp3 - e1mp1_ + delta3_) * sin_theta4_,
                    -pX().py(),
                    -gamma4 * pX().px() * sin_theta4_ + (pX().p() * al3 + e3mp3 - e1mp1_) * gamma4 * cos_theta4_ +
                        mc4_ * delta3_ / (ec4_ + pc4_) - gamma4 * delta3_ * alpha4_);

  CG_DEBUG_LOOP("LPAIR") << "pg = " << pg;

  const auto pt_gam = pg.pt(), p_gam = std::max(std::sqrt(eg * eg - t1()), pg.p() > 0.9 * pt_gam ? pg.p() : -999.);
  const auto cos_phi_gam = pg.px() / pt_gam, sin_phi_gam = pg.py() / pt_gam, sin_theta_gam = pt_gam / p_gam;
  const short theta_sign = pg.pz() > 0. ? 1 : -1;
  const auto cos_theta_gam = theta_sign * std::sqrt(1. - sin_theta_gam * sin_theta_gam);

  const double amap = 0.5 * (m_w4_ - t1() - t2()),
               bmap = 0.5 * std::sqrt((std::pow(m_w4_ - t1() - t2(), 2) - 4. * t1() * t2()) * (1. - 4. * ml2_ / m_w4_)),
               ymap = (amap + bmap) / (amap - bmap), beta = std::pow(ymap, m_x6_);

  // 3D rotation of the first outgoing lepton wrt the CM system
  const auto cos_theta6cm = Limits{-1., 1.}.trim(amap / bmap * (beta - 1.) / (beta + 1.)),
             cos2_theta6cm = cos_theta6cm * cos_theta6cm, sin2_theta6cm = 1. - cos2_theta6cm,
             theta6cm = std::acos(cos_theta6cm);

  // match the Jacobian
  jacobian *= (amap + bmap * cos_theta6cm);
  jacobian *= (amap - bmap * cos_theta6cm);
  jacobian *= 0.5 * std::log(ymap) / amap / bmap;
  if (symmetrise_ &&
      (beams_mode_ == mode::Kinematics::ElasticInelastic || beams_mode_ == mode::Kinematics::InelasticElastic))
    jacobian *= 1.;
  else
    jacobian *= 0.5;

  // First outgoing lepton's 3-momentum in the centre of mass system
  const auto p6cm = Momentum::fromPThetaPhiE(pp6cm, theta6cm, m_phi6_cm_);

  CG_DEBUG_LOOP("LPAIR") << "p3cm6 = " << p6cm;

  const double h1 = p6cm.pz() * sin_theta_gam + p6cm.px() * cos_theta_gam;
  const double pc6z = p6cm.pz() * cos_theta_gam - p6cm.px() * sin_theta_gam;
  const double pc6x = h1 * cos_phi_gam - p6cm.py() * sin_phi_gam;
  const double qcx = 2. * pc6x, qcz = 2. * pc6z;

  const double el6 = (ec4_ * ecm6 + pc4_ * pc6z) / mc4_;
  const double h2 = (ec4_ * pc6z + pc4_ * ecm6) / mc4_;
  CG_DEBUG_LOOP("LPAIR") << "h1 = " << h1 << ", h2 = " << h2;

  // outgoing lepton's kinematics (in the two-photon CM frame)
  const auto pc4 = Momentum::fromPThetaPhiE(pc4_, std::acos(cos_theta4_), 0., ec4_);
  pc(0) = Momentum(+pc6x * cos_theta4_ + h2 * sin_theta4_,
                   p6cm.py() * cos_phi_gam + h1 * sin_phi_gam,
                   -pc6x * sin_theta4_ + h2 * cos_theta4_,
                   el6);
  pc(1) = pc4 - pc(0);
  CG_DEBUG_LOOP("LPAIR") << "Outgoing kinematics\n\t"
                         << " first outgoing lepton: p = " << pc(0).p() << ", E = " << pc(0).energy() << "\n\t"
                         << "second outgoing lepton: p = " << pc(1).p() << ", E = " << pc(1).energy();
  const double phi3 = pX().phi(), cos_phi3 = std::cos(phi3), sin_phi3 = std::sin(phi3);
  const double phi5 = pY().phi(), cos_phi5 = std::cos(phi5), sin_phi5 = std::sin(phi5);

  bb_ = t1() * t2() + (m_w4_ * sin2_theta6cm + 4. * ml2_ * cos2_theta6cm) * p_gam * p_gam;
  q1dq_ = eg * (2. * ecm6 - mc4_) - 2. * p_gam * p6cm.pz();
  q1dq2_ = 0.5 * (m_w4_ - t1() - t2());
  CG_DEBUG_LOOP("LPAIR") << "ecm6 = " << ecm6 << ", mc4 = " << mc4_ << "\n\t"
                         << "eg = " << eg << ", pg = " << p_gam << "\n\t"
                         << "q1dq = " << q1dq_ << ", q1dq2 = " << q1dq2_;

  const double hq = ec4_ * qcz / mc4_;
  const auto qve = Momentum::fromPxPyPzE(+qcx * cos_theta4_ + hq * sin_theta4_,
                                         2. * pc(0).py(),
                                         -qcx * sin_theta4_ + hq * cos_theta4_,
                                         +qcz * pc4_ / mc4_);

  const double c1 = pX().pt() * (qve.px() * sin_phi3 - qve.py() * cos_phi3),
               c2 = pX().pt() * (qve.pz() * ep1_ - qve.energy() * p_cm_),
               c3 = (w31_ * ep1_ * ep1_ + 2. * mA2() * delta3_ * ep1_ - mA2() * delta3_ * delta3_ +
                     pX().pt2() * ep1_ * ep1_) /
                    (pX().pz() * ep1_ + pX().energy() * p_cm_);

  const double b1 = pY().pt() * (qve.px() * sin_phi5 - qve.py() * cos_phi5),
               b2 = pY().pt() * (qve.pz() * ep2_ + qve.energy() * p_cm_),
               b3 = (w52_ * ep2_ * ep2_ + 2. * mB2() * delta5_ * ep2_ - mB2() * delta5_ * delta5_ +
                     pY().pt2() * ep2_ * ep2_) /
                    (pY().pz() * ep2_ - pY().energy() * p_cm_);

  const double r12 = c2 * sin_phi3 + c3 * qve.py(), r13 = -c2 * cos_phi3 - c3 * qve.px();
  const double r22 = b2 * sin_phi5 + b3 * qve.py(), r23 = -b2 * cos_phi5 - b3 * qve.px();

  epsilon_ = p12_ * c1 * b1 + r12 * r22 + r13 * r23;

  gamma5_ = mA2() * c1 * c1 + r12 * r12 + r13 * r13;
  gamma6_ = mB2() * b1 * b1 + r22 * r22 + r23 * r23;

  const double pt3 = pY().pt(), pt5 = pY().pt();
  alpha5_ = -(qve.px() * cos_phi3 + qve.py() * sin_phi3) * pt3 * p1k2_ -
            (ep1_ * qve.energy() - p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
            (delta5_ * qve.pz() + qve.energy() * (p_cm_ + pY().pz())) * c3;
  alpha6_ = -(qve.px() * cos_phi5 + qve.py() * sin_phi5) * pt5 * p2k1_ -
            (ep2_ * qve.energy() + p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
            (delta3_ * qve.pz() - qve.energy() * (p_cm_ - pY().pz())) * b3;

  CG_DEBUG_LOOP("LPAIR") << "alpha5 = " << alpha5_ << "\n\t"
                         << "alpha6 = " << alpha6_;

  ////////////////////////////////////////////////////////////////
  // END of GAMGAMLL subroutine in the FORTRAN version
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // INFO from f.f
  ////////////////////////////////////////////////////////////////

  const Momentum cm = pA() + pB();
  boost_props_.gamma = cm.energy() * inverseSqrtS();
  boost_props_.beta_gamma = cm.pz() * inverseSqrtS();
  CG_DEBUG_LOOP("LPAIR:gmufil") << "sqrt(s)=" << sqrtS() << " GeV, initial two-proton system: " << cm << "\n\t"
                                << "gamma=" << boost_props_.gamma << ", betgam=" << boost_props_.beta_gamma;

  //----- outgoing leptons
  pc(0).betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
  pc(1).betaGammaBoost(boost_props_.gamma, boost_props_.beta_gamma);
  if (!kinematics().cuts().central.contain(event()(Particle::CentralSystem)))
    return 0.;

  const auto peripp = periPP();  // compute the structure functions factors
  if (!utils::positive(peripp))
    return 0.;

  const auto alpha_prod = alphaEM(std::sqrt(-t1())) * alphaEM(std::sqrt(-t2()));
  jacobian *= constb_ * charge_factor_ * alpha_prod * alpha_prod / s();

  CG_DEBUG_LOOP("LPAIR:f") << "Jacobian: " << jacobian << ", str.fun. factor: " << peripp << ".";
  return jacobian * peripp;  // compute the event weight using the Jacobian
}

//---------------------------------------------------------------------------------------------

double LPAIR::periPP() const {
  const double qqq = q1dq_ * q1dq_, qdq = 4. * ml2_ - m_w4_;
  const auto m_em =
      16. *
      Matrix{{(bb_ * (qqq - gamma4_ - qdq * (t1() + t2() + 2. * ml2_)) -
               2. * (t1() + 2. * ml2_) * (t2() + 2. * ml2_) * qqq) *
                  t1() * t2(),
              2. * (-bb_ * (deltas_[1] + gamma6_) - 2. * (t1() + 2. * ml2_) * (sa2_ * qqq + alpha6_ * alpha6_)) * t1()},
             {2. * (-bb_ * (deltas_[3] + gamma5_) - 2. * (t2() + 2. * ml2_) * (sa1_ * qqq + alpha5_ * alpha5_)) * t2(),
              8. * (bb_ * (delta_ * delta_ - gram_) - std::pow(epsilon_ - delta_ * (qdq + q1dq2_), 2) -
                    sa1_ * alpha6_ * alpha6_ - sa2_ * alpha5_ * alpha5_ - sa1_ * sa2_ * qqq)}};

  // compute the electric/magnetic form factors for the two considered parton momenta transfers
  const auto compute_form_factors = [this](bool elastic, double q2, double mi2, double mx2) -> Vector {
    if (elastic) {  // trivial case for elastic photon emission
      const auto ff = (*formfac_)(q2);
      return Vector{ff.FM, ff.FE};
    }
    if (!strfun_)
      throw CG_FATAL("LPAIR:peripp")
          << "Inelastic proton form factors computation requires a structure functions definition!";
    const double xbj = utils::xBj(q2, mi2, mx2);
    if (strfun_->name() == 11 /* SuriYennie */)  // this one requires its own object to deal with FM
      return Vector{strfun_->FM(xbj, q2), strfun_->F2(xbj, q2) * xbj * mp_ / q2};
    return Vector{-2. * strfun_->F1(xbj, q2) / q2, strfun_->F2(xbj, q2) * xbj / q2};
  };
  const double peripp =
      std::pow(t1() * t2() * bb_, -2) *
      (compute_form_factors(kinematics().incomingBeams().positive().elastic(), -t1(), mA2(), mX2()).transposed() *
       m_em * compute_form_factors(kinematics().incomingBeams().negative().elastic(), -t2(), mB2(), mY2()))(0);
  CG_DEBUG_LOOP("LPAIR:peripp") << "bb = " << bb_ << ", qqq = " << qqq << ", qdq = " << qdq << "\n\t"
                                << "e-m matrix = " << m_em << "\n\t"
                                << "=> PeriPP = " << peripp;

  return peripp;
}
REGISTER_PROCESS("lpair", LPAIR);
