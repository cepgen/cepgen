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
      : proc::Process(params),
        opt_(steer<int>("nopt")),
        pair_(steer<ParticleProperties>("pair")),
        symmetrise_(steer<bool>("symmetrise")) {}
  LPAIR(const LPAIR& oth) : proc::Process(oth), opt_(oth.opt_), pair_(oth.pair_), symmetrise_(oth.symmetrise_) {}

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
  void prepareKinematics() override;
  void fillKinematics(bool) override;

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
   * \return Success state of the operation
   */
  bool pickin();

  /// Internal switch for the optimised code version (LPAIR legacy)
  const int opt_;
  const ParticleProperties pair_;
  const bool symmetrise_;

  // mapped variables
  double m_u_t1_{0.};
  double m_u_t2_{0.};
  double m_u_s2_{0.};
  double m_w4_{0.};       ///< squared mass of the two-photon system
  double m_theta4_{0.};   ///< polar angle of the two-photon system
  double m_phi6_cm_{0.};  ///< azimutal angle of the first outgoing lepton
  double m_x6_{0.};

  Limits w_limits_;
  struct Masses {
    double Ml2{0.};  ///< squared mass of the outgoing leptons
    double w12{0.};  ///< \f$\delta_2=m_1^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
    double w31{0.};  ///< \f$\delta_1=m_3^2-m_1^2\f$ as defined in \cite Vermaseren:1982cz
    double w52{0.};  ///< \f$\delta_4=m_5^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
  } masses_;
  double charge_factor_{0.};

  //-- incoming beam particles
  double ep1_{0.};  ///< energy of the first proton-like incoming particle
  double ep2_{0.};  ///< energy of the second proton-like incoming particle
  double p_cm_{0.};

  //-- two-photon system
  double ec4_{0.};         ///< energy of the two-photon system
  double pc4_{0.};         ///< 3-momentum norm of the two-photon system
  double pt4_{0.};         ///< transverse momentum of the two-photon system
  double mc4_{0.};         ///< mass of the two-photon system
  double cos_theta4_{0.};  ///< cosine of the polar angle for the two-photon system
  double sin_theta4_{0.};  ///< sine of the polar angle for the two-photon system

  double p12_{0.};  ///< \f$p_{12} = \frac{1}{2}\left(s-m_{p_1}^2-m_{p_2}^2\right)\f$
  double p1k2_{0.}, p2k1_{0.};
  double p13_{0.};  ///< \f$p_{13} = -\frac{1}{2}\left(t_1-m_{p_1}^2-m_{p_3}^2\right)\f$
  double p14_{0.}, p25_{0.};

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
    double gamma{0.}, betgam{0.};
  } boost_props_;

  double jacobian_{0.};

  /**
   * Define modified variables of integration to avoid peaks integrations (see \cite Vermaseren:1982cz for details)
   * Return a set of two modified variables of integration to maintain the stability of the integrand. These two new variables are :
   * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
   * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
   * \brief Redefine the variables of integration in order to avoid the strong peaking of the integrand
   * \param[in] expo Exponent
   * \param[in] lim Min/maximal value of the variable
   * \param[in] var_name The variable name
   * \return A pair containing the value and the bin width the new variable definition
   * \note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
   *  \f$\mathrm dy_{out}\f$ parameter :
   *  - left unchanged :
   * > `mapw2`, `mapxq`, `mapwx`, `maps2`
   *  - opposite sign :
   * > `mapt1`, `mapt2`
   */
  std::pair<double, double> map(double expo, const Limits& lim, const std::string& var_name = "") {
    const double y = lim.max() / lim.min(), out = lim.min() * std::pow(y, expo), dout = out * log(y);
    CG_DEBUG_LOOP("LPAIR:map") << "Mapping variable \"" << var_name << "\" in range (" << lim << ")"
                               << " (max/min = " << y << ")\n\t"
                               << "exponent = " << expo << " => "
                               << "x = " << out << ", dx = " << dout;
    return {out, dout};
  }
  std::pair<double, double> mapla(double y, double z, int u, const Limits& lim) {
    const double xmb = lim.min() - y - z, xpb = lim.max() - y - z;
    const double c = -4. * y * z;
    const double alp = sqrt(xpb * xpb + c), alm = sqrt(xmb * xmb + c);
    const double am = xmb + alm, ap = xpb + alp;
    const double yy = ap / am, zz = std::pow(yy, u);

    const double out = y + z + 0.5 * (am * zz - c / (am * zz));
    const double ax = std::sqrt(std::pow(out - y - z, 2) + c);
    return {out, ax * log(yy)};
  }
  std::unique_ptr<formfac::Parameterisation> formfac_;
  std::unique_ptr<strfun::Parameterisation> strfun_;
};

//---------------------------------------------------------------------------------------------

void LPAIR::prepareKinematics() {
  masses_.Ml2 = pair_.mass * pair_.mass;
  charge_factor_ = std::pow(pair_.charge / 3., 2);

  formfac_ = FormFactorsFactory::get().build(kinematics().incomingBeams().formFactors());
  strfun_ = StructureFunctionsFactory::get().build(kinematics().incomingBeams().structureFunctions());

  //--- first define the squared mass range for the diphoton/dilepton system
  w_limits_ = kinematics()
                  .cuts()
                  .central.mass_sum.compute([](double ext) { return std::pow(ext, 2); })
                  .truncate(Limits{4. * masses_.Ml2, s()});

  CG_DEBUG_LOOP("LPAIR:prepareKinematics") << "w limits = " << w_limits_ << "\n\t"
                                           << "wmax/wmin = " << w_limits_.max() / w_limits_.min();

  //--- variables mapping
  defineVariable(m_u_t1_, Mapping::linear, {0., 1.}, "u_t1");
  defineVariable(m_u_t2_, Mapping::linear, {0., 1.}, "u_t2");
  defineVariable(m_u_s2_, Mapping::linear, {0., 1.}, "u_s2");
  defineVariable(m_w4_, Mapping::power_law, w_limits_, "w4");
  defineVariable(m_theta4_, Mapping::linear, {0., M_PI}, "theta4");
  defineVariable(m_phi6_cm_, Mapping::linear, {0., 2. * M_PI}, "phi6cm");
  defineVariable(m_x6_, Mapping::linear, {0., 1.}, "x6");

  const double mx0 = mp_ + PDG::get().mass(PDG::piPlus);  // 1.07

  //--- first outgoing beam particle or remnant mass
  if (kinematics().incomingBeams().positive().elastic()) {
    event().oneWithRole(Particle::OutgoingBeam1).setPdgId(event().oneWithRole(Particle::IncomingBeam1).pdgId());
    mX2() = pA().mass2();
  } else {
    const auto wx_lim_ob1 = kinematics()
                                .cuts()
                                .remnants.mx.truncate(Limits{mx0, sqrtS() - mA() - 2. * sqrt(masses_.Ml2)})
                                .compute([](double ext) { return std::pow(ext, 2); });
    defineVariable(mX2(), Mapping::power_law, wx_lim_ob1, "MX2");
  }
  //--- second outgoing beam particle or remnant mass
  if (kinematics().incomingBeams().negative().elastic()) {
    event().oneWithRole(Particle::OutgoingBeam2).setPdgId(event().oneWithRole(Particle::IncomingBeam2).pdgId());
    mY2() = pB().mass2();
  } else {
    const auto wx_lim_ob2 = kinematics()
                                .cuts()
                                .remnants.mx.truncate(Limits{mx0, sqrtS() - mB() - 2. * sqrt(masses_.Ml2)})
                                .compute([](double ext) { return std::pow(ext, 2); });
    defineVariable(mY2(), Mapping::power_law, wx_lim_ob2, "MY2");
  }
}

//---------------------------------------------------------------------------------------------

bool LPAIR::pickin() {
  CG_DEBUG_LOOP("LPAIR") << "Optimised mode? " << opt_;

  jacobian_ = 0.;

  // min(s2) = sigma and sig2 = sigma' in [1]
  const double sig = mc4_ + mY();
  auto s2_range = Limits{sig * sig, s() + mX2() - 2 * mX() * sqrtS()};

  CG_DEBUG_LOOP("LPAIR") << "mc4 = " << mc4_ << "\n\t"
                         << "s2 in range " << s2_range << ".";

  CG_DEBUG_LOOP("LPAIR") << "w1 = " << mA2() << ", w2 = " << mB2() << ", w3 = " << mX2() << ", w4 = " << m_w4_
                         << ", w5 = " << mY2() << ". w31 = " << masses_.w31 << ", w52 = " << masses_.w52
                         << ", w12 = " << masses_.w12 << ".";

  const double ss = s() + masses_.w12, rl1 = ss * ss - 4. * mA2() * s();  // lambda(s, m1**2, m2**2)
  if (!utils::positive(rl1)) {
    CG_DEBUG_LOOP("LPAIR") << "Invalid rl1 = " << rl1 << ".";
    return false;
  }
  sl1_ = std::sqrt(rl1);

  s2_ = 0.;
  double ds2 = 0.;
  if (opt_ == 0) {
    const auto s2 = map(m_u_s2_, s2_range, "s2");
    s2_range.min() = s2_ = s2.first;  // why lower s2 range update?
    ds2 = s2.second;
  }

  CG_DEBUG_LOOP("LPAIR") << "s2 = " << s2_;

  const double sp = s() + mX2() - s2_range.min(), d3 = s2_range.min() - mB2();
  const double rl2 = sp * sp - 4. * s() * mX2();  // lambda(s, m3**2, sigma)
  if (!utils::positive(rl2)) {
    CG_DEBUG_LOOP("LPAIR") << "Invalid rl2 = " << rl2 << ".";
    return false;
  }
  double t1_max = mA2() + mX2() - (ss * sp + sl1_ * std::sqrt(rl2)) / (2. * s());  // definition from eq. (A.4) in [1]
  double t1_min = (masses_.w31 * d3 + (d3 - masses_.w31) * (d3 * mA2() - masses_.w31 * mB2()) / s()) /
                  t1_max;  // definition from eq. (A.5) in [1]
  const auto t1_range = Limits{t1_min, t1_max};

  // ensure the t1 range overlaps with the user-steered Q^2 constraints
  // note: this part was dropped in CDF version
  if (t1_range != t1_range.truncate(-kinematics().cuts().initial.q2))
    return false;

  // definition of first photon propagator
  const auto comp_t1 = map(m_u_t1_, t1_range, "t1");
  t1() = comp_t1.first;
  if (!kinematics().cuts().initial.q2.contains(-t1()))
    return false;
  const double dt1 = -comp_t1.second;  // changes wrt mapt1 : dx->-dx

  CG_DEBUG_LOOP("LPAIR") << "Definition of t1 = " << t1() << " in range " << t1_range << ".";

  deltas_[3] = m_w4_ - t1();

  const double d8 = t1() - mB2();
  const double t13 = t1() - mA2() - mX2();

  sa1_ = -std::pow(t1() - masses_.w31, 2) / 4. + mA2() * t1();
  if (sa1_ >= 0.) {
    CG_WARNING("LPAIR") << "sa1_ = " << sa1_ << " >= 0";
    return false;
  }

  const double sl3 = std::sqrt(-sa1_);

  s2_range.min() = sig * sig;
  double splus;  // computing splus and (s2x=s2max)
  if (mA2() != 0.) {
    const double inv_w1 = 1. / mA2();
    const double sb = mX2() + 0.5 * (s() * (t1() - masses_.w31) + masses_.w12 * t13) * inv_w1;
    const double sd = sl1_ * sl3 * inv_w1;
    const double se =
        (s() * (t1() * (s() + t13 - mB2()) - mB2() * masses_.w31) + mX2() * (masses_.w12 * d8 + mB2() * mX2())) *
        inv_w1;

    if (fabs((sb - sd) / sd) >= 1.) {
      splus = sb - sd;
      s2_range.max() = se / splus;
    } else {
      s2_range.max() = sb + sd;
      splus = se / s2_range.max();
    }
  } else {  // 3
    s2_range.max() =
        (s() * (t1() * (s() + d8 - mX2()) - mB2() * mX2()) + mB2() * mX2() * (mB2() + mX2() - t1())) / (ss * t13);
    splus = s2_range.min();
  }
  // 4
  double s2x = s2_range.max();

  CG_DEBUG_LOOP("LPAIR") << "s2x = s2max = " << s2x;

  if (opt_ < 0) {  // 5
    if (splus > s2_range.min()) {
      s2_range.min() = splus;
      CG_DEBUG_LOOP("LPAIR") << "min(s2) truncated to splus = " << splus;
    }
    const auto s2 = opt_ < -1 ? map(m_u_s2_, s2_range, "s2") : mapla(t1(), mB2(), m_u_s2_, s2_range);  // opt_==-1
    s2_ = s2.first;
    ds2 = s2.second;
    s2x = s2_;
  } else if (opt_ == 0)
    s2x = s2_;  // 6

  CG_DEBUG_LOOP("LPAIR") << "s2x = " << s2x;

  // 7
  const double d6 = m_w4_ - mY2();
  const double r1 = s2x - d8, r2 = s2x - d6;

  const double rl4 = (r1 * r1 - 4. * mB2() * s2x) * (r2 * r2 - 4. * mY2() * s2x);
  if (!utils::positive(rl4)) {
    CG_DEBUG_LOOP("LPAIR") << "Invalid rl4 = " << rl4 << ".";
    return false;
  }
  const double sl4 = std::sqrt(rl4);

  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  const double t2_max = mB2() + mY2() - (r1 * r2 + sl4) / s2x * 0.5,
               t2_min = (masses_.w52 * deltas_[3] +
                         (deltas_[3] - masses_.w52) * (deltas_[3] * mB2() - masses_.w52 * t1()) / s2x) /
                        t2_max;

  // t2, the second photon propagator, is defined here
  const auto comp_t2 = map(m_u_t2_, Limits(t2_min, t2_max), "t2");
  t2() = comp_t2.first;
  if (!kinematics().cuts().initial.q2.contains(-t2()))
    return false;
  const double dt2 = -comp_t2.second;  // changes wrt mapt2 : dx->-dx

  // \f$\delta_6=m_4^2-m_5^2\f$ as defined in Vermaseren's paper
  const double tau = t1() - t2(), r3 = deltas_[3] - t2(), r4 = masses_.w52 - t2();

  CG_DEBUG_LOOP("LPAIR") << "tau = " << tau << ", r1-4 = " << r1 << ", " << r2 << ", " << r3 << ", " << r4;

  const double b = r3 * r4 - 2. * (t1() + mB2()) * t2();
  const double c = t2() * d6 * d8 + (d6 - d8) * (d6 * mB2() - d8 * mY2());

  const double t25 = t2() - mB2() - mY2();

  sa2_ = -0.25 * r4 * r4 + mB2() * t2();
  if (sa2_ >= 0.) {
    CG_WARNING("LPAIR") << "sa2_ = " << sa2_ << " >= 0";
    return false;
  }

  const double sl6 = 2. * std::sqrt(-sa2_);

  gamma4_ = -r3 * r3 / 4. + t1() * t2();
  if (gamma4_ >= 0.) {
    CG_WARNING("LPAIR") << "gamma4 = " << gamma4_ << " >= 0";
    return false;
  }

  const double sl7 = 2. * std::sqrt(-gamma4_);
  const double sl5 = sl6 * sl7;

  double s2p;
  if (fabs((sl5 - b) / sl5) >= 1.) {
    s2p = 0.5 * (sl5 - b) / t2();
    s2_range.min() = c / (t2() * s2p);
  } else {  // 8
    s2_range.min() = 0.5 * (-sl5 - b) / t2();
    s2p = c / (t2() * s2_range.min());
  }
  // 9
  if (opt_ >= 1) {
    const auto s2 = opt_ > 1 ? map(m_u_s2_, s2_range, "s2") : mapla(t1(), mB2(), m_u_s2_, s2_range);
    s2_ = s2.first;
    ds2 = s2.second;
  }

  const double ap = -0.25 * std::pow(s2_ + d8, 2) + s2_ * t1();

  deltas_[0] = 0.25 * (s2_ - s2_range.max()) * (mA2() != 0. ? (splus - s2_) * mA2() : ss * t13);
  deltas_[1] = 0.25 * (s2_ - s2_range.min()) * (s2p - s2_) * t2();

  CG_DEBUG_LOOP("LPAIR") << "\n\t"
                         << "t2       = " << t2() << "\n\t"
                         << "s2       = " << s2_ << "\n\t"
                         << "s2p      = " << s2p << "\n\t"
                         << "splus    = " << splus << "\n\t"
                         << "s2 range = " << s2_range;

  const double yy4 = std::cos(m_theta4_);
  const double dd = deltas_[0] * deltas_[1];
  p12_ = 0.5 * (s() - mA2() - mB2());
  const double st = s2_ - t1() - mB2();
  const double delb = (2. * mB2() * r3 + r4 * st) * (4. * p12_ * t1() - (t1() - masses_.w31) * st) / (16. * ap);

  CG_DEBUG_LOOP("LPAIR") << std::scientific << "dd = " << dd << ", dd1/2 = " << deltas_ << std::fixed;

  if (!utils::positive(dd)) {
    CG_WARNING("LPAIR:pickin") << "Invalid dd = " << dd << ".";
    return false;
  }

  delta_ = delb - 0.5 * yy4 * st * std::sqrt(dd) / ap;
  s1_ = t2() + mA2() + (2. * p12_ * r3 - 4. * delta_) / st;

  if (ap >= 0.) {
    CG_WARNING("LPAIR:pickin") << "ap = " << ap << " >= 0";
    return false;
  }

  jacobian_ = ds2 * dt1 * dt2 * 0.125 * 0.5 / (sl1_ * std::sqrt(-ap));
  if (!utils::positive(jacobian_)) {
    CG_WARNING("LPAIR:pickin") << "Null Jacobian.\n\t"
                               << "D(s2)=" << ds2 << ", D(t1)=" << dt1 << ", D(t2)=" << dt2 << ".";
    return false;
  }

  CG_DEBUG_LOOP("LPAIR:pickin") << "ds2=" << ds2 << ", dt1=" << dt1 << ", dt2=" << dt2 << "\n\t"
                                << "Jacobian=" << std::scientific << jacobian_ << std::fixed;

  gram_ = (1. - yy4 * yy4) * dd / ap;

  p13_ = -0.5 * t13;
  p14_ = 0.5 * (tau + s1_ - mX2());
  p25_ = -0.5 * t25;

  p1k2_ = 0.5 * (s1_ - t2() - mA2());
  p2k1_ = 0.5 * st;

  if (mB2() != 0.) {
    const double inv_w2 = 1. / mB2();
    const double sbb = 0.5 * (s() * (t2() - masses_.w52) - masses_.w12 * t25) * inv_w2 + mY2(),
                 sdd = 0.5 * sl1_ * sl6 * inv_w2,
                 see = (s() * (t2() * (s() + t25 - mA2()) - mA2() * masses_.w52) +
                        mY2() * (mA2() * mY2() - masses_.w12 * (t2() - mA2()))) *
                       inv_w2;
    double s1m, s1p;
    if (sbb * sdd >= 0.) {  // multiplication is more effective than division to check sign+non-null
      s1p = sbb + sdd;
      s1m = see / s1p;
    } else {
      s1m = sbb - sdd;
      s1p = see / s1m;
    }                                                        // 12
    deltas_[2] = -0.25 * mB2() * (s1p - s1_) * (s1m - s1_);  // 13
  } else {                                                   // 14
    const double s1p =
        (s() * (t2() * (s() - mY2() + t2() - mA2()) - mA2() * mY2()) + mA2() * mY2() * (mA2() + mY2() - t2())) /
        (t25 * (s() - masses_.w12));
    deltas_[2] = -0.25 * t25 * (s() - masses_.w12) * (s1p - s1_);
  }
  // 15
  //const double acc3 = (s1p-s1_)/(s1p+s1_);

  const double ssb = t2() + 0.5 * mA2() - r3 * (masses_.w31 - t1()) / t1(), ssd = sl3 * sl7 / t1(),
               sse = (t2() - mA2()) * (m_w4_ - mX2()) +
                     (t2() - m_w4_ + masses_.w31) * ((t2() - mA2()) * mX2() - (m_w4_ - mX2()) * mA2()) / t1();

  double s1pp, s1pm;
  if (ssb / ssd >= 0.) {
    s1pp = ssb + ssd;
    s1pm = sse / s1pp;
  } else {  // 16
    s1pm = ssb - ssd;
    s1pp = sse / s1pm;
  }
  // 17
  deltas_[3] = -0.25 * t1() * (s1_ - s1pp) * (s1_ - s1pm);
  //const double acc4 = ( s1_-s1pm )/( s1_+s1pm );
  deltas_[4] = deltas_[0] + deltas_[2] +
               ((p12_ * (t1() - masses_.w31) * 0.5 - mA2() * p2k1_) * (p2k1_ * (t2() - masses_.w52) - mB2() * r3) -
                delta_ * (2. * p12_ * p2k1_ - mB2() * (t1() - masses_.w31))) /
                   p2k1_;
  if (!utils::positive(deltas_[4])) {
    CG_WARNING("LPAIR") << "Invalid dd5 = " << deltas_[4] << ".";
    return false;
  }

  return true;
}

//---------------------------------------------------------------------------------------------

bool LPAIR::orient() {
  if (!pickin()) {
    CG_DEBUG_LOOP("LPAIR:orient") << "Pickin failed.";
    return false;
  }

  const double re = 0.5 * inverseSqrtS();
  ep1_ = re * (s() + masses_.w12);
  ep2_ = re * (s() - masses_.w12);

  CG_DEBUG_LOOP("LPAIR") << std::scientific << " re = " << re << "\n\t"
                         << "w12 = " << masses_.w12 << std::fixed;
  CG_DEBUG_LOOP("LPAIR") << "Incoming particles' energy = " << ep1_ << ", " << ep2_;

  p_cm_ = re * sl1_;

  delta3_ = re * (s2_ - mX2() + masses_.w12);
  delta5_ = re * (s1_ - mY2() - masses_.w12);

  //----- central two-photon/lepton system

  ec4_ = delta3_ + delta5_;
  if (ec4_ < mc4_) {
    CG_WARNING("LPAIR") << "ec4_ = " << ec4_ << " < mc4_ = " << mc4_ << "\n\t"
                        << "==> delta3 = " << delta3_ << ", delta5 = " << delta5_;
    return false;
  }

  pc4_ = utils::fastSqrtSqDiff(ec4_, mc4_);
  if (pc4_ == 0.) {  // protons' momenta are not along the z-axis
    CG_WARNING("LPAIR") << "pzc4 is null and should not be...";
    return false;
  }

  CG_DEBUG_LOOP("LPAIR") << "Central system's energy: E4 = " << ec4_ << "\n\t"
                         << "               momentum: p4 = " << pc4_ << "\n\t"
                         << "         invariant mass: m4 = " << mc4_ << ".";

  pt4_ = std::sqrt(deltas_[4]) * inverseSqrtS() / p_cm_;
  sin_theta4_ = pt4_ / pc4_;
  if (sin_theta4_ > 1.) {
    CG_WARNING("LPAIR") << "st4 = " << sin_theta4_ << " > 1";
    return false;
  }

  cos_theta4_ = std::sqrt(1. - sin_theta4_ * sin_theta4_);
  if (ep1_ * ec4_ < p14_)
    cos_theta4_ *= -1.;

  alpha4_ = 1. - cos_theta4_;
  beta4_ = 1. + cos_theta4_;

  if (cos_theta4_ < 0.)
    beta4_ = sin_theta4_ * sin_theta4_ / alpha4_;
  else
    alpha4_ = sin_theta4_ * sin_theta4_ / beta4_;

  CG_DEBUG_LOOP("LPAIR") << "cos(theta4) = " << cos_theta4_ << "\t"
                         << "sin(theta4) = " << sin_theta4_ << "\n\t"
                         << "alpha4 = " << alpha4_ << ", beta4 = " << beta4_;

  const double rr = std::sqrt(-gram_) * inverseSqrtS() / (p_cm_ * pt4_);

  //----- outgoing beam states
  const auto prefac = inverseSqrtS() / p_cm_;

  //--- beam 1 -> 3
  const double ep3 = ep1_ - delta3_, pp3 = std::sqrt(ep3 * ep3 - mX2());
  const double pt3 = prefac * std::sqrt(deltas_[0]);

  if (pt3 > pp3) {
    CG_WARNING("LPAIR") << "Invalid momentum for outgoing beam 1.";
    return false;
  }
  if (fabs(rr) > pt3) {
    CG_WARNING("LPAIR") << "Invalid momentum balance for outgoing beam 1.";
    return false;
  }

  pX() = Momentum::fromPThetaPhiE(pp3, -std::asin(pt3 / pp3), std::asin(-rr / pt3), ep3);

  CG_DEBUG_LOOP("LPAIR") << "Positive-z beam state:\n\t" << std::scientific << "energy: E3 = " << ep3
                         << ", pt3 = " << pt3 << "\n\t"
                         << "momentum = " << pX() << ".";

  //--- beam 2 -> 5
  const double ep5 = ep2_ - delta5_, pp5 = std::sqrt(ep5 * ep5 - mY2());
  const double pt5 = prefac * std::sqrt(deltas_[2]);

  if (pt5 > pp5) {
    CG_WARNING("LPAIR") << "Invalid momentum for outgoing beam 2.";
    return false;
  }
  if (fabs(rr) > pt5) {
    CG_WARNING("LPAIR") << "Invalid momentum balance for outgoing beam 2.";
    return false;
  }

  pY() = Momentum::fromPThetaPhiE(pp5, M_PI + std::asin(pt5 / pp5), std::asin(rr / pt5), ep5);

  CG_DEBUG_LOOP("LPAIR") << "Negative-z beam state:\n\t" << std::scientific << "energy: E5 = " << ep5
                         << ", pt5 = " << pt5 << "\n\t"
                         << "momentum = " << pY() << ".";

  //--- mirroring
  const double a1 = pX().px() - pY().px();
  CG_DEBUG_LOOP("LPAIR") << "a1 = " << a1;

  if (std::fabs(pt4_ + pX().px() + pY().px()) < std::fabs(std::fabs(a1) - pt4_)) {
    CG_DEBUG_LOOP("LPAIR") << "|pt4+pt3*cos(phi3)+pt5*cos(phi5)| < | |a1|-pt4 |\n\t"
                           << "pt4 = " << pt4_ << ".";
    return true;
  }
  if (a1 < 0.)
    pY().mirrorX();
  else
    pX().mirrorX();
  return true;
}

//---------------------------------------------------------------------------------------------

double LPAIR::computeWeight() {
  ep1_ = pA().energy();
  ep2_ = pB().energy();
  masses_.w31 = mX2() - mA2();  // mass difference between the first outgoing particle and the first incoming particle
  masses_.w52 = mY2() - mB2();  // mass difference between the second outgoing particle and the second incoming particle
  masses_.w12 = mA2() - mB2();  // mass difference between the two incoming particles
  CG_DEBUG_LOOP("LPAIR") << "sqrt(s) = " << sqrtS() << " GeV\n\t"
                         << "m^2(X) = " << mX2() << " GeV^2, m(X) = " << mX() << " GeV\n\t"
                         << "m^2(Y) = " << mY2() << " GeV^2, m(Y) = " << mY() << " GeV";

  // maximal energy for the CS = beam-beam CM energy - outgoing particles' mass energy
  w_limits_ = w_limits_.truncate(Limits{0., std::pow(sqrtS() - mX() - mY(), 2)});

  mc4_ = std::sqrt(m_w4_);  // compute the two-photon energy for this point
  CG_DEBUG_LOOP("LPAIR") << "Computed value for w4 = " << m_w4_ << " -> mc4 = " << mc4_;

  if (!orient()) {
    CG_DEBUG_LOOP("LPAIR") << "Orient failed.";
    return 0.;
  }
  if (utils::positive(t1()) || utils::positive(t2())) {  // t = -q^2 < 0
    CG_WARNING("LPAIR") << "t1 = " << t1() << ", t2 = " << t2() << ".";
    return 0.;
  }

  const double ecm6 = m_w4_ / (2. * mc4_), pp6cm = std::sqrt(ecm6 * ecm6 - masses_.Ml2);
  const double alpha1 = alphaEM(std::sqrt(-t1())), alpha2 = alphaEM(std::sqrt(-t2()));

  jacobian_ *= pp6cm * constb_ * charge_factor_ * alpha1 * alpha1 * alpha2 * alpha2 / mc4_ / s();

  // Let the most obscure part of this code begin...

  const double e1mp1 = mA2() / (ep1_ + p_cm_);
  const double e3mp3 = mX2() / (pX().energy() + pX().p());

  const double al3 = std::pow(sin(pX().theta()), 2) / (1. + (pX().theta()));

  // 2-photon system kinematics ?!
  const double eg = (m_w4_ + t1() - t2()) / (2. * mc4_);
  double p_gam = std::sqrt(eg * eg - t1());

  const double gamma4 = ec4_ / mc4_;
  const Momentum pg(-pX().px() * cos_theta4_ - (pX().p() * al3 + e3mp3 - e1mp1 + delta3_) * sin_theta4_,
                    -pX().py(),
                    -gamma4 * pX().px() * sin_theta4_ + (pX().p() * al3 + e3mp3 - e1mp1) * gamma4 * cos_theta4_ +
                        mc4_ * delta3_ / (ec4_ + pc4_) - gamma4 * delta3_ * alpha4_);

  CG_DEBUG_LOOP("LPAIR") << "pg = " << pg;

  const double pt_gam = pg.pt(), p_gam_tmp = pg.p();
  if (p_gam_tmp > pt_gam * 0.9 && p_gam_tmp > p_gam)
    p_gam = p_gam_tmp;  //FIXME ???

  // angles for the 2-photon system ?!
  const double cos_phi_gam = pg.px() / pt_gam, sin_phi_gam = pg.py() / pt_gam;
  const double sin_theta_gam = pt_gam / p_gam;

  const int theta_sign = pg.pz() > 0. ? 1 : -1;
  const double cos_theta_gam = theta_sign * std::sqrt(1. - sin_theta_gam * sin_theta_gam);

  const double amap = 0.5 * (m_w4_ - t1() - t2()),
               bmap = 0.5 * std::sqrt((std::pow(m_w4_ - t1() - t2(), 2) - 4. * t1() * t2()) *
                                      (1. - 4. * masses_.Ml2 / m_w4_)),
               ymap = (amap + bmap) / (amap - bmap), beta = std::pow(ymap, 2. * m_x6_ - 1.);
  const auto xx6 = Limits{0., 1.}.trim(0.5 * (1. + amap / bmap * (beta - 1.) / (beta + 1.)));  // in [0., 1.]

  CG_DEBUG_LOOP("LPAIR") << "amap = " << amap << "\n\t"
                         << "bmap = " << bmap << "\n\t"
                         << "ymap = " << ymap << "\n\t"
                         << "beta = " << beta;

  // 3D rotation of the first outgoing lepton wrt the CM system
  const double cos_theta6cm = 1. - 2. * xx6, cos2_theta6cm = cos_theta6cm * cos_theta6cm,
               sin2_theta6cm = 1. - cos2_theta6cm, sin_theta6cm = std::sqrt(sin2_theta6cm),
               theta6cm = std::acos(cos_theta6cm);
  CG_DEBUG_LOOP("LPAIR") << "theta6cm = " << theta6cm << "(cos = " << cos_theta6cm << ", "
                         << "sin = " << sin_theta6cm << ").";

  // match the Jacobian
  jacobian_ *= (amap + bmap * cos_theta6cm);
  jacobian_ *= (amap - bmap * cos_theta6cm);
  jacobian_ /= amap;
  jacobian_ /= bmap;
  jacobian_ *= log(ymap);
  if ((kinematics().incomingBeams().mode() == mode::Kinematics::ElasticInelastic ||
       kinematics().incomingBeams().mode() == mode::Kinematics::InelasticElastic) &&
      symmetrise_)
    jacobian_ *= 1.;
  else
    jacobian_ *= 0.5;
  CG_DEBUG_LOOP("LPAIR") << "Jacobian = " << jacobian_;

  // First outgoing lepton's 3-momentum in the centre of mass system
  auto p6cm = Momentum::fromPThetaPhiE(pp6cm, theta6cm, m_phi6_cm_);

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

  bb_ = t1() * t2() + (m_w4_ * sin2_theta6cm + 4. * masses_.Ml2 * cos2_theta6cm) * p_gam * p_gam;
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
               c3 = (masses_.w31 * ep1_ * ep1_ + 2. * mA2() * delta3_ * ep1_ - mA2() * delta3_ * delta3_ +
                     pX().pt2() * ep1_ * ep1_) /
                    (pX().energy() * p_cm_ + pX().pz() * ep1_);

  const double b1 = pY().pt() * (qve.px() * sin_phi5 - qve.py() * cos_phi5),
               b2 = pY().pt() * (qve.pz() * ep2_ + qve.energy() * p_cm_),
               b3 = (masses_.w52 * ep2_ * ep2_ + 2. * mB2() * delta5_ * ep2_ - mB2() * delta5_ * delta5_ +
                     pY().pt2() * ep2_ * ep2_) /
                    (ep2_ * pY().pz() - pY().energy() * p_cm_);

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
  boost_props_.betgam = cm.pz() * inverseSqrtS();
  CG_DEBUG_LOOP("LPAIR:gmufil") << "sqrt(s)=" << sqrtS() << " GeV, initial two-proton system: " << cm << "\n\t"
                                << "gamma=" << boost_props_.gamma << ", betgam=" << boost_props_.betgam;

  //----- outgoing leptons
  const auto mass_before = (pc(0) + pc(1)).mass();
  pc(0).betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  pc(1).betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  CG_DEBUG_LOOP("LPAIR:gmufil") << "Invariant mass imbalance after beta/gamma boost:"
                                << (pc(0) + pc(1)).mass() - mass_before << ".";
  if (!kinematics().cuts().central.contain(event()(Particle::CentralSystem)))
    return 0.;

  if (const auto peripp = periPP(); peripp > 0.) {  // compute the structure functions factors
    CG_DEBUG_LOOP("LPAIR:f") << "Jacobian: " << jacobian_ << ", str.fun. factor: " << peripp << ".";
    return jacobian_ * peripp;  // compute the event weight using the Jacobian
  }
  return 0.;
}

//---------------------------------------------------------------------------------------------

void LPAIR::fillKinematics(bool) {
  // boost of the incoming beams
  pA() = Momentum(0., 0., +p_cm_, ep1_).betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  pB() = Momentum(0., 0., -p_cm_, ep2_).betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  // boost of the outgoing beams
  pX().betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  pY().betaGammaBoost(boost_props_.gamma, boost_props_.betgam);
  // incoming partons
  q1() = pA() - pX();
  q2() = pB() - pY();

  // randomly rotate all particles
  const short rany = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;
  const double ranphi = rnd_gen_->uniform(0., 2. * M_PI);
  q1().rotatePhi(ranphi, rany);
  q2().rotatePhi(ranphi, rany);
  pc(0).rotatePhi(ranphi, rany);
  pc(1).rotatePhi(ranphi, rany);
  pX().rotatePhi(ranphi, rany);
  pY().rotatePhi(ranphi, rany);
  if (symmetrise_ && rnd_gen_->uniformInt(0, 1) == 1) {
    q1().mirrorZ();
    q2().mirrorZ();
    pc(0).mirrorZ();
    pc(1).mirrorZ();
    pX().mirrorZ();
    pY().mirrorZ();
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
  auto central_system = event()[Particle::CentralSystem];
  const short ransign = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;
  auto& ol1 = central_system[0].get();
  ol1.setChargeSign(+ransign);
  ol1.setStatus(Particle::Status::FinalState);
  auto& ol2 = central_system[1].get();
  ol2.setChargeSign(-ransign);
  ol2.setStatus(Particle::Status::FinalState);
  event().oneWithRole(Particle::Intermediate).setMomentum(pc(0) + pc(1), true);
}

//---------------------------------------------------------------------------------------------

double LPAIR::periPP() const {
  const double qqq = q1dq_ * q1dq_, qdq = 4. * masses_.Ml2 - m_w4_;
  const auto m_em =
      16. *
      Matrix{{(bb_ * (qqq - gamma4_ - qdq * (t1() + t2() + 2. * masses_.Ml2)) -
               2. * (t1() + 2. * masses_.Ml2) * (t2() + 2. * masses_.Ml2) * qqq) *
                  t1() * t2(),
              2. * (-bb_ * (deltas_[1] + gamma6_) - 2. * (t1() + 2. * masses_.Ml2) * (sa2_ * qqq + alpha6_ * alpha6_)) *
                  t1()},
             {2. * (-bb_ * (deltas_[3] + gamma5_) - 2. * (t2() + 2. * masses_.Ml2) * (sa1_ * qqq + alpha5_ * alpha5_)) *
                  t2(),
              8. * (bb_ * (delta_ * delta_ - gram_) - std::pow(epsilon_ - delta_ * (qdq + q1dq2_), 2) -
                    sa1_ * alpha6_ * alpha6_ - sa2_ * alpha5_ * alpha5_ - sa1_ * sa2_ * qqq)}};

  // compute the electric/magnetic form factors for the two considered parton momenta transfers
  const auto compute_form_factors = [this](bool elastic, double q2, double mx2) -> Vector {
    if (elastic) {  // trivial case for elastic photon emission
      const auto ff = (*formfac_)(q2);
      return Vector{ff.FM, ff.FE};
    }
    if (!strfun_)
      throw CG_FATAL("LPAIR:computeFormFactors")
          << "Inelastic proton form factors computation requires a structure functions definition!";
    const double xbj = utils::xBj(q2, mp2_, mx2);
    if (strfun_->name() == 11 /* SuriYennie */)  // this one requires its own object to deal with FM
      return Vector{strfun_->FM(xbj, q2), strfun_->F2(xbj, q2) * xbj * mp_ / q2};
    return Vector{-2. * strfun_->F1(xbj, q2) / q2, strfun_->F2(xbj, q2) * xbj / q2};
  };
  const auto fp1 = compute_form_factors(kinematics().incomingBeams().positive().elastic(), -t1(), mX2()),
             fp2 = compute_form_factors(kinematics().incomingBeams().negative().elastic(), -t2(), mY2());

  const double peripp = std::pow(t1() * t2() * bb_, -2) * (fp1.transposed() * m_em * fp2)(0);
  CG_DEBUG_LOOP("LPAIR:peripp") << "bb = " << bb_ << ", qqq = " << qqq << ", qdq = " << qdq << "\n\t"
                                << "e-m matrix = " << m_em << "\n\t"
                                << "=> PeriPP = " << peripp;

  return peripp;
}
// register process
REGISTER_PROCESS("lpair", LPAIR);
