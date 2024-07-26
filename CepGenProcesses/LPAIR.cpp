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

/// Matrix element for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process as defined in \cite Vermaseren:1982cz
class LPAIR final : public cepgen::proc::Process {
public:
  explicit LPAIR(const ParametersList& params)
      : proc::Process(params),
        pair_(steer<ParticleProperties>("pair")),
        symmetrise_(steer<bool>("symmetrise")),
        randomise_charge_(steer<bool>("randomiseCharge")) {}
  LPAIR(const LPAIR& oth)
      : proc::Process(oth), pair_(oth.pair_), symmetrise_(oth.symmetrise_), randomise_charge_(oth.randomise_charge_) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new LPAIR(*this)); }

  void addEventContent() override {
    proc::Process::setEventContent({{Particle::IncomingBeam1, {PDG::proton}},
                                    {Particle::IncomingBeam2, {PDG::proton}},
                                    {Particle::Parton1, {PDG::photon}},
                                    {Particle::Parton2, {PDG::photon}},
                                    {Particle::OutgoingBeam1, {PDG::proton}},
                                    {Particle::OutgoingBeam2, {PDG::proton}},
                                    {Particle::CentralSystem, {+(spdgid_t)pair_.pdgid, -(spdgid_t)pair_.pdgid}}});
  }

  double computeWeight() override;

  void prepareKinematics() override {
    ml_ = pair_.mass;
    ml2_ = ml_ * ml_;
    charge_factor_ = std::pow(pair_.integerCharge() / 3., 2);
    beams_mode_ = kinematics().incomingBeams().mode();
    pA() = kinematics().incomingBeams().positive().momentum();
    pB() = kinematics().incomingBeams().negative().momentum();
    if (re_ = 0.5 * inverseSqrtS(); !utils::positive(re_))
      throw CG_FATAL("LPAIR:prepareKinematics") << "Invalid centre of mass energy: sqrt(s)=" << sqrtS() << ".";
    w12_ = mA2() - mB2();       // mass difference between the two incoming particles
    ep1_ = re_ * (s() + w12_);  // in centre of mass system (pp != ep)
    ep2_ = re_ * (s() - w12_);
    ss_ = s() + w12_;
    if (const auto rl1 = ss_ * ss_ - 4. * mA2() * s(); rl1 >= 0.)
      sl1_ = std::sqrt(rl1);
    else
      throw CG_FATAL("LPAIR:prepareKinematics") << "Invalid rl1 = " << rl1 << ".";
    p_cm_ = 0.5 * sl1_ * inverseSqrtS();
    mom_prefactor_ = 2. / sl1_;
    p12_ = 0.5 * (s() - mA2() - mB2());
    e1mp1_ = mA2() / (ep1_ + p_cm_);
    {  // definition of boost-to-lab boost variables
      const Momentum cm = pA() + pB();
      gamma_cm_ = cm.energy() * inverseSqrtS();
      beta_gamma_cm_ = cm.pz() * inverseSqrtS();
      CG_DEBUG_LOOP("LPAIR:prepareKinematics")
          << "sqrt(s)=" << sqrtS() << " GeV, initial two-proton system: " << cm << "\n\t"
          << "gamma=" << gamma_cm_ << ", beta*gamma=" << beta_gamma_cm_;
    }

    formfac1_ = FormFactorsFactory::get().build(kinematics().incomingBeams().positive().formFactors());
    formfac2_ = FormFactorsFactory::get().build(kinematics().incomingBeams().negative().formFactors());
    strfun_ = StructureFunctionsFactory::get().build(kinematics().incomingBeams().structureFunctions());
    is_strfun_sy_ = strfun_->name() == "SuriYennie";

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

    mX2() = mA2();
    mY2() = mB2();
    const auto mx_range = [&](double m_in) {
      return kinematics()
          .cuts()
          .remnants.mx.truncate(Limits{mp_ + PDG::get().mass(PDG::piPlus), sqrtS() - m_in - 2. * pair_.mass})
          .compute([](double m) { return m * m; });
    };
    if (beams_mode_ != mode::Kinematics::ElasticElastic)  // first outgoing beam particle or remnant mass
      defineVariable(mX2(), Mapping::power_law, mx_range(mA()), "MX2");
    if (beams_mode_ == mode::Kinematics::InelasticInelastic)  // second outgoing beam particle or remnant mass
      defineVariable(mY2(), Mapping::power_law, mx_range(mB()), "MY2");
    if (symmetrise_ &&
        (beams_mode_ == mode::Kinematics::InelasticElastic || beams_mode_ == mode::Kinematics::ElasticInelastic))
      CG_INFO("LPAIR:prepareKinematics")
          << "Single dissociation kinematics mode was enabled with symmetrisation of the outgoing system.\n\t"
             "The generator-level cross section will be doubled, and beam particles, incoming partons, and central "
             "system will be mirrored in z.";
  }

  void fillKinematics() override {
    // boost of the incoming beams
    pA() = Momentum(0., 0., +p_cm_, ep1_).betaGammaBoost(gamma_cm_, beta_gamma_cm_);
    pB() = Momentum(0., 0., -p_cm_, ep2_).betaGammaBoost(gamma_cm_, beta_gamma_cm_);
    // boost of the outgoing beams
    pX().setMass(mX()).betaGammaBoost(gamma_cm_, beta_gamma_cm_);
    pY().setMass(mY()).betaGammaBoost(gamma_cm_, beta_gamma_cm_);
    // incoming partons
    q1() = pA() - pX();
    q2() = pB() - pY();

    // randomly rotate all particles
    const short rany = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;
    const double ranphi = rnd_gen_->uniform(0., 2. * M_PI);
    for (auto* mom : {&q1(), &q2(), &pX(), &pY(), &pc(0), &pc(1)})
      mom->rotatePhi(ranphi, rany);
    if ((symmetrise_ && rnd_gen_->uniformInt(0, 1) == 1) ||
        beams_mode_ == mode::Kinematics::ElasticInelastic) {  // mirror X/Y and dilepton systems if needed
      std::swap(pX(), pY());
      std::swap(q1(), q2());
      std::swap(pc(0), pc(1));
      for (auto* mom : {&q1(), &q2(), &pX(), &pY(), &pc(0), &pc(1)})
        mom->mirrorZ();
    }
    // first outgoing beam
    event()
        .oneWithRole(Particle::OutgoingBeam1)
        .setStatus(kinematics().incomingBeams().positive().elastic() ? Particle::Status::FinalState
                                                                     : Particle::Status::Unfragmented);
    // second outgoing beam
    event()
        .oneWithRole(Particle::OutgoingBeam2)
        .setStatus(kinematics().incomingBeams().negative().elastic() ? Particle::Status::FinalState
                                                                     : Particle::Status::Unfragmented);

    // central system
    const auto ransign = rnd_gen_->uniformInt(0, 1) == 1;
    if (randomise_charge_) {  // randomise the charge of outgoing system
      event()[Particle::CentralSystem][0].get().setAntiparticle(ransign);
      event()[Particle::CentralSystem][1].get().setAntiparticle(!ransign);
    }
    event()[Particle::CentralSystem][0].get().setStatus(Particle::Status::FinalState);
    event()[Particle::CentralSystem][1].get().setStatus(Particle::Status::FinalState);
  }

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("γγ → l⁺l¯ (LPAIR)");
    desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("Lepton pair considered");
    desc.add<bool>("symmetrise", false).setDescription("Symmetrise along z the central system?");
    desc.add<bool>("randomiseCharge", true).setDescription("randomise the charges of the two central fermions?");
    return desc;
  }

private:
  static constexpr double constb_ = 0.5 * M_1_PI * M_1_PI * M_1_PI;
  /// Calculate energies and momenta of full event content, in the CM system
  bool orient();
  /// Compute the squared matrix element squared for the \f$\gamma\gamma\rightarrow\ell^{+}\ell^{-}\f$ process
  /// \return Convolution of the form factor or structure functions with the squared central two-photons matrix element (for a pair of spin\f$-\frac{1}{2}-\f$point particles)
  /// \note Its expression is of the form:
  ///     \f[M = \frac{1}{4bt_1 t_2}\sum_{i=1}^2\sum_{j=1}^2 u_i v_j t_{ij} = \frac{1}{4}\frac{u_1 v_1 t_{11}+u_2 v_1 t_{21}+u_1 v_2 t_{12}+u_2 v_2 t_{22}}{t_1 t_2 b}\f]
  ///   where \f$b\f$ = \a bb_ is defined in \a computeWeight as:
  ///     \f[b = t_1 t_2+\left(w_{\gamma\gamma}\sin^2{\theta^{\rm CM}_6}+4m_\ell\cos^2{\theta^{\rm CM}_6}\right) p_g^2\f]
  double periPP() const {
    const auto qdq = 4. * ml2_ - m_w4_;
    const auto m_em =
        Matrix{{(bb_ * (q2dq_ - gamma4_ - qdq * (t1() + t2() + 2. * ml2_)) -
                 2. * (t1() + 2. * ml2_) * (t2() + 2. * ml2_) * q2dq_) *
                    t1() * t2(),
                2. * (-bb_ * (deltas1_[1] + gamma6_) - 2. * (t1() + 2. * ml2_) * (sa2_ * q2dq_ + alpha6_ * alpha6_)) *
                    t1()},
               {2. * (-bb_ * (deltas2_[1] + gamma5_) - 2. * (t2() + 2. * ml2_) * (sa1_ * q2dq_ + alpha5_ * alpha5_)) *
                    t2(),
                8. * (bb_ * (delta_ * delta_ - gram_) -
                      std::pow(epsilon_ - delta_ * (qdq + 0.5 * (m_w4_ - t1() - t2())), 2) - sa1_ * alpha6_ * alpha6_ -
                      sa2_ * alpha5_ * alpha5_ - sa1_ * sa2_ * q2dq_)}} *
        std::pow(4. / t1() / t2() / bb_, 2);

    // compute the electric/magnetic form factors for the two considered parton momenta transfers
    const auto compute_form_factors =
        [this](formfac::Parameterisation& formfac, bool elastic, double q2, double mi2, double mx2) -> Vector {
      if (elastic) {  // trivial case for elastic photon emission
        const auto ff = formfac(q2);
        return Vector{ff.FM, ff.FE};
      }
      if (!strfun_)
        throw CG_FATAL("LPAIR:peripp")
            << "Inelastic proton form factors computation requires a structure functions definition!";
      const double xbj = utils::xBj(q2, mi2, mx2);
      if (is_strfun_sy_)  // this one requires its own object to deal with FM
        return Vector{strfun_->FM(xbj, q2), strfun_->F2(xbj, q2) * xbj * mp_ / q2};
      return Vector{-2. * strfun_->F1(xbj, q2) / q2, strfun_->F2(xbj, q2) * xbj / q2};
    };
    const auto u1 = beams_mode_ == mode::Kinematics::ElasticInelastic
                        ? compute_form_factors(*formfac2_, false, -t1(), mA2(), mX2())
                        : compute_form_factors(
                              *formfac1_, kinematics().incomingBeams().positive().elastic(), -t1(), mA2(), mX2()),
               u2 = beams_mode_ == mode::Kinematics::ElasticInelastic
                        ? compute_form_factors(*formfac1_, true, -t2(), mB2(), mY2())
                        : compute_form_factors(
                              *formfac2_, kinematics().incomingBeams().negative().elastic(), -t2(), mB2(), mY2());
    const auto peripp = (u1.transposed() * m_em * u2)(0);
    CG_DEBUG_LOOP("LPAIR:peripp") << "bb = " << bb_ << ", qqq = " << q2dq_ << ", qdq = " << qdq << "\n\t"
                                  << "e-m matrix=\n"
                                  << m_em << "\n\t"
                                  << "u1-2: " << u1 << ", " << u2 << " -> PeriPP = " << peripp << ".";
    return peripp;
  }
  /// Describe the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in terms of Lorentz-invariant variables.
  /// \note These variables (along with others) will then be fed into the \a periPP method (thus are essential for the evaluation of the full matrix element).
  /// \return Value of the Jacobian after the operation
  double pickin();

  const ParticleProperties pair_;
  const bool symmetrise_;
  const bool randomise_charge_;

  ///////////////////////////////////////////////////////////////////
  // variables computed at phase space definition
  double ml_{0.};   ///< mass of the outgoing leptons
  double ml2_{0.};  ///< squared mass of the outgoing leptons
  double charge_factor_{0.};
  mode::Kinematics beams_mode_{mode::Kinematics::invalid};
  double re_{0.};
  double ep1_{0.};  ///< energy of the first proton-like incoming particle
  double ep2_{0.};  ///< energy of the second proton-like incoming particle
  double w12_{0.};  ///< \f$\delta_2=m_1^2-m_2^2\f$ as defined in \cite Vermaseren:1982cz
  double ss_{0.};
  double p12_{0.};  ///< \f$p_{12} = \frac{1}{2}\left(s-m_{p_1}^2-m_{p_2}^2\right)\f$
  double sl1_{0.};
  double e1mp1_{0.};
  double p_cm_{0.}, mom_prefactor_{0.};
  double gamma_cm_{0.}, beta_gamma_cm_{0.};

  std::unique_ptr<formfac::Parameterisation> formfac1_, formfac2_;
  std::unique_ptr<strfun::Parameterisation> strfun_;
  bool is_strfun_sy_{false};

  // mapped variables
  double m_u_t1_{0.};  ///< \f$t_1\f$, first parton normalised virtuality
  double m_u_t2_{0.};  ///< \f$t_2\f$, second parton normalised virtuality
  double m_u_s2_{0.};  ///< \f$s_2\f$
  // yy4 = \f$\cos\left(\pi x_3\right)\f$
  double m_w4_{0.};       ///< \f$w_4\f$, squared invariant mass of the two-parton system
  double m_theta4_{0.};   ///< polar angle of the two-photon system
  double m_phi6_cm_{0.};  ///< \f$\phi_6^{\rm CM}\f$, azimutal angle of the first outgoing lepton
  /// xx6 = \f$\frac{1}{2}\left(1-\cos\theta^{\rm CM}_6\right)\f$ definition (3D rotation of the first outgoing lepton with respect to the two-photon centre-of-mass system).
  /// \note If the \a nm_ optimisation flag is set this angle coefficient value becomes
  ///     \f[\frac{1}{2}\left(\frac{a_{\rm map}}{b_{\rm map}}\frac{\beta-1}{\beta+1}+1\right)\f] with
  ///     \f$a_{\rm map}=\frac{1}{2}\left(w_4-t_1-t_2\right)\f$, \f$b_{\rm map}=\frac{1}{2}\sqrt{\left(\left(w_4-t_1-t_2\right)^2-4t_1t_2\right)\left(1-4\frac{w_6}{w_4}\right)}\f$,
  ///   and
  ///     \f$\beta=\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)^{2x_5-1}\f$
  ///   and the Jacobian element is scaled by a factor
  ///     \f$\frac{1}{2}\frac{\left(a_{\rm map}^2-b_{\rm map}^2\cos^2\theta^{\rm CM}_6\right)}{a_{\rm map}b_{\rm map}}\log\left(\frac{a_{\rm map}+b_{\rm map}}{a_{\rm map}-b_{\rm map}}\right)\f$
  double m_x6_{0.};

  ///////////////////////////////////////////////////////////////////
  // variables computed for each phase space point computation
  double s1_{0.}, s2_{0.};
  double sa1_{0.}, sa2_{0.};
  double p1k2_{0.}, p2k1_{0.};
  double ec4_{0.};         ///< central system energy
  double pc4_{0.};         ///< central system 3-momentum norm
  double pt4_{0.};         ///< central system transverse momentum
  double mc4_{0.};         ///< central system invariant mass
  double cos_theta4_{0.};  ///< central system polar angle cosine
  double sin_theta4_{0.};  ///< central system polar angle sine
  double q2dq_{0.};
  double epsilon_{0.};
  double alpha4_{0.}, beta4_{0.}, gamma4_{0.};
  double alpha5_{0.}, gamma5_{0.}, alpha6_{0.}, gamma6_{0.};
  double bb_{0.};
  double gram_{0.};
  double dd5_{0.};
  std::array<double, 2> deltas1_, deltas2_;
  /// Invariant used to tame divergences in the matrix element computation
  /// \note Defined as \f[\Delta = \left(p_1\cdot p_2\right)\left(q_1\cdot q_2\right)-\left(p_1\cdot q_2\right)\left(p_2\cdot q_1\right)\f]
  ///   with \f$p_i, q_i\f$ the 4-momenta associated to the incoming proton-like particle and to the photon emitted from it.
  double delta_{0.};
  double eph1_{0.}, eph2_{0.};
};

//---------------------------------------------------------------------------------------------

double LPAIR::pickin() {
  /// Define modified variables of integration to avoid peaks integrations
  /// \return two modified variables of integration to maintain the stability of the integrand:
  ///   - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
  ///   - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$,
  ///     the new variable's differential form
  /// \param[in] expo Exponent
  /// \param[in] lim Min/maximal value of the variable
  const auto map_expo = [](double expo, const Limits& lim) {
    const double y = lim.max() / lim.min(), out = lim.min() * std::pow(y, expo);
    return std::make_pair(out, out * std::log(y));
  };
  const auto [s2_val, s2_width] =
      map_expo(m_u_s2_, Limits{mc4_ + mY(), sqrtS() - mX()}.compute([](double lim) { return lim * lim; }));
  s2_ = s2_val;
  if (s2_width <= 0.)
    return 0.;

  const double sp = s() + mX2() - s2_, d3 = s2_ - mB2();
  const double rl2 = sp * sp - 4. * s() * mX2();  // lambda(s, m3**2, sigma)
  if (!utils::positive(rl2)) {
    CG_WARNING("LPAIR:pickin") << "Invalid rl2 = " << rl2 << ".";
    return 0.;
  }
  const auto w31 = mX2() - mA2();
  // definition from eq. (A.4) and (A.5) in [1]
  auto t1_max = mA2() + mX2() - 0.5 * (ss_ * sp + sl1_ * std::sqrt(rl2)) / s(),
       t1_min = (w31 * d3 + (d3 - w31) * (d3 * mA2() - w31 * mB2()) / s()) / t1_max;
  t1_max = std::max(t1_max, -kinematics().cuts().initial.q2.at(0).max());
  t1_min = std::min(t1_min, -kinematics().cuts().initial.q2.at(0).min());
  const auto t1_limits = Limits{t1_min, t1_max};
  CG_DEBUG_LOOP("LPAIR:pickin") << "t1 in range: " << t1_limits << ".";
  const auto [t1_val, t1_width] = map_expo(m_u_t1_, t1_limits);  // definition of the first photon propagator (t1 < 0)
  t1() = t1_val;
  if (t1_width >= 0.)
    return 0.;

  const auto r1 = s2_ - t1() + mB2(), r2 = s2_ - m_w4_ + mY2(),
             rl4 = (r1 * r1 - 4. * s2_ * mB2()) * (r2 * r2 - 4. * s2_ * mY2());
  if (!utils::positive(rl4)) {
    CG_WARNING("LPAIR:pickin") << "Invalid rl4 = " << rl4 << ".";
    return 0.;
  }

  const auto d4 = m_w4_ - t1();
  const auto w52 = mY2() - mB2();
  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  auto t2_max = mB2() + mY2() - 0.5 * (r1 * r2 + std::sqrt(rl4)) / s2_,
       t2_min = (w52 * d4 + (d4 - w52) * (d4 * mB2() - w52 * t1()) / s2_) / t2_max;
  t2_max = std::max(t2_max, -kinematics().cuts().initial.q2.at(1).max());
  t2_min = std::min(t2_min, -kinematics().cuts().initial.q2.at(1).min());
  const auto t2_limits = Limits{t2_min, t2_max};
  CG_DEBUG_LOOP("LPAIR:pickin") << "t2 in range: " << t2_limits << ".";
  const auto [t2_val, t2_width] = map_expo(m_u_t2_, t2_limits);  // definition of the second photon propagator (t2 < 0)
  t2() = t2_val;
  if (t2_width >= 0.)
    return 0.;

  const auto r3 = m_w4_ - t1() - t2();
  if (gamma4_ = t1() * t2() - 0.25 * r3 * r3; gamma4_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "gamma4 = " << gamma4_ << " >= 0";
    return 0.;
  }

  const auto compute_deltas = [this](double var,
                                     short sign,
                                     double t_1,
                                     double mi2_1,
                                     double mf2_1,
                                     double t_2,
                                     double mi2_2,
                                     double mf2_2,
                                     std::array<double, 2>& deltas) -> std::tuple<double, double> {
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
    if (mi2_1 == 0.) {
      var_max =
          (s() * (t_1 * (s() + del1 - mf2_1) - mi2_2 * mf2_1) + mi2_2 * mf2_1 * (mf2_1 - del1)) / ((s() + w12_) * del2);
      deltas[0] = -0.25 * (var_max - var) * ss_ * del2;
    } else {
      const auto inv_w1 = 1. / mi2_1;
      const auto sb = mf2_1 + 0.5 * (s() * (t_1 - m2diff) + w12_ * del2) * inv_w1,
                 sd = sl1_ * std::sqrt(-sa_1) * inv_w1,
                 se = (s() * (t_1 * (s() + del2 - mi2_2) - mi2_2 * m2diff) + mf2_1 * (mi2_2 * mf2_1 + w12_ * del1)) *
                      inv_w1;
      std::tie(var_pm, var_max) = compute_boundaries(sb, sd, se);
      deltas[0] = -0.25 * (var_max - var) * (var_pm - var) * mi2_1;
    }
    {
      const auto inv_t = 1. / t_2;
      const auto sb = mi2_2 + t_1 - 0.5 * (m_w4_ - t_1 - t_2) * (mf2_2 - mi2_2 - t_2) * inv_t,
                 sd = 2. * sign * std::sqrt(sa_2 * gamma4_) * inv_t,
                 se = del3 * del1 + (del3 - del1) * (del3 * mi2_2 - del1 * mf2_2) * inv_t;
      std::tie(var_mp, var_min) = compute_boundaries(sb, sd, se);
      deltas[1] = -0.25 * (var_min - var) * (var_mp - var) * t_2;
    }
    return std::make_tuple(sa_1, 0.5 * (var - t_1 - mi2_2));
  };

  if (std::tie(sa1_, p2k1_) = compute_deltas(s2_, -1, t1(), mA2(), mX2(), t2(), mB2(), mY2(), deltas1_); sa1_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "sa1_ = " << sa1_ << " >= 0";
    return 0.;
  }

  const auto dd = deltas1_[0] * deltas1_[1];
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
           ((mB2() * r3 + 0.5 * (w52 - t2()) * st) * (p12_ * t1() - 0.25 * (t1() - w31) * st) -
            std::cos(m_theta4_) * st * std::sqrt(dd)) *
           inv_ap;

  s1_ = t2() + mA2() + 2. * (p12_ * r3 - 2. * delta_) / st;

  const auto jacobian = s2_width * t1_width * t2_width * 0.125 * 0.5 / (sl1_ * std::sqrt(-ap));
  if (!utils::positive(jacobian)) {
    CG_WARNING("LPAIR:pickin") << "Null Jacobian.\n\t"
                               << "ds2=" << s2_width << ", dt1=" << t1_width << ", dt2=" << t2_width << ".";
    return 0.;
  }
  CG_DEBUG_LOOP("LPAIR:pickin") << "s1=" << s1_ << ", s2=" << s2_ << ", ds2=" << s2_width << ", t1=" << t1()
                                << ", dt1=" << t1_width << ", t2=" << t2() << ", dt2=" << t2_width << "\n\t"
                                << "Jacobian=" << std::scientific << jacobian
                                << ", LPAIR original dj=" << (jacobian * M_PI * M_PI * 2.) << std::fixed;

  gram_ = std::pow(std::sin(m_theta4_), 2) * dd * inv_ap;

  if (std::tie(sa2_, p1k2_) = compute_deltas(s1_, +1, t2(), mB2(), mY2(), t1(), mA2(), mX2(), deltas2_); sa2_ >= 0.) {
    CG_WARNING("LPAIR:pickin") << "sa2_ = " << sa2_ << " >= 0";
    return 0.;
  }
  CG_DEBUG_LOOP("LPAIR:pickin") << std::scientific << "deltas = " << deltas1_ << ", " << deltas2_ << std::fixed;

  if (dd5_ = deltas1_[0] + deltas2_[0] +
             ((p12_ * (t1() - w31) * 0.5 - mA2() * p2k1_) * (p2k1_ * (t2() - w52) - mB2() * r3) -
              delta_ * (2. * p12_ * p2k1_ - mB2() * (t1() - w31))) /
                 p2k1_;
      !utils::positive(dd5_)) {
    CG_WARNING("LPAIR:pickin") << "Invalid dd5=" << dd5_ << ", with all deltas=" << deltas1_ << ", " << deltas2_ << ".";
    return 0.;
  }

  return jacobian;
}

//---------------------------------------------------------------------------------------------

bool LPAIR::orient() {
  eph1_ = re_ * (s2_ - mX2() + w12_);  // de3 in original LPAIR
  eph2_ = re_ * (s1_ - mY2() - w12_);  // de5 un original LPAIR

  //----- central two-photon/lepton system
  if (ec4_ = eph1_ + eph2_; ec4_ < mc4_) {
    CG_WARNING("LPAIR:orient") << "ec4_ = " << ec4_ << " < mc4_ = " << mc4_ << "==> photon energies: " << eph1_ << ", "
                               << eph2_ << ".";
    return false;
  }
  if (pc4_ = utils::fastSqrtSqDiff(ec4_, mc4_); pc4_ == 0.) {  // protons' momenta are not along the z-axis
    CG_WARNING("LPAIR:orient") << "pzc4 is null and should not be...";
    return false;
  }

  CG_DEBUG_LOOP("LPAIR:orient") << "Central system's energy: E4 = " << ec4_ << "\n\t"
                                << "               momentum: p4 = " << pc4_ << "\n\t"
                                << "         invariant mass: m4 = " << mc4_ << ".";

  pt4_ = mom_prefactor_ * std::sqrt(dd5_);
  if (sin_theta4_ = pt4_ / pc4_; !Limits{-1., 1.}.contains(sin_theta4_)) {
    CG_WARNING("LPAIR:orient") << "Invalid sin(theta4): " << sin_theta4_ << ".";
    return false;
  }
  const auto p14 = +0.5 * (s1_ + t1() - t2() - mX2());
  cos_theta4_ = std::sqrt(1. - sin_theta4_ * sin_theta4_) * (ep1_ * ec4_ < p14 ? -1. : 1.);
  const auto sin2_theta4 = sin_theta4_ * sin_theta4_;
  alpha4_ = 1. - cos_theta4_;
  beta4_ = 1. + cos_theta4_;
  if (cos_theta4_ < 0.)
    beta4_ = sin2_theta4 / alpha4_;
  else
    alpha4_ = sin2_theta4 / beta4_;

  CG_DEBUG_LOOP("LPAIR:orient") << "cos(theta4) = " << cos_theta4_ << "\t"
                                << "sin(theta4) = " << sin_theta4_ << "\n\t"
                                << "alpha4 = " << alpha4_ << ", beta4 = " << beta4_;

  //----- outgoing beam states
  const auto rr = mom_prefactor_ * std::sqrt(-gram_) / pt4_;

  //--- beam 1 -> 3
  const auto ep3 = ep1_ - eph1_, pp3 = utils::fastSqrtSqDiff(ep3, mX()), pt3 = mom_prefactor_ * std::sqrt(deltas1_[0]);
  if (pt3 > pp3) {
    CG_DEBUG("LPAIR:orient") << "Invalid momentum for outgoing beam 1: pt=" << pt3 << ", p=" << pp3 << ".";
    return false;
  }
  if (pt3 < rr) {
    CG_DEBUG("LPAIR:orient") << "Invalid momentum balance for outgoing beam 1.";
    return false;
  }
  pX() = Momentum::fromPThetaPhiE(pp3, -std::asin(pt3 / pp3), std::asin(-rr / pt3), ep3);
  CG_DEBUG_LOOP("LPAIR:orient") << "Positive-z beam state:\n\t" << std::scientific << "energy: E3 = " << ep3
                                << ", pt3 = " << pt3 << "\n\t"
                                << "momentum = " << pX() << ".";

  //--- beam 2 -> 5
  const auto ep5 = ep2_ - eph2_, pp5 = utils::fastSqrtSqDiff(ep5, mY()), pt5 = mom_prefactor_ * std::sqrt(deltas2_[0]);
  if (pt5 > pp5) {
    CG_DEBUG("LPAIR:orient") << "Invalid momentum for outgoing beam 2: pt=" << pt5 << ", p=" << pp5 << ".";
    return false;
  }
  if (pt5 < rr) {
    CG_DEBUG("LPAIR:orient") << "Invalid momentum balance for outgoing beam 2.";
    return false;
  }
  pY() = Momentum::fromPThetaPhiE(pp5, M_PI + std::asin(pt5 / pp5), std::asin(+rr / pt5), ep5);
  CG_DEBUG_LOOP("LPAIR:orient") << "Negative-z beam state:\n\t" << std::scientific << "energy: E5 = " << ep5
                                << ", pt5 = " << pt5 << "\n\t"
                                << "momentum = " << pY() << ".";

  // x-axis mirroring
  if (const double a1 = pX().px() - pY().px();
      std::fabs(pt4_ + pX().px() + pY().px()) >= std::fabs(std::fabs(a1) - pt4_)) {
    CG_DEBUG_LOOP("LPAIR:orient") << "|pt4+pt3*cos(phi3)+pt5*cos(phi5)| < | |a1|-pt4 | ; pt4 = " << pt4_ << ".";
    if (a1 < 0.)
      pY().mirrorX();
    else
      pX().mirrorX();
  }
  return true;
}

//---------------------------------------------------------------------------------------------

double LPAIR::computeWeight() {
  if (mc4_ = std::sqrt(m_w4_); !utils::positive(mc4_))  // compute the two-photon energy for this point
    return 0.;

  CG_DEBUG_LOOP("LPAIR:weight") << "Masses dump:\n\t"
                                << "m1 = " << mA() << ", m2 = " << mB() << ", m3 = " << mX() << ", m4 = " << mc4_
                                << ", m5 = " << mY() << ".\n\t"
                                << "w1 = " << mA2() << ", w2 = " << mB2() << ", w3 = " << mX2() << ", w4 = " << m_w4_
                                << ", w5 = " << mY2() << ".";

  auto jacobian = pickin();
  if (!utils::positive(jacobian)) {
    CG_DEBUG_LOOP("LPAIR:weight") << "Pickin failed.";
    return 0.;
  }
  if (!orient()) {
    CG_DEBUG_LOOP("LPAIR:weight") << "Orient failed.";
    return 0.;
  }

  const double ecm6 = m_w4_ / (2. * mc4_), pp6cm = utils::fastSqrtSqDiff(ecm6, ml_);

  jacobian *= pp6cm / mc4_;

  // Let the most obscure part of this code begin...

  const double e3mp3 = mX2() / (pX().energy() + pX().p());
  const double theta_x = pX().theta(), al3 = std::pow(std::sin(theta_x), 2) / (1. + theta_x);

  // 2-photon system kinematics ?!
  const double eg = (m_w4_ + t1() - t2()) / (2. * mc4_);

  const double gamma4 = ec4_ / mc4_;
  const Momentum pg(-pX().px() * cos_theta4_ - (pX().p() * al3 + e3mp3 - e1mp1_ + eph1_) * sin_theta4_,
                    -pX().py(),
                    -gamma4 * pX().px() * sin_theta4_ + (pX().p() * al3 + e3mp3 - e1mp1_) * gamma4 * cos_theta4_ +
                        mc4_ * eph1_ / (ec4_ + pc4_) - gamma4 * eph1_ * alpha4_);

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
             theta6cm = M_PI - std::acos(cos_theta6cm);

  // match the Jacobian
  jacobian *= (amap + bmap * cos_theta6cm);
  jacobian *= (amap - bmap * cos_theta6cm);
  jacobian *= 0.5 * std::log(ymap) / amap / bmap;
  if (symmetrise_ &&
      (beams_mode_ == mode::Kinematics::ElasticInelastic || beams_mode_ == mode::Kinematics::InelasticElastic))
    jacobian *= 1.;
  else
    jacobian *= 0.5;

  const auto p6cm = Momentum::fromPThetaPhiE(pp6cm, theta6cm, m_phi6_cm_);  // 1st outgoing lepton 3-mom. in CoM system

  const double h1 = p6cm.pz() * sin_theta_gam + p6cm.px() * cos_theta_gam;
  const double pc6z = p6cm.pz() * cos_theta_gam - p6cm.px() * sin_theta_gam;
  const double pc6x = h1 * cos_phi_gam - p6cm.py() * sin_phi_gam;
  const double qcx = 2. * pc6x, qcz = 2. * pc6z;

  const double el6 = (ec4_ * ecm6 + pc4_ * pc6z) / mc4_;
  const double h2 = (ec4_ * pc6z + pc4_ * ecm6) / mc4_;

  // outgoing leptons' kinematics (in the two-photon CM frame)
  const auto pc4 = Momentum::fromPThetaPhiE(pc4_, std::acos(cos_theta4_), 0., ec4_);
  pc(0) = Momentum(+pc6x * cos_theta4_ + h2 * sin_theta4_,
                   p6cm.py() * cos_phi_gam + h1 * sin_phi_gam,
                   -pc6x * sin_theta4_ + h2 * cos_theta4_,
                   el6);
  pc(1) = pc4 - pc(0);
  CG_DEBUG_LOOP("LPAIR") << "Outgoing kinematics\n\t"
                         << " first outgoing lepton: p = " << pc(0) << "\n\t"
                         << "second outgoing lepton: p = " << pc(1) << ".";

  bb_ = t1() * t2() + (m_w4_ * sin2_theta6cm + 4. * ml2_ * cos2_theta6cm) * p_gam * p_gam;
  q2dq_ = std::pow(eg * (2. * ecm6 - mc4_) - 2. * p_gam * p6cm.pz(), 2);

  const auto hq = ec4_ * qcz / mc4_;
  const auto qve = Momentum::fromPxPyPzE(+qcx * cos_theta4_ + hq * sin_theta4_,
                                         2. * pc(0).py(),
                                         -qcx * sin_theta4_ + hq * cos_theta4_,
                                         +qcz * pc4_ / mc4_);

  for (size_t i = 0; i < 2; ++i)  // boost outgoing leptons' kinematics into lab frame
    pc(i).betaGammaBoost(gamma_cm_, beta_gamma_cm_);
  if (!kinematics().cuts().central.contain(event()(Particle::CentralSystem)))  // cuts on outgoing leptons
    return 0.;

  {  // preparation for the periPP call
    const auto compute_coeffs =
        [&qve](double e_in, double m_in, const Momentum& pout, double ene_pho, double pcm, double& gamma)
        -> std::tuple<double, double, double, double, double, double, double> {
      const auto phi_out = pout.phi(), cos_phi_out = std::cos(phi_out), sin_phi_out = std::sin(phi_out);
      const auto e2_in = e_in * e_in, m2_in = m_in * m_in, pt_out = pout.pt();
      const auto c1 = pt_out * (qve.px() * sin_phi_out - qve.py() * cos_phi_out),
                 c2 = pt_out * (qve.pz() * e_in - qve.energy() * pcm),
                 c3 = ((pout.mass2() - m2_in) * e2_in + 2. * m2_in * ene_pho * e_in - m2_in * ene_pho * ene_pho +
                       pt_out * pt_out * e2_in) /
                      (pout.pz() * e2_in + pout.energy() * pcm);
      const auto r2 = c2 * sin_phi_out + c3 * qve.py(), r3 = -c2 * cos_phi_out - c3 * qve.px();
      gamma = m2_in * c1 * c1 + r2 * r2 + r3 * r3;
      return std::make_tuple(cos_phi_out, sin_phi_out, c1, c2, c3, r2, r3);
    };
    const auto [cos_phi3, sin_phi3, c1, c2, c3, r12, r13] = compute_coeffs(ep1_, mA(), pX(), eph1_, +p_cm_, gamma5_);
    const auto [cos_phi5, sin_phi5, b1, b2, b3, r22, r23] = compute_coeffs(ep2_, mB(), pY(), eph2_, -p_cm_, gamma6_);
    const auto pt3 = pX().pt(), pt5 = pY().pt();
    alpha5_ = -(qve.px() * cos_phi3 + qve.py() * sin_phi3) * pt3 * p1k2_ -
              (ep1_ * qve.energy() - p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
              (eph2_ * qve.pz() + qve.energy() * (p_cm_ + pY().pz())) * c3;
    alpha6_ = -(qve.px() * cos_phi5 + qve.py() * sin_phi5) * pt5 * p2k1_ -
              (ep2_ * qve.energy() + p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
              (eph1_ * qve.pz() - qve.energy() * (p_cm_ - pY().pz())) * b3;
    epsilon_ = p12_ * c1 * b1 + r12 * r22 + r13 * r23;
  }
  const auto peripp = periPP();  // compute the structure functions factors
  if (!utils::positive(peripp))
    return 0.;

  const auto alpha_prod = alphaEM(std::sqrt(-t1())) * alphaEM(std::sqrt(-t2()));
  jacobian *= constb_ * charge_factor_ * alpha_prod * alpha_prod / s();

  CG_DEBUG_LOOP("LPAIR:f") << "Jacobian: " << jacobian << ", str.fun. factor: " << peripp << ".";
  return jacobian * peripp;  // compute the event weight using the Jacobian
}
REGISTER_PROCESS("lpair", LPAIR);
