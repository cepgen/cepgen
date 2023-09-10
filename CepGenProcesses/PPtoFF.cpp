/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"

using namespace cepgen;

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process using \f$k_{\rm T}\f$-factorization approach
class PPtoFF final : public cepgen::proc::Process2to4 {
public:
  explicit PPtoFF(const ParametersList& params)
      : cepgen::proc::Process2to4(params, params.get<ParticleProperties>("pair").pdgid),
        method_(steerAs<int, Mode>("method")),
        osp_(steer<ParametersList>("offShellParameters")) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PPtoFF(*this)); }

  static ParametersDescription description() {
    auto desc = Process2to4::description();
    desc.setDescription("γγ → f⁺f¯ (kt-factor.)");
    desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("type of central particles emitted");
    desc.addAs<int, Mode>("method", Mode::offShell)
        .setDescription("Matrix element computation method (0 = on-shell, 1 = off-shell)");
    desc.add("offShellParameters", OffShellParameters::description());
    return desc;
  }

private:
  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override {
    switch (method_) {
      case Mode::onShell:
        return onShellME();
      case Mode::offShell:
        return offShellME();
      default:
        throw CG_FATAL("PPtoFF") << "Invalid ME calculation method (" << (int)method_ << ")!";
    }
  }
  double couplingPrefactor(double q_1, double q_2) const {
    double prefactor = g2_prefactor_ * g2_prefactor_;
    if (event().oneWithRole(Particle::Parton1).pdgId() == PDG::gluon)
      prefactor *= 0.5 * alphaS(q_1);
    else
      prefactor *= alphaEM(q_1);
    if (event().oneWithRole(Particle::Parton2).pdgId() == PDG::gluon)
      prefactor *= 0.5 * alphaS(q_2);
    else
      prefactor *= alphaEM(q_2);
    return prefactor;
  }

  /// Rapidity range for the outgoing fermions
  double onShellME() const;
  double offShellME() const;

  const enum class Mode { onShell = 0, offShell = 1 } method_;

  //--- parameters for off-shell matrix element
  struct OffShellParameters : SteeredObject<OffShellParameters> {
    explicit OffShellParameters(const ParametersList& params) : SteeredObject(params) {
      (*this)
          .add("mat1", mat1)
          .add("mat2", mat2)
          .add("termLL", term_ll)
          .add("termLT", term_lt)
          .add("termTT", term_tt1)
          .add("termtt", term_tt2);
    }

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add("mat1", 1).setDescription("symmetry factor for the first incoming photon");
      desc.add("mat2", 1).setDescription("symmetry factor for the second incoming photon");
      desc.add("termLL", 1).setDescription("fully longitudinal relative weight");
      desc.add("termLT", 1).setDescription("cross-polarisation relative weight");
      desc.add("termTT", 1).setDescription("fully transverse relative weight");
      desc.add("termtt", 1).setDescription("fully transverse relative weight");
      return desc;
    }

    int mat1{0}, mat2{0};
    int term_ll{0}, term_lt{0}, term_tt1{0}, term_tt2{0};
  } osp_;

  double mf2_{0.};
  /// Prefactor for the alpha(S/EM) coupling
  double g2_prefactor_{0.};
};

void PPtoFF::prepareProcessKinematics() {
  const auto& cs_prop = PDG::get()(produced_parts_.at(0));
  if (!cs_prop.fermion || cs_prop.charge == 0.)
    throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion pair selected: " << cs_prop << ".";

  mf2_ = cs_prop.mass * cs_prop.mass;
  g2_prefactor_ = 4. * M_PI;

  CG_DEBUG("PPtoFF:prepare") << "Incoming beams: mp(1/2) = " << mA() << "/" << mB() << ".\n\t"
                             << "Produced particles: " << cs_prop_ << ".\n\t"
                             << "ME computation method: " << (int)method_ << ".";

  if (!kinematics().cuts().central.pt_diff.valid())
    kinematics().cuts().central.pt_diff = {0., 50.};  // tighter cut for fermions

  for (const auto& role : {Particle::Parton1, Particle::Parton2})
    switch (event().oneWithRole(role).pdgId()) {
      case PDG::gluon:
        break;
      case PDG::photon:
        g2_prefactor_ *= std::pow(cs_prop.charge / 3., 2);  // electromagnetic coupling
        break;
      default:
        throw CG_FATAL("PPtoFF:prepare") << "Only photon & gluon partons are supported!";
    }
}

double PPtoFF::onShellME() const {
  const double s_hat = shat(), t_hat = that(), u_hat = uhat();
  CG_DEBUG_LOOP("PPtoFF:onShell") << "shat: " << s_hat << ", that: " << t_hat << ", uhat: " << u_hat << ".";

  if (t_hat == mf2_ || u_hat == mf2_)
    return 0.;
  const auto q = std::sqrt(t_hat), g2_sq = couplingPrefactor(q, q);

  const double mf4 = mf2_ * mf2_, mf8 = mf4 * mf4;

  double out = 6. * mf8;
  out += -3. * mf4 * t_hat * t_hat;
  out += -14. * mf4 * t_hat * u_hat;
  out += -3. * mf4 * u_hat * u_hat;
  out += 1. * mf2_ * t_hat * t_hat * t_hat;
  out += 7. * mf2_ * t_hat * t_hat * u_hat;
  out += 7. * mf2_ * t_hat * u_hat * u_hat;
  out += 1. * mf2_ * u_hat * u_hat * u_hat;
  out += -1. * t_hat * t_hat * t_hat * u_hat;
  out += -1. * t_hat * u_hat * u_hat * u_hat;
  return -2. * out * g2_sq * std::pow((mf2_ - t_hat) * (mf2_ - u_hat), -2);
}

double PPtoFF::offShellME() const {
  //--- positive polarisation
  const auto alpha1 = amt1_ / sqrtS() * exp(y_c1_), alpha2 = amt2_ / sqrtS() * exp(y_c2_);
  const auto x1 = alpha1 + alpha2, z1p = alpha1 / x1, z1m = alpha2 / x1, z1 = z1p * z1m;
  const auto ak1 = Momentum(z1m * pc(0) - z1p * pc(1)).setPz(0.);
  const auto ph_p1 = ak1 + z1p * q2(), ph_m1 = ak1 - z1m * q2();
  const auto qt1 = q1().p(), inv_qt1 = 1. / qt1;
  const auto t1abs = (qt1 * qt1 + x1 * (mX2() - mA2()) + x1 * x1 * mA2()) / (1. - x1);
  const auto eps12 = mf2_ + z1 * t1abs;
  const auto kp1 = 1. / (ph_p1.pt2() + eps12), km1 = 1. / (ph_m1.pt2() + eps12);

  const auto phi1 = Momentum(kp1 * ph_p1 - km1 * ph_m1).setPz(0.).setEnergy(kp1 - km1);
  const auto dot1 = phi1.threeProduct(q1()) * inv_qt1, cross1 = phi1.crossProduct(q1()) * inv_qt1;

  const auto aux2_1 = osp_.term_ll * (mf2_ + 4. * z1 * z1 * t1abs) * phi1.energy2() +
                      osp_.term_tt1 * ((z1p * z1p + z1m * z1m) * (dot1 * dot1 + cross1 * cross1)) +
                      osp_.term_tt2 * (cross1 * cross1 - dot1 * dot1) -
                      osp_.term_lt * 4. * z1 * (z1p - z1m) * phi1.energy() * qt1 * dot1;
  const auto amat2_1 = 2. * aux2_1 * z1 * qt1 * qt1;

  //--- negative polarisation
  const auto beta1 = amt1_ / sqrtS() * exp(-y_c1_), beta2 = amt2_ / sqrtS() * exp(-y_c2_);
  const auto x2 = beta1 + beta2, z2p = beta1 / x2, z2m = beta2 / x2, z2 = z2p * z2m;
  const auto ak2 = Momentum(z2m * pc(0) - z2p * pc(1)).setPz(0.);
  const auto ph_p2 = ak2 + z2p * q1(), ph_m2 = ak2 - z2m * q1();
  const auto qt2 = q2().p(), inv_qt2 = 1. / qt2;
  const auto t2abs = (qt2 * qt2 + x2 * (mY2() - mB2()) + x2 * x2 * mB2()) / (1. - x2);
  const auto eps22 = mf2_ + z2 * t2abs;
  const auto kp2 = 1. / (ph_p2.pt2() + eps22), km2 = 1. / (ph_m2.pt2() + eps22);

  const auto phi2 = Momentum(kp2 * ph_p2 - km2 * ph_m2).setPz(0.).setEnergy(kp2 - km2);
  const auto dot2 = phi2.threeProduct(q2()) * inv_qt2, cross2 = phi2.crossProduct(q2()) * inv_qt2;

  const auto aux2_2 = osp_.term_ll * (mf2_ + 4. * z2 * z2 * t2abs) * phi2.energy2() +
                      osp_.term_tt1 * ((z2p * z2p + z2m * z2m) * (dot2 * dot2 + cross2 * cross2)) +
                      osp_.term_tt2 * (cross2 * cross2 - dot2 * dot2) -
                      osp_.term_lt * 4. * z2 * (z2p - z2m) * phi2.energy() * qt2 * dot2;
  const auto amat2_2 = 2. * aux2_2 * z2 * qt2 * qt2;

  //--- symmetrisation
  const auto amat2 = 0.5 * (osp_.mat1 * amat2_1 + osp_.mat2 * amat2_2);
  if (amat2 <= 0.)
    return 0.;

  CG_DEBUG_LOOP("PPtoFF:offShell") << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
                                   << "z(1/2) = " << z1 << " / " << z2 << "\n\t"
                                   << "q(1/2)t2 = " << q1().pt2() << " / " << q2().pt2() << "\n\t"
                                   << "amat2(1/2) = " << amat2_1 << " / " << amat2_2 << "\n\t"
                                   << "amat2 = " << amat2 << ".";

  const double tmax = pow(std::max(amt1_, amt2_), 2);
  const double q1val = std::sqrt(std::max(eps12, tmax)), q2val = std::sqrt(std::max(eps22, tmax));
  const auto g2_sq = couplingPrefactor(q1val, q2val);

  return g2_sq * std::pow(x1 * x2 * s() * inv_qt1 * inv_qt2, 2) * amat2;
}
// register process
REGISTER_PROCESS("pptoff", PPtoFF);
