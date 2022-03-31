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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace proc {
    const Limits Process2to4::x_limits_{0., 1.};

    Process2to4::Process2to4(const ParametersList& params, std::array<pdgid_t, 2> partons, pdgid_t cs_id)
        : KTProcess(params, partons, {cs_id, cs_id}), cs_prop_(PDG::get()(cs_id)), single_limits_(params) {}

    Process2to4::Process2to4(const Process2to4& proc)
        : KTProcess(proc), cs_prop_(proc.cs_prop_), single_limits_(proc.single_limits_) {}

    void Process2to4::setCuts(const cuts::Central& single) { single_limits_ = single; }

    void Process2to4::preparePhaseSpace() {
      if (cs_prop_.pdgid == PDG::invalid)  // ensure the central particles properties are correctly initialised
        cs_prop_ = PDG::get()(steer<ParticleProperties>("pair").pdgid);
      {
        const auto& beamA = event_->oneWithRole(Particle::IncomingBeam1);
        pA_ = beamA.momentum();
        mA2_ = beamA.mass2();
      }
      {
        const auto& beamB = event_->oneWithRole(Particle::IncomingBeam2);
        pB_ = beamB.momentum();
        mB2_ = beamB.mass2();
      }
      CG_DEBUG_LOOP("2to4:incoming") << "incoming particles:\n"
                                     << "  pA = " << pA_ << ", mA2 = " << mA2_ << "\n"
                                     << "  pB = " << pB_ << ", mB2 = " << mB2_ << ".";

      ww_ = 0.5 * (1. + sqrt(1. - 4. * sqrt(mA2_ * mB2_) / s_));

      defineVariable(
          y_c1_, Mapping::linear, kin_.cuts().central.rapidity_single(), {-6., 6.}, "First outgoing particle rapidity");
      defineVariable(
          y_c2_, Mapping::linear, kin_.cuts().central.rapidity_single(), {-6., 6.}, "Second outgoing particle rapidity");
      defineVariable(pt_diff_,
                     Mapping::linear,
                     kin_.cuts().central.pt_diff(),
                     {0., 500.},
                     "Final state particles transverse momentum difference");
      defineVariable(phi_pt_diff_,
                     Mapping::linear,
                     kin_.cuts().central.phi_diff(),
                     {0., 2. * M_PI},
                     "Final state particles azimuthal angle difference");

      prepareProcessKinematics();
    }

    double Process2to4::computeKTFactorisedMatrixElement() {
      //--- transverse kinematics of initial partons
      const auto qt_1 = Momentum::fromPtEtaPhiE(qt1_, 0., phi_qt1_);
      if (fabs(qt_1.pt() - qt1_) > NUM_LIMITS)
        throw CG_FATAL("Process2to4") << "|qt1|=" << qt1_ << " != qt1.pt()=" << qt_1.pt() << ", qt1=" << qt_1 << ".";

      const auto qt_2 = Momentum::fromPtEtaPhiE(qt2_, 0., phi_qt2_);
      if (fabs(qt_2.pt() - qt2_) > NUM_LIMITS)
        throw CG_FATAL("Process2to4") << "|qt2|=" << qt1_ << " != qt2.pt()=" << qt_2.pt() << ", qt2=" << qt_2 << ".";

      //--- two-parton system (in transverse plane)
      const auto qt_sum = qt_1 + qt_2;

      CG_DEBUG_LOOP("2to4:me") << "q(1/2)x = " << qt_1.px() << " / " << qt_2.px() << "\n\t"
                               << "q(1/2)y = " << qt_1.py() << " / " << qt_2.py() << "\n\t"
                               << "sum(qt) = " << qt_sum;

      //--- transverse kinematics of outgoing central system
      const auto pt_diff = Momentum::fromPtEtaPhiE(pt_diff_, 0., phi_pt_diff_);
      if (fabs(pt_diff.pt() - pt_diff_) > NUM_LIMITS)
        throw CG_FATAL("Process2to4") << "|dpt|=" << pt_diff_ << " != dpt.pt()=" << pt_diff.pt() << ", dpt=" << pt_diff
                                      << ".";

      const auto pt_c1 = 0.5 * (qt_sum + pt_diff);
      const auto pt_c2 = 0.5 * (qt_sum - pt_diff);
      const double p1t = pt_c1.pt(), p2t = pt_c2.pt();

      CG_DEBUG_LOOP("2to4:me") << "diff(pt) = " << pt_diff << "\n\t"
                               << "p(1/2)x = " << pt_c1.px() << " / " << pt_c2.px() << "\n\t"
                               << "p(1/2)y = " << pt_c1.py() << " / " << pt_c2.py() << "\n\t"
                               << "p(1/2)t = " << p1t << " / " << p2t;

      //--- window in rapidity distance
      if (!kin_.cuts().central.rapidity_diff().contains(fabs(y_c1_ - y_c2_)))
        return 0.;

      //--- apply the pt cut already at this stage (remains unchanged)
      if (!kin_.cuts().central.pt_single().contains(p1t))
        return 0.;
      if (!kin_.cuts().central.pt_single().contains(p2t))
        return 0.;
      if (!single_limits_.pt_single().contains(p1t))
        return 0.;
      if (!single_limits_.pt_single().contains(p2t))
        return 0.;

      //--- window in transverse momentum difference
      if (!kin_.cuts().central.pt_diff().contains(fabs(p1t - p2t)))
        return 0.;

      //--- transverse mass for the two central particles
      amt1_ = std::hypot(p1t, cs_prop_.mass);
      amt2_ = std::hypot(p2t, cs_prop_.mass);

      //--- window in central system invariant mass
      const double invm = sqrt(amt1_ * amt1_ + amt2_ * amt2_ + 2. * amt1_ * amt2_ * cosh(y_c1_ - y_c2_) - qt_sum.pt2());
      if (!kin_.cuts().central.mass_sum().contains(invm))
        return 0.;

      //--- auxiliary quantities

      const double alpha1 = amt1_ / sqs_ * exp(y_c1_), beta1 = amt1_ / sqs_ * exp(-y_c1_);
      const double alpha2 = amt2_ / sqs_ * exp(y_c2_), beta2 = amt2_ / sqs_ * exp(-y_c2_);
      x1_ = alpha1 + alpha2;
      x2_ = beta1 + beta2;

      CG_DEBUG_LOOP("2to4:sudakov") << "Sudakov parameters:\n\t"
                                    << "  alpha(1/2) = " << alpha1 << " / " << alpha2 << "\n\t"
                                    << "   beta(1/2) = " << beta1 << " / " << beta2 << ".";

      //--- sanity check for x_i values
      if (!x_limits_.contains(x1_) || !x_limits_.contains(x2_))
        return 0.;

      //--- additional conditions for energy-momentum conservation

      const double s1_eff = x1_ * s_ - qt1_ * qt1_, s2_eff = x2_ * s_ - qt2_ * qt2_;

      CG_DEBUG_LOOP("2to4:central") << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
                                    << "central system invariant mass = " << invm << " GeV";

      if (kin_.incomingBeams().positive().fragmented() && (sqrt(s2_eff) <= sqrt(mX2_) + invm))
        return 0.;
      if (kin_.incomingBeams().negative().fragmented() && (sqrt(s1_eff) <= sqrt(mY2_) + invm))
        return 0.;

      //--- four-momenta of the outgoing protons (or remnants)

      const double px_plus = (1. - x1_) * pA_.p() * M_SQRT2;
      const double py_minus = (1. - x2_) * pB_.p() * M_SQRT2;
      const double px_minus = (mX2_ + qt1_ * qt1_) * 0.5 / px_plus;
      const double py_plus = (mY2_ + qt2_ * qt2_) * 0.5 / py_minus;
      // warning! sign of pz??

      CG_DEBUG_LOOP("2to4:pxy") << "px± = " << px_plus << " / " << px_minus << "\n\t"
                                << "py± = " << py_plus << " / " << py_minus << ".";

      pX_ = -Momentum(qt_1).setPz((px_plus - px_minus) * M_SQRT1_2).setEnergy((px_plus + px_minus) * M_SQRT1_2);

      pY_ = -Momentum(qt_2).setPz((py_plus - py_minus) * M_SQRT1_2).setEnergy((py_plus + py_minus) * M_SQRT1_2);

      CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << pX_ << ", mass = " << pX_.mass() << "\n\t"
                                     << "Second remnant: " << pY_ << ", mass = " << pY_.mass() << ".";

      if (fabs(pX_.mass2() - mX2_) > NUM_LIMITS) {
        CG_WARNING("2to4:px") << "Invalid X system squared mass: " << pX_.mass2() << "/" << mX2_ << ".";
        return 0.;
      }
      if (fabs(pY_.mass2() - mY2_) > NUM_LIMITS) {
        CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << pY_.mass2() << "/" << mY2_ << ".";
        return 0.;
      }

      //--- four-momenta of the intermediate partons

      const double norm = 1. / ww_ / ww_ / s_;
      const double tau1 = norm * qt1_ * qt1_ / x1_ / x1_;
      q1_ =
          Momentum(qt_1).setPz(+0.5 * x1_ * ww_ * sqs_ * (1. - tau1)).setEnergy(+0.5 * x1_ * ww_ * sqs_ * (1. + tau1));

      const double tau2 = norm * qt2_ * qt2_ / x2_ / x2_;
      q2_ =
          Momentum(qt_2).setPz(-0.5 * x2_ * ww_ * sqs_ * (1. - tau2)).setEnergy(+0.5 * x2_ * ww_ * sqs_ * (1. + tau2));

      CG_DEBUG_LOOP("2to4:partons") << "First parton:  " << q1_ << ", mass2 = " << q1_.mass2() << "\n\t"
                                    << "Second parton: " << q2_ << ", mass2 = " << q2_.mass2() << ".";

      //--- four-momenta of the outgoing central particles

      p_c1_ = (pt_c1 + alpha1 * pA_ + beta1 * pB_).setEnergy(alpha1 * pA_.energy() + beta1 * pB_.energy());
      p_c2_ = (pt_c2 + alpha2 * pA_ + beta2 * pB_).setEnergy(alpha2 * pA_.energy() + beta2 * pB_.energy());

      CG_DEBUG_LOOP("2to4:central") << "First central particle:  " << p_c1_ << ", mass = " << p_c1_.mass() << "\n\t"
                                    << "Second central particle: " << p_c2_ << ", mass = " << p_c2_.mass() << ".";

      //--- compute the central 2-to-2 matrix element

      const double amat2 = computeCentralMatrixElement();
      if (amat2 <= 0.)  // skip computing the fluxes if no contribution
        return 0.;

      //=================================================================
      // factor 1/4 from jacobian of transformations
      // factors 1/pi and 1/pi due to integration over
      //     d^2(kappa_1)d^2(kappa_2) instead of d(kappa_1^2)d(kappa_2^2)
      //=================================================================

      return amat2 * pow(4. * x1_ * x2_ * s_ * M_PI, -2) * 0.25 * constants::GEVM2_TO_PB * pt_diff_ * qt1_ * qt2_;
    }

    void Process2to4::fillCentralParticlesKinematics() {
      //--- randomise the charge of outgoing system
      short sign = (drand() > 0.5) ? +1 : -1;

      //--- first outgoing central particle
      auto& oc1 = (*event_)[Particle::CentralSystem][0];
      oc1.setPdgId(cs_prop_.pdgid, +sign);
      oc1.setStatus(Particle::Status::Undecayed);
      oc1.setMomentum(p_c1_);

      //--- second outgoing central particle
      auto& oc2 = (*event_)[Particle::CentralSystem][1];
      oc2.setPdgId(cs_prop_.pdgid, -sign);
      oc2.setStatus(Particle::Status::Undecayed);
      oc2.setMomentum(p_c2_);
    }

    //----- utilities

    double Process2to4::shat() const { return (q1_ + q2_).mass2(); }

    double Process2to4::that() const {
      const double that1 = (q1_ - p_c1_).mass2();
      const double that2 = (q2_ - p_c2_).mass2();
      return 0.5 * (that1 + that2);
    }

    double Process2to4::uhat() const {
      const double uhat1 = (q1_ - p_c2_).mass2();
      const double uhat2 = (q2_ - p_c1_).mass2();
      return 0.5 * (uhat1 + uhat2);
    }

    ParametersDescription Process2to4::description() {
      auto desc = KTProcess::description();
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
