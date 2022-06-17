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
#include "CepGen/Event/Event.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenProcesses/LPAIR.h"

namespace cepgen {
  namespace proc {
    LPAIR::LPAIR(const ParametersList& params)
        : Process(params, true),
          n_opt_(steer<int>("nopt")),
          pair_(steer<int>("pair")),
          symmetrise_(steer<bool>("symmetrise")),
          rnd_phi_(0., 2. * M_PI),
          rnd_side_(0, 1) {
      if (params_.has<ParticleProperties>("pair"))
        pair_ = steer<ParticleProperties>("pair").pdgid;
    }

    LPAIR::LPAIR(const LPAIR& proc)
        : Process(proc),
          n_opt_(proc.n_opt_),
          pair_(proc.pair_),
          symmetrise_(proc.symmetrise_),
          rnd_phi_(proc.rnd_phi_),
          rnd_side_(proc.rnd_side_) {}

    //---------------------------------------------------------------------------------------------

    void LPAIR::addEventContent() {
      Process::setEventContent({{Particle::IncomingBeam1, PDG::proton},
                                {Particle::IncomingBeam2, PDG::proton},
                                {Particle::Parton1, PDG::photon},
                                {Particle::Parton2, PDG::photon}},
                               {{Particle::OutgoingBeam1, {PDG::proton}},
                                {Particle::OutgoingBeam2, {PDG::proton}},
                                {Particle::CentralSystem, {pair_, pair_}}});
    }

    //---------------------------------------------------------------------------------------------

    void LPAIR::prepareKinematics() {
      masses_.Ml2 = (*event_)(Particle::CentralSystem)[0].mass2();

      //--- first define the squared mass range for the diphoton/dilepton system
      const auto& mll_limits = kin_.cuts().central.mass_sum();
      w_limits_ = Limits(mll_limits.hasMin() ? std::pow(mll_limits.min(), 2) : 4. * masses_.Ml2,
                         mll_limits.hasMax() ? std::pow(mll_limits.max(), 2) : s_);

      CG_DEBUG_LOOP("LPAIR:prepareKinematics") << "w limits = " << w_limits_ << "\n\t"
                                               << "wmax/wmin = " << w_limits_.max() / w_limits_.min();

      p1_lab_ = (*event_)(Particle::IncomingBeam1)[0].momentum();
      p2_lab_ = (*event_)(Particle::IncomingBeam2)[0].momentum();

      const double mx0 = mp_ + PDG::get().mass(PDG::piPlus);  // 1.07
      const double min_wx = pow(std::max(mx0, kin_.cuts().remnants.mx().min()), 2);
      const Limits wx_lim_ob1(
          min_wx, pow(std::min(sqs_ - p1_lab_.mass() - 2. * sqrt(masses_.Ml2), kin_.cuts().remnants.mx().max()), 2));
      const Limits wx_lim_ob2(
          min_wx, pow(std::min(sqs_ - p2_lab_.mass() - 2. * sqrt(masses_.Ml2), kin_.cuts().remnants.mx().max()), 2));

      //--- variables mapping

      defineVariable(u_t1_, Mapping::linear, {0., 1.}, {0., 1.}, "u_t1");
      defineVariable(u_t2_, Mapping::linear, {0., 1.}, {0., 1.}, "u_t2");
      defineVariable(u_s2_, Mapping::linear, {0., 1.}, {0., 1.}, "u_s2");
      defineVariable(w4_, Mapping::power_law, w_limits_, w_limits_, "w4");
      defineVariable(theta4_, Mapping::linear, {0., M_PI}, {0., M_PI}, "theta4");
      defineVariable(phi6_cm_, Mapping::linear, {0., 2. * M_PI}, {0., 2. * M_PI}, "phi6cm");
      defineVariable(x6_, Mapping::linear, {0., 1.}, {0., 1.}, "x6");

      //--- first outgoing beam particle or remnant mass
      switch (kin_.incomingBeams().positive().mode()) {
        case Beam::Mode::PointLikeFermion:
        case Beam::Mode::ProtonElastic:
          event_->oneWithRole(Particle::OutgoingBeam1).setPdgId(event_->oneWithRole(Particle::IncomingBeam1).pdgId());
          mX2_ = p1_lab_.mass2();
          break;
        case Beam::Mode::ProtonInelastic:
          defineVariable(mX2_, Mapping::power_law, wx_lim_ob1, wx_lim_ob1, "MX2");
          break;
        default:
          throw CG_FATAL("LPAIR:kinematics")
              << "Invalid mode for beam 1: " << kin_.incomingBeams().positive().mode() << " is not supported!";
      }
      //--- second outgoing beam particle or remnant mass
      switch (kin_.incomingBeams().negative().mode()) {
        case Beam::Mode::PointLikeFermion:
        case Beam::Mode::ProtonElastic:
          event_->oneWithRole(Particle::OutgoingBeam2).setPdgId(event_->oneWithRole(Particle::IncomingBeam2).pdgId());
          mY2_ = p2_lab_.mass2();
          break;
        case Beam::Mode::ProtonInelastic:
          defineVariable(mY2_, Mapping::power_law, wx_lim_ob2, wx_lim_ob2, "MY2");
          break;
        default:
          throw CG_FATAL("LPAIR:kinematics")
              << "Invalid mode for beam 2: " << kin_.incomingBeams().negative().mode() << " is not supported!";
      }
    }

    //---------------------------------------------------------------------------------------------

    bool LPAIR::pickin() {
      CG_DEBUG_LOOP("LPAIR") << "Optimised mode? " << n_opt_;

      jacobian_ = 0.;

      // min(s2) = sigma and sig2 = sigma' in [1]
      const double sig = mc4_ + sqrt(mY2_);
      Limits s2_range(sig * sig, s_ + mX2_ - 2 * sqrt(mX2_) * sqs_);

      CG_DEBUG_LOOP("LPAIR") << "mc4 = " << mc4_ << "\n\t"
                             << "s2 in range " << s2_range << ".";

      const double d6 = w4_ - mY2_;

      CG_DEBUG_LOOP("LPAIR") << "w1 = " << mA2_ << ", w2 = " << mB2_ << ", w3 = " << mX2_ << ", w4 = " << w4_
                             << ", w5 = " << mY2_ << ". w31 = " << masses_.w31 << ", w52 = " << masses_.w52
                             << ", w12 = " << masses_.w12 << ".";

      const double ss = s_ + masses_.w12;

      const double rl1 = ss * ss - 4. * mA2_ * s_;  // lambda(s, m1**2, m2**2)
      if (rl1 <= 0.) {
        CG_DEBUG_LOOP("LPAIR") << "rl1 = " << rl1 << " <= 0";
        return false;
      }
      sl1_ = sqrt(rl1);

      s2_ = 0.;
      double ds2 = 0.;
      if (n_opt_ == 0) {
        const auto s2 = map(u_s2_, s2_range, "s2");
        s2_range.min() = s2_ = s2.first;  // why lower s2 range update?
        ds2 = s2.second;
      }

      CG_DEBUG_LOOP("LPAIR") << "s2 = " << s2_;

      const double sp = s_ + mX2_ - s2_range.min(), d3 = s2_range.min() - mB2_;
      const double rl2 = sp * sp - 4. * s_ * mX2_;  // lambda(s, m3**2, sigma)
      if (rl2 <= 0.) {
        CG_DEBUG_LOOP("LPAIR") << "rl2 = " << rl2 << " <= 0";
        return false;
      }
      const double sl2 = sqrt(rl2);

      double t1_max = mA2_ + mX2_ - (ss * sp + sl1_ * sl2) / (2. * s_);  // definition from eq. (A.4) in [1]
      double t1_min = (masses_.w31 * d3 + (d3 - masses_.w31) * (d3 * mA2_ - masses_.w31 * mB2_) / s_) /
                      t1_max;  // definition from eq. (A.5) in [1]

      // FIXME dropped in CDF version
      if (t1_max > -kin_.cuts().initial.q2().min()) {
        CG_DEBUG_LOOP("LPAIR") << "t1max = " << t1_max << " > -q2min = " << -kin_.cuts().initial.q2().min();
        return false;
      }
      if (t1_min < -kin_.cuts().initial.q2().max() && kin_.cuts().initial.q2().hasMax()) {
        CG_DEBUG_LOOP("LPAIR") << "t1min = " << t1_min << " < -q2max = " << -kin_.cuts().initial.q2().max();
        return false;
      }
      if (t1_max < -kin_.cuts().initial.q2().max() && kin_.cuts().initial.q2().hasMax())
        t1_max = -kin_.cuts().initial.q2().max();
      if (t1_min > -kin_.cuts().initial.q2().min() && kin_.cuts().initial.q2().hasMin())
        t1_min = -kin_.cuts().initial.q2().min();
      /////

      // t1, the first photon propagator, is defined here
      const Limits t1_range(t1_min, t1_max);
      const auto t1 = map(u_t1_, t1_range, "t1");
      t1_ = t1.first;
      const double dt1 = -t1.second;  // changes wrt mapt1 : dx->-dx

      CG_DEBUG_LOOP("LPAIR") << "Definition of t1 = " << t1_ << " in range " << t1_range << ".";

      deltas_[3] = w4_ - t1_;

      const double d8 = t1_ - mB2_;
      const double t13 = t1_ - mA2_ - mX2_;

      sa1_ = -pow(t1_ - masses_.w31, 2) / 4. + mA2_ * t1_;
      if (sa1_ >= 0.) {
        CG_WARNING("LPAIR") << "sa1_ = " << sa1_ << " >= 0";
        return false;
      }

      const double sl3 = sqrt(-sa1_);

      s2_range.min() = sig * sig;
      // one computes splus and (s2x=s2max)
      double splus;
      if (mA2_ != 0.) {
        const double inv_w1 = 1. / mA2_;
        const double sb = mX2_ + 0.5 * (s_ * (t1_ - masses_.w31) + masses_.w12 * t13) * inv_w1;
        const double sd = sl1_ * sl3 * inv_w1;
        const double se =
            (s_ * (t1_ * (s_ + t13 - mB2_) - mB2_ * masses_.w31) + mX2_ * (masses_.w12 * d8 + mB2_ * mX2_)) * inv_w1;

        if (fabs((sb - sd) / sd) >= 1.) {
          splus = sb - sd;
          s2_range.max() = se / splus;
        } else {
          s2_range.max() = sb + sd;
          splus = se / s2_range.max();
        }
      } else {  // 3
        s2_range.max() = (s_ * (t1_ * (s_ + d8 - mX2_) - mB2_ * mX2_) + mB2_ * mX2_ * (mB2_ + mX2_ - t1_)) / (ss * t13);
        splus = s2_range.min();
      }
      // 4
      double s2x = s2_range.max();

      CG_DEBUG_LOOP("LPAIR") << "s2x = s2max = " << s2x;

      if (n_opt_ < 0) {  // 5
        if (splus > s2_range.min()) {
          s2_range.min() = splus;
          CG_DEBUG_LOOP("LPAIR") << "min(s2) truncated to splus = " << splus;
        }
        const auto s2 = n_opt_ < -1 ? map(u_s2_, s2_range, "s2") : mapla(t1_, mB2_, u_s2_, s2_range);  // n_opt_==-1
        s2_ = s2.first;
        ds2 = s2.second;
        s2x = s2_;
      } else if (n_opt_ == 0)
        s2x = s2_;  // 6

      CG_DEBUG_LOOP("LPAIR") << "s2x = " << s2x;

      // 7
      const double r1 = s2x - d8, r2 = s2x - d6;

      const double rl4 = (r1 * r1 - 4. * mB2_ * s2x) * (r2 * r2 - 4. * mY2_ * s2x);
      if (rl4 <= 0.) {
        CG_DEBUG_LOOP("LPAIR") << "rl4 = " << rl4 << " <= 0";
        return false;
      }
      const double sl4 = sqrt(rl4);

      // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
      const double t2_max = mB2_ + mY2_ - (r1 * r2 + sl4) / s2x * 0.5,
                   t2_min = (masses_.w52 * deltas_[3] +
                             (deltas_[3] - masses_.w52) * (deltas_[3] * mB2_ - masses_.w52 * t1_) / s2x) /
                            t2_max;

      // t2, the second photon propagator, is defined here
      const auto t2 = map(u_t2_, Limits(t2_min, t2_max), "t2");
      t2_ = t2.first;
      const double dt2 = -t2.second;  // changes wrt mapt2 : dx->-dx

      // \f$\delta_6=m_4^2-m_5^2\f$ as defined in Vermaseren's paper
      const double tau = t1_ - t2_, r3 = deltas_[3] - t2_, r4 = masses_.w52 - t2_;

      CG_DEBUG_LOOP("LPAIR") << "tau = " << tau << ", r1-4 = " << r1 << ", " << r2 << ", " << r3 << ", " << r4;

      const double b = r3 * r4 - 2. * (t1_ + mB2_) * t2_;
      const double c = t2_ * d6 * d8 + (d6 - d8) * (d6 * mB2_ - d8 * mY2_);

      const double t25 = t2_ - mB2_ - mY2_;

      sa2_ = -0.25 * r4 * r4 + mB2_ * t2_;
      if (sa2_ >= 0.) {
        CG_WARNING("LPAIR") << "sa2_ = " << sa2_ << " >= 0";
        return false;
      }

      const double sl6 = 2. * sqrt(-sa2_);

      g4_ = -r3 * r3 / 4. + t1_ * t2_;
      if (g4_ >= 0.) {
        CG_WARNING("LPAIR") << "g4_ = " << g4_ << " >= 0";
        return false;
      }

      const double sl7 = 2. * sqrt(-g4_);
      const double sl5 = sl6 * sl7;

      double s2p;
      if (fabs((sl5 - b) / sl5) >= 1.) {
        s2p = 0.5 * (sl5 - b) / t2_;
        s2_range.min() = c / (t2_ * s2p);
      } else {  // 8
        s2_range.min() = 0.5 * (-sl5 - b) / t2_;
        s2p = c / (t2_ * s2_range.min());
      }
      // 9
      if (n_opt_ >= 1) {
        const auto s2 = n_opt_ > 1 ? map(u_s2_, s2_range, "s2") : mapla(t1_, mB2_, u_s2_, s2_range);
        s2_ = s2.first;
        ds2 = s2.second;
      }

      const double ap = -0.25 * pow(s2_ + d8, 2) + s2_ * t1_;

      deltas_[0] = 0.25 * (s2_ - s2_range.max()) * (mA2_ != 0. ? (splus - s2_) * mA2_ : ss * t13);
      deltas_[1] = 0.25 * (s2_ - s2_range.min()) * (s2p - s2_) * t2_;

      CG_DEBUG_LOOP("LPAIR") << "\n\t"
                             << "t2       = " << t2_ << "\n\t"
                             << "s2       = " << s2_ << "\n\t"
                             << "s2p      = " << s2p << "\n\t"
                             << "splus    = " << splus << "\n\t"
                             << "s2 range = " << s2_range;

      const double yy4 = cos(theta4_);
      const double dd = deltas_[0] * deltas_[1];
      p12_ = 0.5 * (s_ - mA2_ - mB2_);
      const double st = s2_ - t1_ - mB2_;
      const double delb = (2. * mB2_ * r3 + r4 * st) * (4. * p12_ * t1_ - (t1_ - masses_.w31) * st) / (16. * ap);

      CG_DEBUG_LOOP("LPAIR") << std::scientific << "dd = " << dd << ", "
                             << "dd1 = " << deltas_[0] << ", "
                             << "dd2 = " << deltas_[1] << std::fixed;

      if (dd <= 0.) {
        CG_WARNING("LPAIR:pickin") << "dd = " << dd << " <= 0.";
        return false;
      }

      delta_ = delb - yy4 * st * sqrt(dd) / ap * 0.5;
      s1_ = t2_ + mA2_ + (2. * p12_ * r3 - 4. * delta_) / st;

      if (ap >= 0.) {
        CG_WARNING("LPAIR:pickin") << "ap = " << ap << " >= 0";
        return false;
      }

      jacobian_ = ds2 * dt1 * dt2 * 0.125 * 0.5 / (sl1_ * sqrt(-ap));
      if (jacobian_ == 0.) {
        CG_WARNING("LPAIR:pickin") << "Null Jacobian.\n\t"
                                   << "D(s2)=" << ds2 << ", D(t1)=" << dt1 << ", D(t2)=" << dt2 << ".";
        return false;
      }

      CG_DEBUG_LOOP("LPAIR:pickin") << "ds2=" << ds2 << ", dt1=" << dt1 << ", dt2=" << dt2 << "\n\t"
                                    << "Jacobian=" << std::scientific << jacobian_ << std::fixed;

      gram_ = (1. - yy4 * yy4) * dd / ap;

      p13_ = -0.5 * t13;
      p14_ = 0.5 * (tau + s1_ - mX2_);
      p25_ = -0.5 * t25;

      p1k2_ = 0.5 * (s1_ - t2_ - mA2_);
      p2k1_ = 0.5 * st;

      if (mB2_ != 0.) {
        const double inv_w2 = 1. / mB2_;
        const double sbb = 0.5 * (s_ * (t2_ - masses_.w52) - masses_.w12 * t25) * inv_w2 + mY2_,
                     sdd = 0.5 * sl1_ * sl6 * inv_w2,
                     see = (s_ * (t2_ * (s_ + t25 - mA2_) - mA2_ * masses_.w52) +
                            mY2_ * (mA2_ * mY2_ - masses_.w12 * (t2_ - mA2_))) *
                           inv_w2;
        double s1m, s1p;
        if (sbb * sdd >= 0.) {  // multiplication is more effective than division to check sign+non-null
          s1p = sbb + sdd;
          s1m = see / s1p;
        } else {
          s1m = sbb - sdd;
          s1p = see / s1m;
        }                                                       // 12
        deltas_[2] = -0.25 * mB2_ * (s1p - s1_) * (s1m - s1_);  // 13
      } else {                                                  // 14
        const double s1p = (s_ * (t2_ * (s_ - mY2_ + t2_ - mA2_) - mA2_ * mY2_) + mA2_ * mY2_ * (mA2_ + mY2_ - t2_)) /
                           (t25 * (s_ - masses_.w12));
        deltas_[2] = -0.25 * t25 * (s_ - masses_.w12) * (s1p - s1_);
      }
      // 15
      //const double acc3 = (s1p-s1_)/(s1p+s1_);

      const double ssb = t2_ + 0.5 * mA2_ - r3 * (masses_.w31 - t1_) / t1_, ssd = sl3 * sl7 / t1_,
                   sse = (t2_ - mA2_) * (w4_ - mX2_) +
                         (t2_ - w4_ + masses_.w31) * ((t2_ - mA2_) * mX2_ - (w4_ - mX2_) * mA2_) / t1_;

      double s1pp, s1pm;
      if (ssb / ssd >= 0.) {
        s1pp = ssb + ssd;
        s1pm = sse / s1pp;
      } else {  // 16
        s1pm = ssb - ssd;
        s1pp = sse / s1pm;
      }
      // 17
      deltas_[3] = -0.25 * t1_ * (s1_ - s1pp) * (s1_ - s1pm);
      //const double acc4 = ( s1_-s1pm )/( s1_+s1pm );
      deltas_[4] = deltas_[0] + deltas_[2] +
                   ((p12_ * (t1_ - masses_.w31) * 0.5 - mA2_ * p2k1_) * (p2k1_ * (t2_ - masses_.w52) - mB2_ * r3) -
                    delta_ * (2. * p12_ * p2k1_ - mB2_ * (t1_ - masses_.w31))) /
                       p2k1_;
      if (deltas_[4] < 0.) {
        CG_WARNING("LPAIR") << "dd5 = " << deltas_[4] << " < 0";
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

      const double re = 0.5 / sqs_;
      ep1_ = re * (s_ + masses_.w12);
      ep2_ = re * (s_ - masses_.w12);

      CG_DEBUG_LOOP("LPAIR") << std::scientific << " re = " << re << "\n\t"
                             << "w12 = " << masses_.w12 << std::fixed;
      CG_DEBUG_LOOP("LPAIR") << "Incoming particles' energy = " << ep1_ << ", " << ep2_;

      p_cm_ = re * sl1_;

      de3_ = re * (s2_ - mX2_ + masses_.w12);
      de5_ = re * (s1_ - mY2_ - masses_.w12);

      //----- central two-photon/lepton system

      ec4_ = de3_ + de5_;
      if (ec4_ < mc4_) {
        CG_WARNING("LPAIR") << "ec4_ = " << ec4_ << " < mc4_ = " << mc4_ << "\n\t"
                            << "==> de3 = " << de3_ << ", de5 = " << de5_;
        return false;
      }

      // What if the protons' momenta are not along the z-axis?
      pc4_ = sqrt(ec4_ * ec4_ - mc4_ * mc4_);
      if (pc4_ == 0.) {
        CG_WARNING("LPAIR") << "pzc4 is null and should not be...";
        return false;
      }

      CG_DEBUG_LOOP("LPAIR") << "Central system's energy: E4 = " << ec4_ << "\n\t"
                             << "               momentum: p4 = " << pc4_ << "\n\t"
                             << "         invariant mass: m4 = " << mc4_ << ".";

      pt4_ = sqrt(deltas_[4] / s_) / p_cm_;
      sin_theta4_ = pt4_ / pc4_;

      if (sin_theta4_ > 1.) {
        CG_WARNING("LPAIR") << "st4 = " << sin_theta4_ << " > 1";
        return false;
      }

      cos_theta4_ = sqrt(1. - sin_theta4_ * sin_theta4_);
      if (ep1_ * ec4_ < p14_)
        cos_theta4_ *= -1.;

      al4_ = 1. - cos_theta4_;
      be4_ = 1. + cos_theta4_;

      if (cos_theta4_ < 0.)
        be4_ = sin_theta4_ * sin_theta4_ / al4_;
      else
        al4_ = sin_theta4_ * sin_theta4_ / be4_;

      CG_DEBUG_LOOP("LPAIR") << "cos(theta4) = " << cos_theta4_ << "\t"
                             << "sin(theta4) = " << sin_theta4_ << "\n\t"
                             << "al4 = " << al4_ << ", be4 = " << be4_;

      const double rr = sqrt(-gram_ / s_) / (p_cm_ * pt4_);

      //----- outgoing beam states

      //--- beam 1 -> 3
      const double ep3 = ep1_ - de3_, pp3 = sqrt(ep3 * ep3 - mX2_);
      const double pt3 = sqrt(deltas_[0] / s_) / p_cm_;

      if (pt3 > pp3) {
        CG_WARNING("LPAIR") << "Invalid momentum for outgoing beam 1.";
        return false;
      }
      if (fabs(rr) > pt3) {
        CG_WARNING("LPAIR") << "Invalid momentum balance for outgoing beam 1.";
        return false;
      }

      p3_lab_ = Momentum::fromPThetaPhiE(pp3, -asin(pt3 / pp3), asin(-rr / pt3), ep3);

      CG_DEBUG_LOOP("LPAIR") << "Positive-z beam state:\n\t" << std::scientific << "energy: E3 = " << ep3
                             << ", pt3 = " << pt3 << "\n\t"
                             << "momentum = " << p3_lab_ << ".";

      //--- beam 2 -> 5
      const double ep5 = ep2_ - de5_, pp5 = sqrt(ep5 * ep5 - mY2_);
      const double pt5 = sqrt(deltas_[2] / s_) / p_cm_;

      if (pt5 > pp5) {
        CG_WARNING("LPAIR") << "Invalid momentum for outgoing beam 2.";
        return false;
      }
      if (fabs(rr) > pt5) {
        CG_WARNING("LPAIR") << "Invalid momentum balance for outgoing beam 2.";
        return false;
      }

      p5_lab_ = Momentum::fromPThetaPhiE(pp5, M_PI + asin(pt5 / pp5), asin(rr / pt5), ep5);

      CG_DEBUG_LOOP("LPAIR") << "Negative-z beam state:\n\t" << std::scientific << "energy: E5 = " << ep5
                             << ", pt5 = " << pt5 << "\n\t"
                             << "momentum = " << p5_lab_ << ".";

      //--- mirroring
      const double a1 = p3_lab_.px() - p5_lab_.px();

      CG_DEBUG_LOOP("LPAIR") << "a1 = " << a1;

      if (fabs(pt4_ + p3_lab_.px() + p5_lab_.px()) < fabs(fabs(a1) - pt4_)) {
        CG_DEBUG_LOOP("LPAIR") << "|pt4+pt3*cos(phi3)+pt5*cos(phi5)| < | |a1|-pt4 |\n\t"
                               << "pt4 = " << pt4_ << ".";
        return true;
      }
      if (a1 < 0.)
        p5_lab_.mirrorX();
      else
        p3_lab_.mirrorX();
      return true;
    }

    //---------------------------------------------------------------------------------------------

    double LPAIR::computeWeight() {
      ep1_ = (*event_)(Particle::IncomingBeam1)[0].energy();
      ep2_ = (*event_)(Particle::IncomingBeam2)[0].energy();
      // Mass difference between the first outgoing particle and the first incoming particle
      masses_.w31 = mX2_ - mA2_;
      // Mass difference between the second outgoing particle and the second incoming particle
      masses_.w52 = mY2_ - mB2_;
      // Mass difference between the two incoming particles
      masses_.w12 = mA2_ - mB2_;
      // Mass difference between the central two-photons system and the second outgoing particle

      const double mx = sqrt(mX2_), my = sqrt(mY2_);
      CG_DEBUG_LOOP("LPAIR") << "sqrt(s) = " << sqs_ << " GeV\n\t"
                             << "m^2(X) = " << mX2_ << " GeV^2, m(X) = " << mx << " GeV\n\t"
                             << "m^2(Y) = " << mY2_ << " GeV^2, m(Y) = " << my << " GeV";

      // Maximal energy for the central system is beam-beam CM energy with the outgoing particles' mass energy subtracted (or wmax if specified)
      w_limits_.max() = std::min(pow(sqs_ - mx - my, 2), w_limits_.max());

      // compute the two-photon energy for this point
      mc4_ = sqrt(w4_);

      CG_DEBUG_LOOP("LPAIR") << "Computed value for w4 = " << w4_ << " -> mc4 = " << mc4_;

      if (!orient()) {
        CG_DEBUG_LOOP("LPAIR") << "Orient failed.";
        return 0.;
      }

      if (t1_ > 0.) {
        CG_WARNING("LPAIR") << "t1 = " << t1_ << " > 0";
        return 0.;
      }
      if (t2_ > 0.) {
        CG_WARNING("LPAIR") << "t2 = " << t2_ << " > 0";
        return 0.;
      }

      const double ecm6 = w4_ / (2. * mc4_), pp6cm = sqrt(ecm6 * ecm6 - masses_.Ml2);
      const double alpha1 = (*alphaem_)(std::sqrt(-t1_)), alpha2 = (*alphaem_)(std::sqrt(-t2_));

      jacobian_ *= pp6cm * constb_ * alpha1 * alpha1 * alpha2 * alpha2 / mc4_ / s_;

      // Let the most obscure part of this code begin...

      const double e1mp1 = mA2_ / (ep1_ + p_cm_);
      const double e3mp3 = mX2_ / (p3_lab_.energy() + p3_lab_.p());

      const double al3 = pow(sin(p3_lab_.theta()), 2) / (1. + (p3_lab_.theta()));

      // 2-photon system kinematics ?!
      const double eg = (w4_ + t1_ - t2_) / (2. * mc4_);
      double p_gam = sqrt(eg * eg - t1_);

      /*const double pgx = -p3_lab_.px()*cos_theta4_-sin_theta4_*( de3_-e1mp1 + e3mp3 + p3_lab_.p()*al3 ),
                   pgy = -p3_lab_.py(),
                   pgz = mc4_*de3_/( ec4_+pc4_ )-ec4_*de3_*al4_/mc4_-p3_lab_.px()*ec4_*sin_theta4_/mc4_+ec4_*cos_theta4_/mc4_*( p3_lab_.p()*al3+e3mp3-e1mp1 );*/
      const double gam4 = ec4_ / mc4_;
      const Momentum pg(-p3_lab_.px() * cos_theta4_ - (p3_lab_.p() * al3 + e3mp3 - e1mp1 + de3_) * sin_theta4_,
                        -p3_lab_.py(),
                        -gam4 * p3_lab_.px() * sin_theta4_ + (p3_lab_.p() * al3 + e3mp3 - e1mp1) * gam4 * cos_theta4_ +
                            mc4_ * de3_ / (ec4_ + pc4_) - gam4 * de3_ * al4_);

      CG_DEBUG_LOOP("LPAIR") << "pg = " << pg;

      const double pt_gam = pg.pt(), p_gam_tmp = pg.p();
      if (p_gam_tmp > pt_gam * 0.9 && p_gam_tmp > p_gam)
        p_gam = p_gam_tmp;  //FIXME ???

      // angles for the 2-photon system ?!
      const double cos_phi_gam = pg.px() / pt_gam, sin_phi_gam = pg.py() / pt_gam;
      const double sin_theta_gam = pt_gam / p_gam;

      const int theta_sign = pg.pz() > 0. ? 1 : -1;
      const double cos_theta_gam = theta_sign * sqrt(1. - sin_theta_gam * sin_theta_gam);

      const double amap = 0.5 * (w4_ - t1_ - t2_),
                   bmap = 0.5 * sqrt((pow(w4_ - t1_ - t2_, 2) - 4. * t1_ * t2_) * (1. - 4. * masses_.Ml2 / w4_)),
                   ymap = (amap + bmap) / (amap - bmap), beta = pow(ymap, 2. * x6_ - 1.);
      double xx6 = 0.5 * (1. + amap / bmap * (beta - 1.) / (beta + 1.));
      xx6 = std::max(0., std::min(xx6, 1.));  // xx6 in [0., 1.]

      CG_DEBUG_LOOP("LPAIR") << "amap = " << amap << "\n\t"
                             << "bmap = " << bmap << "\n\t"
                             << "ymap = " << ymap << "\n\t"
                             << "beta = " << beta;

      // 3D rotation of the first outgoing lepton wrt the CM system
      const double theta6cm = acos(1. - 2. * xx6);

      // match the Jacobian
      jacobian_ *= (amap + bmap * cos(theta6cm));
      jacobian_ *= (amap - bmap * cos(theta6cm));
      jacobian_ /= amap;
      jacobian_ /= bmap;
      jacobian_ *= log(ymap);
      if ((kin_.incomingBeams().mode() == mode::Kinematics::ElasticInelastic ||
           kin_.incomingBeams().mode() == mode::Kinematics::InelasticElastic) &&
          symmetrise_)
        jacobian_ *= 1.;
      else
        jacobian_ *= 0.5;

      CG_DEBUG_LOOP("LPAIR") << "Jacobian = " << jacobian_;

      CG_DEBUG_LOOP("LPAIR") << "ctcm6 = " << cos(theta6cm) << "\n\t"
                             << "stcm6 = " << sin(theta6cm);

      // First outgoing lepton's 3-momentum in the centre of mass system
      auto p6cm = Momentum::fromPThetaPhiE(pp6cm, theta6cm, phi6_cm_);

      CG_DEBUG_LOOP("LPAIR") << "p3cm6 = " << p6cm;

      const double h1 = p6cm.pz() * sin_theta_gam + p6cm.px() * cos_theta_gam;
      const double pc6z = p6cm.pz() * cos_theta_gam - p6cm.px() * sin_theta_gam;
      const double pc6x = h1 * cos_phi_gam - p6cm.py() * sin_phi_gam;

      const double qcx = 2. * pc6x, qcz = 2. * pc6z;

      const double el6 = (ec4_ * ecm6 + pc4_ * pc6z) / mc4_;
      const double h2 = (ec4_ * pc6z + pc4_ * ecm6) / mc4_;

      CG_DEBUG_LOOP("LPAIR") << "h1 = " << h1 << ", h2 = " << h2;

      // first outgoing lepton's kinematics
      p6_cm_ = Momentum(+pc6x * cos_theta4_ + h2 * sin_theta4_,
                        p6cm.py() * cos_phi_gam + h1 * sin_phi_gam,
                        -pc6x * sin_theta4_ + h2 * cos_theta4_,
                        el6);

      CG_DEBUG_LOOP("LPAIR") << "p6(cm) = " << p6_cm_;

      const double hq = ec4_ * qcz / mc4_;

      const Momentum qve(+qcx * cos_theta4_ + hq * sin_theta4_,
                         2. * p6_cm_.py(),
                         -qcx * sin_theta4_ + hq * cos_theta4_,
                         +qcz * pc4_ / mc4_  // energy
      );

      // second outgoing lepton's kinematics
      p7_cm_ = Momentum::fromPThetaPhiE(pc4_, acos(cos_theta4_), 0., ec4_) - p6_cm_;

      CG_DEBUG_LOOP("LPAIR") << "Outgoing kinematics\n\t"
                             << " first outgoing lepton: p = " << p6_cm_.p() << ", E = " << p6_cm_.energy() << "\n\t"
                             << "second outgoing lepton: p = " << p7_cm_.p() << ", E = " << p7_cm_.energy();

      q1dq_ = eg * (2. * ecm6 - mc4_) - 2. * p_gam * p6cm.pz();
      q1dq2_ = 0.5 * (w4_ - t1_ - t2_);

      CG_DEBUG_LOOP("LPAIR") << "ecm6 = " << ecm6 << ", mc4 = " << mc4_ << "\n\t"
                             << "eg = " << eg << ", pg = " << p_gam << "\n\t"
                             << "q1dq = " << q1dq_ << ", q1dq2 = " << q1dq2_;

      const double phi3 = p3_lab_.phi(), cos_phi3 = cos(phi3), sin_phi3 = sin(phi3);
      const double phi5 = p5_lab_.phi(), cos_phi5 = cos(phi5), sin_phi5 = sin(phi5);

      bb_ = t1_ * t2_ + (w4_ * pow(sin(theta6cm), 2) + 4. * masses_.Ml2 * pow(cos(theta6cm), 2)) * p_gam * p_gam;

      const double c1 = p3_lab_.pt() * (qve.px() * sin_phi3 - qve.py() * cos_phi3),
                   c2 = p3_lab_.pt() * (qve.pz() * ep1_ - qve.energy() * p_cm_),
                   c3 = (masses_.w31 * ep1_ * ep1_ + 2. * mA2_ * de3_ * ep1_ - mA2_ * de3_ * de3_ +
                         p3_lab_.pt2() * ep1_ * ep1_) /
                        (p3_lab_.energy() * p_cm_ + p3_lab_.pz() * ep1_);

      const double b1 = p5_lab_.pt() * (qve.px() * sin_phi5 - qve.py() * cos_phi5),
                   b2 = p5_lab_.pt() * (qve.pz() * ep2_ + qve.energy() * p_cm_),
                   b3 = (masses_.w52 * ep2_ * ep2_ + 2. * mB2_ * de5_ * ep2_ - mB2_ * de5_ * de5_ +
                         p5_lab_.pt2() * ep2_ * ep2_) /
                        (ep2_ * p5_lab_.pz() - p5_lab_.energy() * p_cm_);

      const double r12 = c2 * sin_phi3 + c3 * qve.py(), r13 = -c2 * cos_phi3 - c3 * qve.px();

      const double r22 = b2 * sin_phi5 + b3 * qve.py(), r23 = -b2 * cos_phi5 - b3 * qve.px();

      epsi_ = p12_ * c1 * b1 + r12 * r22 + r13 * r23;

      g5_ = mA2_ * c1 * c1 + r12 * r12 + r13 * r13;
      g6_ = mB2_ * b1 * b1 + r22 * r22 + r23 * r23;

      const double pt3 = p3_lab_.pt(), pt5 = p5_lab_.pt();
      a5_ = -(qve.px() * cos_phi3 + qve.py() * sin_phi3) * pt3 * p1k2_ -
            (ep1_ * qve.energy() - p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
            (de5_ * qve.pz() + qve.energy() * (p_cm_ + p5_lab_.pz())) * c3;
      a6_ = -(qve.px() * cos_phi5 + qve.py() * sin_phi5) * pt5 * p2k1_ -
            (ep2_ * qve.energy() + p_cm_ * qve.pz()) * (cos_phi3 * cos_phi5 + sin_phi3 * sin_phi5) * pt3 * pt5 +
            (de3_ * qve.pz() - qve.energy() * (p_cm_ - p3_lab_.pz())) * b3;

      CG_DEBUG_LOOP("LPAIR") << "a5 = " << a5_ << "\n\t"
                             << "a6 = " << a6_;

      ////////////////////////////////////////////////////////////////
      // END of GAMGAMLL subroutine in the FORTRAN version
      ////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////
      // INFO from f.f
      ////////////////////////////////////////////////////////////////

      const Momentum cm = p1_lab_ + p2_lab_;
      const double gamma = cm.energy() / sqs_, betgam = cm.pz() / sqs_;

      //--- kinematics computation for both leptons

      CG_DEBUG_LOOP("LPAIR:gmufil") << "unboosted P(l1)=" << p6_cm_ << "\n\t"
                                    << "unboosted P(l2)=" << p7_cm_;

      const double mass_before = (p6_cm_ + p7_cm_).mass();

      p6_cm_.betaGammaBoost(gamma, betgam);
      p7_cm_.betaGammaBoost(gamma, betgam);

      CG_DEBUG_LOOP("LPAIR:gmufil") << "Invariant mass difference from boost:" << (p6_cm_ + p7_cm_).mass() - mass_before
                                    << ".";

      //--- cut on mass of final hadronic system (MX/Y)

      if (kin_.cuts().remnants.mx().valid()) {
        if (kin_.incomingBeams().positive().mode() == Beam::Mode::ProtonInelastic &&
            !kin_.cuts().remnants.mx().contains(mx))
          return 0.;
        if (kin_.incomingBeams().negative().mode() == Beam::Mode::ProtonInelastic &&
            !kin_.cuts().remnants.mx().contains(my))
          return 0.;
      }

      //--- cut on the proton's Q2 (first photon propagator T1)

      if (!kin_.cuts().initial.q2().contains(-t1_))
        return 0.;

      //----- cuts on the individual leptons

      if (kin_.cuts().central.pt_single().valid()) {
        const Limits& pt_limits = kin_.cuts().central.pt_single();
        if (!pt_limits.contains(p6_cm_.pt()) || !pt_limits.contains(p7_cm_.pt()))
          return 0.;
      }

      if (kin_.cuts().central.energy_single().valid()) {
        const Limits& energy_limits = kin_.cuts().central.energy_single();
        if (!energy_limits.contains(p6_cm_.energy()) || !energy_limits.contains(p7_cm_.energy()))
          return 0.;
      }

      if (kin_.cuts().central.eta_single().valid()) {
        const Limits& eta_limits = kin_.cuts().central.eta_single();
        if (!eta_limits.contains(p6_cm_.eta()) || !eta_limits.contains(p7_cm_.eta()))
          return 0.;
      }

      //--- compute the structure functions factors

      jacobian_ *= periPP();

      CG_DEBUG_LOOP("LPAIR:f") << "Jacobian: " << jacobian_;

      //--- compute the event weight using the Jacobian

      return constants::GEVM2_TO_PB * jacobian_;
    }

    //---------------------------------------------------------------------------------------------

    void LPAIR::fillKinematics(bool) {
      const Momentum cm = event_->oneWithRole(Particle::IncomingBeam1).momentum() +
                          event_->oneWithRole(Particle::IncomingBeam2).momentum();

      const double gamma = cm.energy() / sqs_, betgam = cm.pz() / sqs_;

      CG_DEBUG_LOOP("LPAIR:gmufil") << "sqrt(s)=" << sqs_ << " GeV, initial two-proton system: " << cm << "\n\t"
                                    << "gamma=" << gamma << ", betgam=" << betgam;

      auto plab_ip1 = Momentum(0., 0., p_cm_, ep1_).betaGammaBoost(gamma, betgam);
      auto plab_ip2 = Momentum(0., 0., -p_cm_, ep2_).betaGammaBoost(gamma, betgam);

      CG_DEBUG_LOOP("LPAIR:gmufil") << "unboosted PX=" << p3_lab_ << "\n\t"
                                    << "unboosted PY=" << p5_lab_;

      p3_lab_.betaGammaBoost(gamma, betgam);
      p5_lab_.betaGammaBoost(gamma, betgam);

      CG_DEBUG_LOOP("LPAIR:gmufil") << "boosted PX=" << p3_lab_ << "\n\t"
                                    << "boosted PY=" << p5_lab_ << "\n\t"
                                    << "boosted P(l1)=" << p6_cm_ << "\n\t"
                                    << "boosted P(l2)=" << p7_cm_;

      //----- parameterise a random rotation around z-axis
      const short rany = rnd_side_(rnd_gen_) == 1 ? 1 : -1, ransign = rnd_side_(rnd_gen_) == 1 ? 1 : -1;
      const double ranphi = rnd_phi_(rnd_gen_);
      const short ranz = symmetrise_ ? (rnd_side_(rnd_gen_) == 1 ? 1 : -1) : 1;

      Momentum plab_ph1 = plab_ip1 - p3_lab_;
      plab_ph1.rotatePhi(ranphi, rany);

      Momentum plab_ph2 = plab_ip2 - p5_lab_;
      plab_ph2.rotatePhi(ranphi, rany);

      p3_lab_.rotatePhi(ranphi, rany);
      p5_lab_.rotatePhi(ranphi, rany);
      p6_cm_.rotatePhi(ranphi, rany);
      p7_cm_.rotatePhi(ranphi, rany);

      CG_DEBUG_LOOP("LPAIR:gmufil") << "boosted+rotated PX=" << p3_lab_ << "\n\t"
                                    << "boosted+rotated PY=" << p5_lab_ << "\n\t"
                                    << "boosted+rotated P(l1)=" << p6_cm_ << "\n\t"
                                    << "boosted+rotated P(l2)=" << p7_cm_;

      //----- incoming protons
      event_->oneWithRole(Particle::IncomingBeam1).setMomentum(plab_ip1);
      event_->oneWithRole(Particle::IncomingBeam2).setMomentum(plab_ip2);

      //----- first outgoing proton
      auto& op1 = event_->oneWithRole(Particle::OutgoingBeam1);
      p3_lab_.setPz(p3_lab_.pz() * ranz);
      op1.setMomentum(p3_lab_);
      if (kin_.incomingBeams().positive().fragmented()) {
        op1.setStatus(Particle::Status::Unfragmented);  // fragmenting remnants
        op1.setMass(sqrt(mX2_));
      } else
        op1.setStatus(Particle::Status::FinalState);  // stable proton

      //----- second outgoing proton
      auto& op2 = event_->oneWithRole(Particle::OutgoingBeam2);
      p5_lab_.setPz(p5_lab_.pz() * ranz);
      op2.setMomentum(p5_lab_);
      if (kin_.incomingBeams().negative().fragmented()) {
        op2.setStatus(Particle::Status::Unfragmented);  // fragmenting remnants
        op2.setMass(sqrt(mY2_));
      } else
        op2.setStatus(Particle::Status::FinalState);  // stable proton

      //----- first incoming photon
      auto& ph1 = event_->oneWithRole(Particle::Parton1);
      plab_ph1.setPz(plab_ph1.pz() * ranz);
      ph1.setMomentum(plab_ph1);

      //----- second incoming photon
      auto& ph2 = event_->oneWithRole(Particle::Parton2);
      plab_ph2.setPz(plab_ph2.pz() * ranz);
      ph2.setMomentum(plab_ph2);

      auto central_system = (*event_)[Particle::CentralSystem];

      //----- first outgoing lepton
      auto& ol1 = central_system[0].get();
      ol1.setChargeSign(+ransign);
      p6_cm_.setPz(p6_cm_.pz() * ranz);
      ol1.setMomentum(p6_cm_);
      ol1.setStatus(Particle::Status::FinalState);

      //----- second outgoing lepton
      auto& ol2 = central_system[1].get();
      ol2.setChargeSign(-ransign);
      p7_cm_.setPz(p7_cm_.pz() * ranz);
      ol2.setMomentum(p7_cm_);
      ol2.setStatus(Particle::Status::FinalState);

      //----- intermediate two-lepton system
      event_->oneWithRole(Particle::Intermediate).setMomentum(p6_cm_ + p7_cm_);
    }

    //---------------------------------------------------------------------------------------------

    double LPAIR::periPP() const {
      const double qqq = q1dq_ * q1dq_, qdq = 4. * masses_.Ml2 - w4_;
      const double t11 = 64. *
                         (bb_ * (qqq - g4_ - qdq * (t1_ + t2_ + 2. * masses_.Ml2)) -
                          2. * (t1_ + 2. * masses_.Ml2) * (t2_ + 2. * masses_.Ml2) * qqq) *
                         t1_ * t2_;  // magnetic-magnetic
      const double t12 = 128. * (-bb_ * (deltas_[1] + g6_) - 2. * (t1_ + 2. * masses_.Ml2) * (sa2_ * qqq + a6_ * a6_)) *
                         t1_;  // electric-magnetic
      const double t21 = 128. * (-bb_ * (deltas_[3] + g5_) - 2. * (t2_ + 2. * masses_.Ml2) * (sa1_ * qqq + a5_ * a5_)) *
                         t2_;  // magnetic-electric
      const double t22 = 512. * (bb_ * (delta_ * delta_ - gram_) - pow(epsi_ - delta_ * (qdq + q1dq2_), 2) -
                                 sa1_ * a6_ * a6_ - sa2_ * a5_ * a5_ - sa1_ * sa2_ * qqq);  // electric-electric

      auto* ff = kin_.incomingBeams().formFactors();
      auto* sf = kin_.incomingBeams().structureFunctions();
      //--- compute the electric/magnetic form factors for the two
      //    considered parton momenta transfers
      const auto fp1 = kin_.incomingBeams().positive().flux(-t1_, mX2_, ff, sf);
      const auto fp2 = kin_.incomingBeams().negative().flux(-t2_, mY2_, ff, sf);

      const double peripp =
          (fp1.FM * fp2.FM * t11 + fp1.FE * fp2.FM * t21 + fp1.FM * fp2.FE * t12 + fp1.FE * fp2.FE * t22) /
          pow(2. * t1_ * t2_ * bb_, 2);

      CG_DEBUG_LOOP("LPAIR:peripp") << "bb = " << bb_ << ", qqq = " << qqq << ", qdq = " << qdq << "\n\t"
                                    << "t11 = " << t11 << "\t"
                                    << "t12 = " << t12 << "\n\t"
                                    << "t21 = " << t21 << "\t"
                                    << "t22 = " << t22 << "\n\t"
                                    << "=> PeriPP = " << peripp;

      return peripp;
    }

    std::pair<double, double> LPAIR::map(double expo, const Limits& lim, const std::string& var_name_) {
      const double y = lim.max() / lim.min(), out = lim.min() * pow(y, expo), dout = out * log(y);
      CG_DEBUG_LOOP("LPAIR:map") << "Mapping variable \"" << var_name_ << "\" in range (" << lim << ")"
                                 << " (max/min = " << y << ")\n\t"
                                 << "exponent = " << expo << " => "
                                 << "x = " << out << ", dx = " << dout;
      return {out, dout};
    }

    std::pair<double, double> LPAIR::mapla(double y, double z, int u, const Limits& lim) {
      const double xmb = lim.min() - y - z, xpb = lim.max() - y - z;
      const double c = -4. * y * z;
      const double alp = sqrt(xpb * xpb + c), alm = sqrt(xmb * xmb + c);
      const double am = xmb + alm, ap = xpb + alp;
      const double yy = ap / am, zz = pow(yy, u);

      const double out = y + z + 0.5 * (am * zz - c / (am * zz));
      const double ax = sqrt(pow(out - y - z, 2) + c);
      return {out, ax * log(yy)};
    }

    ParametersDescription LPAIR::description() {
      auto desc = Process::description();
      desc.setDescription("γγ → l⁺l¯ (LPAIR)");
      desc.add<int>("nopt", 0).setDescription("Optimised mode? (inherited from LPAIR, by default disabled = 0)");
      desc.add<int>("pair", (int)PDG::muon).setDescription("Lepton pair considered");
      desc.add<bool>("symmetrise", false).setDescription("Symmetrise along z the central system?");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
// register process
REGISTER_PROCESS("lpair", LPAIR)
