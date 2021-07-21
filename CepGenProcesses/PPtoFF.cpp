#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/AlphaS.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/Process2to4.h"

namespace cepgen {
  namespace proc {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process using \f$k_{\rm T}\f$-factorization approach
    class PPtoFF : public Process2to4 {
    public:
      PPtoFF(const ParametersList& params = ParametersList());
      ProcessPtr clone() const override { return ProcessPtr(new PPtoFF(*this)); }
      static std::string description() { return "ɣɣ → f⁺f¯ (kt-factor.)"; }

    private:
      void prepareProcessKinematics() override;
      double computeCentralMatrixElement() const override;

      /// Rapidity range for the outgoing fermions
      double onShellME() const;
      double offShellME() const;

      const enum class Mode { onShell = 0, offShell = 1, offShellLegacy = 2 } method_;

      ParametersList alphas_params_;
      std::shared_ptr<AlphaS> alphas_;

      bool gluon1_, gluon2_;
      double prefactor_;

      //--- parameters for off-shell matrix element
      unsigned short p_mat1_, p_mat2_;
      unsigned short p_term_ll_, p_term_lt_, p_term_tt1_, p_term_tt2_;

      double mf2_;
      short qf3_;
      unsigned short colf_;
    };

    PPtoFF::PPtoFF(const ParametersList& params)
        : Process2to4(params, {PDG::photon, PDG::photon}, params.get<ParticleProperties>("pair").pdgid),
          method_((Mode)params.get<int>("method", (int)Mode::offShell)),
          alphas_params_(params.get<ParametersList>("alphaS", ParametersList().setName<std::string>("pegasus"))),
          gluon1_(false),
          gluon2_(false),
          prefactor_(1.),
          p_mat1_(0),
          p_mat2_(0),
          p_term_ll_(0),
          p_term_lt_(0),
          p_term_tt1_(0),
          p_term_tt2_(0),
          mf2_(0.),
          qf3_(0),
          colf_(0) {
      if (method_ == Mode::offShell || method_ == Mode::offShellLegacy) {  // off-shell matrix element
        const auto& ofp = params.get<ParametersList>("offShellParameters");
        p_mat1_ = ofp.get<int>("mat1", method_ == Mode::offShell ? 1 : 2);
        p_mat2_ = ofp.get<int>("mat2", method_ == Mode::offShell ? 1 : 0);
        p_term_ll_ = ofp.get<int>("termLL", 1);
        p_term_lt_ = ofp.get<int>("termLT", 1);
        p_term_tt1_ = ofp.get<int>("termTT", 1);
        p_term_tt2_ = ofp.get<int>("termtt", 1);
      }
    }

    void PPtoFF::prepareProcessKinematics() {
      if (!cs_prop_.fermion || cs_prop_.charge == 0.)
        throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion pair selected: " << cs_prop_.description << " ("
                                         << (int)cs_prop_.pdgid << ")!";

      mf2_ = cs_prop_.mass * cs_prop_.mass;
      qf3_ = cs_prop_.charge;
      colf_ = cs_prop_.colours;
      prefactor_ = 1.;

      CG_DEBUG("PPtoFF:prepare") << "Produced particles: " << cs_prop_.description << " ("
                                 << "mass = " << cs_prop_.mass << " GeV, "
                                 << "charge = " << std::setprecision(2) << qf3_ / 3. << " e)\n\t"
                                 << "matrix element computation method: " << (int)method_ << ".";

      if (!kin_.cuts.central.pt_diff().valid())
        kin_.cuts.central.pt_diff() = {0., 50.};  // tighter cut for fermions

      CG_DEBUG("PPtoFF:prepare") << "Incoming state:\n\t"
                                 << "mp(1/2) = " << sqrt(mA2_) << "/" << sqrt(mB2_) << ".";

      switch (event_->oneWithRole(Particle::Parton1).pdgId()) {
        case PDG::gluon:
          gluon1_ = true;
          prefactor_ *= 4. * M_PI;
          break;
        case PDG::photon:
          prefactor_ *= constants::G_EM_SQ * pow(qf3_, 2) / 9.;
          break;
        default:
          throw CG_FATAL("PPtoFF:prepare") << "Only photon & gluon partons are supported!";
      }
      switch (event_->oneWithRole(Particle::Parton2).pdgId()) {
        case PDG::gluon:
          gluon2_ = true;
          prefactor_ *= 4. * M_PI;
          break;
        case PDG::photon:
          prefactor_ *= constants::G_EM_SQ * pow(qf3_, 2) / 9.;
          break;
        default:
          throw CG_FATAL("PPtoFF:prepare") << "Only photon & gluon partons are supported!";
      }
      if (gluon1_ || gluon2_)
        // at least one gluon; need to initialise the alpha(s) evolution algorithm
        alphas_ = AlphaSFactory::get().build(alphas_params_);
    }

    double PPtoFF::computeCentralMatrixElement() const {
      double mat_el;
      switch (method_) {
        case Mode::onShell:
          mat_el = onShellME();
          break;
        case Mode::offShell:
        case Mode::offShellLegacy:
          mat_el = offShellME();
          break;
        default:
          throw CG_FATAL("PPtoFF") << "Invalid ME calculation method (" << (int)method_ << ")!";
      }
      CG_DEBUG_LOOP("PPtoFF:ME") << "prefactor: " << prefactor_ << "\n\t"
                                 << "matrix element: " << mat_el << ".";

      return mat_el;
    }

    double PPtoFF::onShellME() const {
      if (gluon1_ || gluon2_)
        throw CG_FATAL("PPtoFF:onShell") << "On-shell matrix element only compatible with photon-photon mode!";

      const double s_hat = shat(), t_hat = that(), u_hat = uhat();
      CG_DEBUG_LOOP("PPtoFF:onShell") << "shat: " << s_hat << ", that: " << t_hat << ", uhat: " << u_hat << ".";

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
      return -2. * out / (pow((mf2_ - t_hat) * (mf2_ - u_hat), 2));
    }

    double PPtoFF::offShellME() const {
      const double alpha1 = amt1_ / sqs_ * exp(y_c1_), beta1 = amt1_ / sqs_ * exp(-y_c1_);
      const double alpha2 = amt2_ / sqs_ * exp(y_c2_), beta2 = amt2_ / sqs_ * exp(-y_c2_);
      const double x1 = alpha1 + alpha2, x2 = beta1 + beta2;
      const double z1p = alpha1 / x1, z1m = alpha2 / x1, z1 = z1p * z1m;
      const double z2p = beta1 / x2, z2m = beta2 / x2, z2 = z2p * z2m;

      CG_DEBUG_LOOP("2to4:zeta") << "amt(1/2) = " << amt1_ << " / " << amt2_ << "\n\t"
                                 << "z(1/2)p = " << z1p << " / " << z2p << ", z1 = " << z1 << "\n\t"
                                 << "z(1/2)m = " << z1m << " / " << z2m << ", z2 = " << z2 << ".";

      //--- positive-z photon kinematics
      const Momentum ak1 = (z1m * p_c1_ - z1p * p_c2_).setPz(0.);
      const Momentum ph_p1 = ak1 + z1p * q2_, ph_m1 = ak1 - z1m * q2_;
      const double t1abs = (q1_.pt2() + x1 * (mX2_ - mA2_) + x1 * x1 * mA2_) / (1. - x1);
      const double eps12 = mf2_ + z1 * t1abs;
      const double kp1 = 1. / (ph_p1.pt2() + eps12);
      const double km1 = 1. / (ph_m1.pt2() + eps12);

      const Momentum phi1 = (kp1 * ph_p1 - km1 * ph_m1).setPz(0.).setEnergy(kp1 - km1);
      const double dot1 = phi1.threeProduct(q1_) / qt1_;
      const double cross1 = phi1.crossProduct(q1_) / qt1_;

      //--- negative-z photon kinematics
      const Momentum ak2 = (z2m * p_c1_ - z2p * p_c2_).setPz(0.);
      const Momentum ph_p2 = ak2 + z2p * q1_, ph_m2 = ak2 - z2m * q1_;
      const double t2abs = (q2_.pt2() + x2 * (mY2_ - mB2_) + x2 * x2 * mB2_) / (1. - x2);
      const double eps22 = mf2_ + z2 * t2abs;
      const double kp2 = 1. / (ph_p2.pt2() + eps22);
      const double km2 = 1. / (ph_m2.pt2() + eps22);

      const Momentum phi2 = (kp2 * ph_p2 - km2 * ph_m2).setPz(0.).setEnergy(kp2 - km2);
      const double dot2 = phi2.threeProduct(q2_) / qt2_;
      const double cross2 = phi2.crossProduct(q2_) / qt2_;

      CG_DEBUG_LOOP("PPtoFF:offShell") << "Photon kinematics:\n\t"
                                       << "q1 = " << q1_ << ", q1t = " << qt1_ << ", q1t2 = " << q1_.pt2() << "\n\t"
                                       << "q2 = " << q2_ << ", q2t = " << qt2_ << ", q2t2 = " << q2_.pt2() << "\n\t"
                                       << "phi1 = " << phi1 << ", phi2 = " << phi2 << "\n\t"
                                       << "t1abs = " << t1abs << ", t2abs = " << t2abs << "\n\t"
                                       << "x1 = " << x1 << ", x2 = " << x2 << "\n\t"
                                       << "(dot):   " << dot1 << " / " << dot2 << "\n\t"
                                       << "(cross): " << cross1 << " / " << cross2 << ".";

      const double aux2_1 = p_term_ll_ * (mf2_ + 4. * z1 * z1 * t1abs) * phi1.energy2() +
                            p_term_tt1_ * ((z1p * z1p + z1m * z1m) * (dot1 * dot1 + cross1 * cross1)) +
                            p_term_tt2_ * (cross1 * cross1 - dot1 * dot1) -
                            p_term_lt_ * 4. * z1 * (z1p - z1m) * phi1.energy() * qt1_ * dot1;

      const double aux2_2 = p_term_ll_ * (mf2_ + 4. * z2 * z2 * t2abs) * phi2.energy2() +
                            p_term_tt1_ * ((z2p * z2p + z2m * z2m) * (dot2 * dot2 + cross2 * cross2)) +
                            p_term_tt2_ * (cross2 * cross2 - dot2 * dot2) -
                            p_term_lt_ * 4. * z2 * (z2p - z2m) * phi2.energy() * qt2_ * dot2;

      //=================================================================
      //     convention of matrix element as in our kt-factorization
      //     for heavy flavours
      //=================================================================

      // Marta's version
      double amat2_1 = 2. * aux2_1 * z1, amat2_2 = 2. * aux2_2 * z2;
      const double inv_q1t2 = 1. / q1_.pt2(), inv_q2t2 = 1. / q2_.pt2();
      if (method_ == Mode::offShell) {
        amat2_1 *= inv_q2t2;
        amat2_2 *= inv_q1t2;
      } else if (method_ == Mode::offShellLegacy) {
        amat2_1 *= (t1abs * inv_q1t2 * inv_q2t2) * (t2abs * inv_q2t2);
        amat2_2 *= (t2abs * inv_q1t2 * inv_q2t2);
      }

      //=================================================================
      //     symmetrization
      //=================================================================

      double amat2 = 0.5 * prefactor_ * pow(x1 * x2 * s_, 2) * (p_mat1_ * amat2_1 + p_mat2_ * amat2_2);

      const double tmax = pow(std::max(amt1_, amt2_), 2);
      if (gluon1_)
        amat2 *= 0.5 * (*alphas_)(sqrt(std::max(eps12, tmax)));
      if (gluon2_)
        amat2 *= 0.5 * (*alphas_)(sqrt(std::max(eps22, tmax)));

      CG_DEBUG_LOOP("PPtoFF:offShell") << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
                                       << "z(1/2) = " << z1 << " / " << z2 << "\n\t"
                                       << "q(1/2)t2 = " << q1_.pt2() << " / " << q2_.pt2() << "\n\t"
                                       << "amat2(1/2) = " << amat2_1 << " / " << amat2_2 << "\n\t"
                                       << "amat2 = " << amat2 << ".";

      return amat2;
    }
  }  // namespace proc
}  // namespace cepgen
// register process
REGISTER_PROCESS("pptoff", PPtoFF)
