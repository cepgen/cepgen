#include "CepGen/Processes/KTProcess.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace proc {
    KTProcess::KTProcess(const ParametersList& params,
                         const std::array<pdgid_t, 2>& partons,
                         const std::vector<pdgid_t>& central)
        : Process(params, true),
          qt1_(0.),
          phi_qt1_(0.),
          qt2_(0.),
          phi_qt2_(0.),
          kIntermediateParts(partons),
          kProducedParts(central) {}

    void KTProcess::addEventContent() {
      Process::setEventContent(
          {// incoming state
           {Particle::IncomingBeam1, PDG::proton},
           {Particle::IncomingBeam2, PDG::proton},
           {Particle::Parton1, kIntermediateParts[0]},
           {Particle::Parton2, kIntermediateParts[1]}},
          {// outgoing state
           {Particle::OutgoingBeam1, {PDG::proton}},
           {Particle::OutgoingBeam2, {PDG::proton}},
           {Particle::CentralSystem, kProducedParts}});
      setExtraContent();
    }

    void KTProcess::prepareKinematics() {
      //============================================================================================
      // try to extrapolate unintegrated fluxes from kinematics mode
      //============================================================================================

      //----- ensure the first incoming flux is compatible with the kinematics mode
      const KTFlux& flux1 = (KTFlux)kin_.incomingBeams().positive().kt_flux;
      const HeavyIon hi1(kin_.incomingBeams().positive().pdg);
      switch (kin_.incomingBeams().positive().mode) {
        case mode::Beam::ProtonElastic:
          if (flux1 != KTFlux::P_Photon_Elastic && flux1 != KTFlux::P_Photon_Elastic_Budnev &&
              flux1 != KTFlux::HI_Photon_Elastic && flux1 != KTFlux::P_Gluon_KMR) {
            kin_.incomingBeams().positive().kt_flux = hi1 ? KTFlux::HI_Photon_Elastic : KTFlux::P_Photon_Elastic_Budnev;
            CG_DEBUG("KTProcess:kinematics") << "KT flux for positive-z incoming parton set to \""
                                             << kin_.incomingBeams().positive().kt_flux << "\".";
          }
          break;
        case mode::Beam::ProtonInelastic:
          if (flux1 != KTFlux::P_Photon_Inelastic && flux1 != KTFlux::P_Photon_Inelastic_Budnev) {
            if (hi1)
              throw CG_FATAL("KTProcess:kinematics") << "Inelastic photon emission from HI not yet supported!";
            kin_.incomingBeams().positive().kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
            CG_INFO("KTProcess:kinematics") << "KT flux for positive-z incoming parton set to \""
                                            << kin_.incomingBeams().positive().kt_flux << "\".";
          }
          break;
        default:
          throw CG_FATAL("KTProcess:kinematics")
              << "Invalid positive-z beam mode for KT process: " << kin_.incomingBeams().positive().mode << "!";
      }

      //----- ensure the second incoming flux is compatible with the kinematics mode
      const KTFlux& flux2 = (KTFlux)kin_.incomingBeams().negative().kt_flux;
      const HeavyIon hi2(kin_.incomingBeams().negative().pdg);
      switch (kin_.incomingBeams().negative().mode) {
        case mode::Beam::ProtonElastic:
          if (flux2 != KTFlux::P_Photon_Elastic && flux2 != KTFlux::P_Photon_Elastic_Budnev &&
              flux2 != KTFlux::HI_Photon_Elastic && flux2 != KTFlux::P_Gluon_KMR) {
            kin_.incomingBeams().negative().kt_flux = hi2 ? KTFlux::HI_Photon_Elastic : KTFlux::P_Photon_Elastic_Budnev;
            CG_DEBUG("KTProcess:kinematics") << "KT flux for negative-z incoming parton set to \""
                                             << kin_.incomingBeams().negative().kt_flux << "\".";
          }
          break;
        case mode::Beam::ProtonInelastic:
          if (flux2 != KTFlux::P_Photon_Inelastic && flux2 != KTFlux::P_Photon_Inelastic_Budnev) {
            if (hi2)
              throw CG_FATAL("KTProcess:kinematics") << "Inelastic photon emission from HI not yet supported!";
            kin_.incomingBeams().negative().kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
            CG_INFO("KTProcess:kinematics") << "KT flux for negative-z incoming parton set to \""
                                            << kin_.incomingBeams().negative().kt_flux << "\".";
          }
          break;
        default:
          throw CG_FATAL("KTProcess:kinematics")
              << "Invalid negative-z beam mode for KT process: " << kin_.incomingBeams().negative().mode << "!";
      }

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      defineVariable(
          qt1_, Mapping::exponential, kin_.cuts().initial.qt(), {1.e-10, 500.}, "First incoming parton virtuality");
      defineVariable(
          qt2_, Mapping::exponential, kin_.cuts().initial.qt(), {1.e-10, 500.}, "Second incoming parton virtuality");
      defineVariable(phi_qt1_,
                     Mapping::linear,
                     kin_.cuts().initial.phi_qt(),
                     {0., 2. * M_PI},
                     "First incoming parton azimuthal angle");
      defineVariable(phi_qt2_,
                     Mapping::linear,
                     kin_.cuts().initial.phi_qt(),
                     {0., 2. * M_PI},
                     "Second incoming parton azimuthal angle");

      //============================================================================================
      // register the incoming partons
      //============================================================================================

      switch (kin_.incomingBeams().positive().kt_flux) {
        case KTFlux::P_Gluon_KMR:
          event_->oneWithRole(Particle::Parton1).setPdgId((pdgid_t)PDG::gluon);
          break;
        case KTFlux::P_Photon_Elastic:
        case KTFlux::P_Photon_Elastic_Budnev:
        case KTFlux::P_Photon_Inelastic:
        case KTFlux::P_Photon_Inelastic_Budnev:
        case KTFlux::HI_Photon_Elastic:
          event_->oneWithRole(Particle::Parton1).setPdgId((pdgid_t)PDG::photon);
          break;
        case KTFlux::invalid:
        default:
          throw CG_FATAL("KTProcess:kinematics")
              << "Invalid flux for 2nd incoming parton: " << kin_.incomingBeams().positive().kt_flux << "!";
      }
      switch (kin_.incomingBeams().negative().kt_flux) {
        case KTFlux::P_Gluon_KMR:
          event_->oneWithRole(Particle::Parton2).setPdgId((pdgid_t)PDG::gluon);
          break;
        case KTFlux::P_Photon_Elastic:
        case KTFlux::P_Photon_Elastic_Budnev:
        case KTFlux::P_Photon_Inelastic:
        case KTFlux::P_Photon_Inelastic_Budnev:
        case KTFlux::HI_Photon_Elastic:
          event_->oneWithRole(Particle::Parton2).setPdgId((pdgid_t)PDG::photon);
          break;
        case KTFlux::invalid:
        default:
          throw CG_FATAL("KTProcess:kinematics")
              << "Invalid flux for 2nd incoming parton: " << kin_.incomingBeams().negative().kt_flux << "!";
      }

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2_ = event_->oneWithRole(Particle::IncomingBeam1).mass2();
      mY2_ = event_->oneWithRole(Particle::IncomingBeam2).mass2();
      if (kin_.incomingBeams().positive().mode == mode::Beam::ProtonInelastic)
        defineVariable(
            mX2_, Mapping::square, kin_.cuts().remnants.mx(), {1.07, 1000.}, "Positive z proton remnant squared mass");
      if (kin_.incomingBeams().negative().mode == mode::Beam::ProtonInelastic)
        defineVariable(
            mY2_, Mapping::square, kin_.cuts().remnants.mx(), {1.07, 1000.}, "Negative z proton remnant squared mass");
    }

    double KTProcess::computeWeight() { return std::max(0., computeKTFactorisedMatrixElement()); }

    void KTProcess::fillKinematics(bool) {
      fillCentralParticlesKinematics();  // process-dependent!
      fillPrimaryParticlesKinematics();
    }

    void KTProcess::fillPrimaryParticlesKinematics() {
      //============================================================================================
      //     outgoing protons
      //============================================================================================

      Particle& op1 = event_->oneWithRole(Particle::OutgoingBeam1);
      op1.setMomentum(pX_);
      if (kin_.incomingBeams().positive().mode == mode::Beam::ProtonElastic)
        op1.setStatus(Particle::Status::FinalState);
      else if (kin_.incomingBeams().positive().mode == mode::Beam::ProtonInelastic)
        op1.setStatus(Particle::Status::Unfragmented).setMass(sqrt(mX2_));
      else
        throw CG_FATAL("KTProcess") << "This kT factorisation process is intended for p-on-p collisions! Aborting.";

      Particle& op2 = event_->oneWithRole(Particle::OutgoingBeam2);
      op2.setMomentum(pY_);
      if (kin_.incomingBeams().negative().mode == mode::Beam::ProtonElastic)
        op2.setStatus(Particle::Status::FinalState);
      else if (kin_.incomingBeams().negative().mode == mode::Beam::ProtonInelastic)
        op2.setStatus(Particle::Status::Unfragmented).setMass(sqrt(mY2_));
      else
        throw CG_FATAL("KTProcess") << "This kT factorisation process is intended for p-on-p collisions! Aborting.";

      //============================================================================================
      //     incoming partons (photons, pomerons, ...)
      //============================================================================================

      Particle& g1 = event_->oneWithRole(Particle::Parton1);
      g1.setMomentum(event_->oneWithRole(Particle::IncomingBeam1).momentum() - pX_, true);

      Particle& g2 = event_->oneWithRole(Particle::Parton2);
      g2.setMomentum(event_->oneWithRole(Particle::IncomingBeam2).momentum() - pY_, true);

      //============================================================================================
      //     two-parton system
      //============================================================================================

      event_->oneWithRole(Particle::Intermediate).setMomentum(g1.momentum() + g2.momentum());
    }
  }  // namespace proc
}  // namespace cepgen
