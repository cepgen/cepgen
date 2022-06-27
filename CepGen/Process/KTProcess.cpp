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
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/KTProcess.h"

namespace cepgen {
  namespace proc {
    KTProcess::KTProcess(const ParametersList& params,
                         const std::array<pdgid_t, 2>& partons,
                         const std::vector<pdgid_t>& central)
        : Process(params, true), intermediate_parts_(partons), produced_parts_(central) {}

    KTProcess::KTProcess(const KTProcess& proc)
        : Process(proc),
          log_qt_limits_(proc.log_qt_limits_),
          phi_qt_limits_(proc.phi_qt_limits_),
          mx_limits_(proc.mx_limits_),
          intermediate_parts_(proc.intermediate_parts_),
          produced_parts_(proc.produced_parts_) {}

    void KTProcess::addEventContent() {
      Process::setEventContent(
          {// incoming state
           {Particle::IncomingBeam1, PDG::proton},
           {Particle::IncomingBeam2, PDG::proton},
           {Particle::Parton1, intermediate_parts_[0]},
           {Particle::Parton2, intermediate_parts_[1]}},
          {// outgoing state
           {Particle::OutgoingBeam1, {PDG::proton}},
           {Particle::OutgoingBeam2, {PDG::proton}},
           {Particle::CentralSystem, produced_parts_}});
      setExtraContent();
      CG_DEBUG("KTProcess:addEventContent") << "Addition of:\n\t"
                                            << "Intermediate partons: " << intermediate_parts_ << "\n\t"
                                            << "Produced system: " << produced_parts_ << ".\n\t" << *event_;
    }

    void KTProcess::prepareKinematics() {
      const auto& ib1 = event_->oneWithRole(Particle::IncomingBeam1);
      const auto& ib2 = event_->oneWithRole(Particle::IncomingBeam2);
      auto& p1 = event_->oneWithRole(Particle::Parton1);
      auto& p2 = event_->oneWithRole(Particle::Parton2);

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

      p1.setPdgId(kin_.incomingBeams().positive().daughterId());
      p2.setPdgId(kin_.incomingBeams().negative().daughterId());

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2_ = ib1.mass2();
      mY2_ = ib2.mass2();
      if (kin_.incomingBeams().positive().fragmented())
        defineVariable(
            mX2_, Mapping::square, kin_.cuts().remnants.mx(), {1.07, 1000.}, "Positive z proton remnant squared mass");
      if (kin_.incomingBeams().negative().fragmented())
        defineVariable(
            mY2_, Mapping::square, kin_.cuts().remnants.mx(), {1.07, 1000.}, "Negative z proton remnant squared mass");
    }

    double KTProcess::computeWeight() {
      const auto cent_me = computeKTFactorisedMatrixElement();
      if (cent_me <= 0)
        return 0.;

      //--- compute fluxes according to modelling specified in parameters card

      auto* ff = kin_.incomingBeams().formFactors();
      auto* sf = kin_.incomingBeams().structureFunctions();
      const double q2_1 = qt1_ * qt1_, q2_2 = qt2_ * qt2_;
      const double f1 = kin_.incomingBeams().positive().ktFlux(x1_, q2_1, mX2_, ff, sf) * M_1_PI,
                   f2 = kin_.incomingBeams().negative().ktFlux(x2_, q2_2, mY2_, ff, sf) * M_1_PI;

      return f1 * f2 * cent_me;
    }

    void KTProcess::fillKinematics(bool) {
      fillCentralParticlesKinematics();  // process-dependent!
      fillPrimaryParticlesKinematics();
    }

    void KTProcess::fillPrimaryParticlesKinematics() {
      const auto& ib1 = event_->oneWithRole(Particle::IncomingBeam1);
      const auto& ib2 = event_->oneWithRole(Particle::IncomingBeam2);
      auto& ob1 = event_->oneWithRole(Particle::OutgoingBeam1);
      auto& ob2 = event_->oneWithRole(Particle::OutgoingBeam2);
      auto& p1 = event_->oneWithRole(Particle::Parton1);
      auto& p2 = event_->oneWithRole(Particle::Parton2);
      auto& cm = event_->oneWithRole(Particle::Intermediate);

      // beam systems
      ob1.setMomentum(pX_);
      if (kin_.incomingBeams().positive().fragmented())
        ob1.setMass(sqrt(mX2_));
      ob2.setMomentum(pY_);
      if (kin_.incomingBeams().negative().fragmented())
        ob2.setMass(sqrt(mY2_));

      // parton systems
      p1.setMomentum(ib1.momentum() - pX_, true);
      p2.setMomentum(ib2.momentum() - pY_, true);


      // two-parton system
      cm.setMomentum(p1.momentum() + p2.momentum());
    }

    ParametersDescription KTProcess::description() {
      auto desc = Process::description();
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
