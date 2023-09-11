/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Fortran/KTStructures.h"
#include "CepGen/Process/FortranKTProcess.h"
#include "CepGen/Utils/Message.h"

namespace {
  extern "C" {
  extern cepgen::ktblock::Constants constants_;
  extern cepgen::ktblock::GenParameters genparams_;
  extern cepgen::ktblock::KTKinematics ktkin_;
  extern cepgen::ktblock::KinCuts kincuts_;
  extern cepgen::ktblock::EventKinematics evtkin_;
  }
}  // namespace

extern "C" {
/// Print the full list of parameters in the runtime process parameters collection
void cepgen_list_params_() {
  CG_LOG << "\t" << cepgen::ParametersDescription(cepgen::proc::FortranKTProcess::kProcParameters).describe(1);
}

/// Retrieve an integer process parameter from runtime parameters collection
/// \param[in] pname Parameter name string
/// \param[in] def Default parameter value if not found in collection
int cepgen_param_int_(char* pname, int& def) {
  //--- first check if the "integer" is a particle id
  if (cepgen::proc::FortranKTProcess::kProcParameters.has<cepgen::ParticleProperties>(pname))
    return cepgen::proc::FortranKTProcess::kProcParameters.get<cepgen::ParticleProperties>(pname).pdgid;
  if (cepgen::proc::FortranKTProcess::kProcParameters.has<unsigned long long>(pname)) {
    unsigned long long ulong_def = def;
    return cepgen::proc::FortranKTProcess::kProcParameters.get<unsigned long long>(pname, ulong_def);
  }
  //--- if not, proceed with retrieving the integer value
  return cepgen::proc::FortranKTProcess::kProcParameters.get<int>(pname, def);
}

/// Retrieve a double precision floating point process parameter from runtime parameters collection
/// \param[in] pname Parameter name string
/// \param[in] def Default parameter value if not found in collection
double cepgen_param_real_(char* pname, double& def) {
  return cepgen::proc::FortranKTProcess::kProcParameters.get<double>(pname, def);
}
}

namespace cepgen {
  namespace proc {
    ParametersList FortranKTProcess::kProcParameters;  ///< List of parameters to steer the process

    FortranKTProcess::FortranKTProcess(const ParametersList& params, std::function<double(void)> func)
        : KTProcess(params, {PDG::muon, PDG::muon}), func_(func) {
      constants_.m_p = Process::mp_;
      constants_.units = constants::GEVM2_TO_PB;
      constants_.pi = M_PI;
    }

    void FortranKTProcess::preparePhaseSpace() {
      const auto lim_rap = kinematics().cuts().central.rapidity_single.truncate(Limits{-6., 6.});
      const auto lim_pt_diff = kinematics().cuts().central.pt_diff.truncate(Limits{0., 50.});
      const auto lim_phi_diff = kinematics().cuts().central.phi_diff.truncate(Limits{0., 2. * M_PI});

      defineVariable(y1_, Mapping::linear, lim_rap, "First central particle rapidity");
      defineVariable(y2_, Mapping::linear, lim_rap, "Second central particle rapidity");
      defineVariable(
          pt_diff_, Mapping::linear, lim_pt_diff, "Transverse momentum difference between central particles");
      defineVariable(phi_pt_diff_, Mapping::linear, lim_phi_diff, "Central particles azimuthal angle difference");

      //===========================================================================================
      // feed phase space cuts to the common block
      //===========================================================================================

      // export the limits into external variables
      auto save_lim = [](const Limits& lim, int& on, double& min, double& max) {
        on = lim.valid();
        min = max = 0.;
        if (lim.hasMin())
          min = lim.min();
        max = lim.hasMax() ? lim.max() : 9999.999;
      };

      save_lim(kinematics().cuts().central.pt_single, kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max);
      save_lim(kinematics().cuts().central.energy_single, kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max);
      save_lim(kinematics().cuts().central.eta_single, kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max);
      save_lim(kinematics().cuts().central.mass_sum, kincuts_.iinvm, kincuts_.invm_min, kincuts_.invm_max);
      save_lim(kinematics().cuts().central.pt_sum, kincuts_.iptsum, kincuts_.ptsum_min, kincuts_.ptsum_max);
      save_lim(kinematics().cuts().central.rapidity_diff, kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max);

      //===========================================================================================
      // feed run parameters to the common block
      //===========================================================================================

      genparams_.icontri = (int)kinematics().incomingBeams().mode();

      //-------------------------------------------------------------------------------------------
      // incoming beams information
      //-------------------------------------------------------------------------------------------

      //--- positive-z incoming beam
      genparams_.inp1 = kinematics().incomingBeams().positive().momentum().pz();
      //--- check if first incoming beam is a heavy ion
      if (HeavyIon::isHI(kinematics().incomingBeams().positive().pdgId())) {
        const auto in1 = HeavyIon::fromPdgId(kinematics().incomingBeams().positive().pdgId());
        genparams_.a_nuc1 = in1.A;
        genparams_.z_nuc1 = (unsigned short)in1.Z;
        if (genparams_.z_nuc1 > 1) {
          event().oneWithRole(Particle::IncomingBeam1).setPdgId((pdgid_t)in1);
          event().oneWithRole(Particle::OutgoingBeam1).setPdgId((pdgid_t)in1);
        }
      } else
        genparams_.a_nuc1 = genparams_.z_nuc1 = 1;

      //--- negative-z incoming beam
      genparams_.inp2 = kinematics().incomingBeams().negative().momentum().pz();
      //--- check if second incoming beam is a heavy ion
      if (HeavyIon::isHI(kinematics().incomingBeams().negative().pdgId())) {
        const auto in2 = HeavyIon::fromPdgId(kinematics().incomingBeams().negative().pdgId());
        genparams_.a_nuc2 = in2.A;
        genparams_.z_nuc2 = (unsigned short)in2.Z;
        if (genparams_.z_nuc2 > 1) {
          event().oneWithRole(Particle::IncomingBeam2).setPdgId((pdgid_t)in2);
          event().oneWithRole(Particle::OutgoingBeam2).setPdgId((pdgid_t)in2);
        }
      } else
        genparams_.a_nuc2 = genparams_.z_nuc2 = 1;

      //-------------------------------------------------------------------------------------------
      // intermediate partons information
      //-------------------------------------------------------------------------------------------

      //FIXME
      //genparams_.iflux1 = (int)kinematics().incomingBeams().positive().ktFlux();
      //genparams_.iflux2 = (int)kinematics().incomingBeams().negative().ktFlux();
    }

    double FortranKTProcess::computeKTFactorisedMatrixElement() {
      //--- set all kinematics variables for this phase space point
      ktkin_.q1t = qt1_;
      ktkin_.q2t = qt2_;
      ktkin_.phiq1t = phi_qt1_;
      ktkin_.phiq2t = phi_qt2_;
      ktkin_.y1 = y1_;
      ktkin_.y2 = y2_;
      ktkin_.ptdiff = pt_diff_;
      ktkin_.phiptdiff = phi_pt_diff_;
      ktkin_.m_x = mX();
      ktkin_.m_y = mY();

      //--- compute the event weight
      return func_();
    }

    void FortranKTProcess::fillCentralParticlesKinematics() {
      //===========================================================================================
      // outgoing beam remnants
      //===========================================================================================

      pX() = Momentum(evtkin_.px);
      pY() = Momentum(evtkin_.py);
      // express these momenta per nucleon
      pX() *= 1. / genparams_.a_nuc1;
      pY() *= 1. / genparams_.a_nuc2;

      //===========================================================================================
      // intermediate partons
      //===========================================================================================

      q1() = pA() - pX();
      q2() = pB() - pY();
      event().oneWithRole(Particle::Intermediate).setMomentum(q1() + q2());

      //===========================================================================================
      // central system
      //===========================================================================================

      auto oc = event()[Particle::CentralSystem];  // retrieve all references
                                                   // to central system particles
      for (int i = 0; i < evtkin_.nout; ++i) {
        auto& p = oc[i].get();  // retrieve a reference to the specific particle
        p.setPdgId((long)evtkin_.pdg[i]);
        p.setStatus(Particle::Status::FinalState);
        p.setMomentum(Momentum(evtkin_.pc[i]));
      }
    }
  }  // namespace proc
}  // namespace cepgen
