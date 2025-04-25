/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Fortran/KTStructures.h"
#include "CepGen/Process/FortranFactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace cepgen::proc;

extern "C" {
extern ktblock::Constants constants_;
extern ktblock::GenParameters genparams_;
extern ktblock::KTKinematics ktkin_;
extern ktblock::KinCuts kincuts_;
extern ktblock::EventKinematics evtkin_;

/// Print the full list of parameters in the runtime process parameters collection
void cepgen_list_params_() {
  CG_LOG << "\t" << ParametersDescription(FortranFactorisedProcess::kProcParameters).describe(1);
}

/// Retrieve an integer process parameter from runtime parameters collection
/// \param[in] pname Parameter name string
/// \param[in] def Default parameter value if not found in collection
int cepgen_param_int_(char* pname, int& def) {
  if (const std::string parameter_name(pname); !parameter_name.empty()) {
    // check if the "integer" is a particle id
    if (FortranFactorisedProcess::kProcParameters.has<ParticleProperties>(parameter_name))
      return FortranFactorisedProcess::kProcParameters.get<ParticleProperties>(pname).pdgid;
    if (FortranFactorisedProcess::kProcParameters.has<unsigned long long>(parameter_name))
      return FortranFactorisedProcess::kProcParameters.get<unsigned long long>(parameter_name,
                                                                               static_cast<unsigned long long>(def));
    // if not, proceed with retrieving the integer value
    return FortranFactorisedProcess::kProcParameters.get<int>(parameter_name, def);
  }
  return def;
}

/// Retrieve a double precision floating point process parameter from runtime parameters collection
/// \param[in] pname Parameter name string
/// \param[in] def Default parameter value if not found in collection
double cepgen_param_real_(char* pname, double& def) {
  if (const std::string parameter_name(pname); !parameter_name.empty())
    return FortranFactorisedProcess::kProcParameters.get<double>(parameter_name, def);
  return def;
}
}

auto make_pdgids_pair = [](pdgid_t pair) {
  return spdgids_t{static_cast<spdgid_t>(pair), -static_cast<spdgid_t>(pair)};
};

ParametersList FortranFactorisedProcess::kProcParameters;  ///< List of parameters to steer the process

FortranFactorisedProcess::FortranFactorisedProcess(const ParametersList& params,
                                                   const std::function<double(void)>& func)
    : FactorisedProcess(params), func_(func) {
  constants_.m_p = mp_;
  constants_.units = constants::GEVM2_TO_PB;
  constants_.pi = M_PI;
  if (!func_)
    throw CG_FATAL("FortranFactorisedProcess")
        << "Invalid (Fortran) function given to the factorised process interface.";
  if (phase_space_generator_->central().empty())  // event content was not yet initialised
    setCentral(make_pdgids_pair(steer<ParticleProperties>("pair").pdgid));
}

void FortranFactorisedProcess::prepareFactorisedPhaseSpace() {
  const auto lim_rap = kinematics().cuts().central.rapidity_single.truncate(Limits{-6., 6.});
  defineVariable(m_y1_, Mapping::linear, lim_rap, "y1", "First central particle rapidity");
  defineVariable(m_y2_, Mapping::linear, lim_rap, "y2", "Second central particle rapidity");

  const auto lim_pt_diff = kinematics().cuts().central.pt_diff.truncate(Limits{0., 50.});
  defineVariable(
      m_pt_diff_, Mapping::linear, lim_pt_diff, "pt_diff", "Central particles transverse momentum difference");

  const auto lim_phi_diff = kinematics().cuts().central.phi_diff.truncate(Limits{0., 2. * M_PI});
  defineVariable(
      m_phi_pt_diff_, Mapping::linear, lim_phi_diff, "phi_diff", "Central particles azimuthal angle difference");

  // feed phase space cuts to the common block: export the limits into external variables
  auto save_lim = [](const Limits& lim, int& on, double& min, double& max) {
    on = lim.valid();
    min = lim.hasMin() ? lim.min() : -9999.999;
    max = lim.hasMax() ? lim.max() : +9999.999;
  };
  save_lim(kinematics().cuts().central.pt_single, kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max);
  save_lim(kinematics().cuts().central.energy_single, kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max);
  save_lim(kinematics().cuts().central.eta_single, kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max);
  save_lim(kinematics().cuts().central.mass_sum, kincuts_.iinvm, kincuts_.invm_min, kincuts_.invm_max);
  save_lim(kinematics().cuts().central.pt_sum, kincuts_.iptsum, kincuts_.ptsum_min, kincuts_.ptsum_max);
  save_lim(kinematics().cuts().central.rapidity_diff, kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max);

  // feed run parameters to the common block
  genparams_.icontri = static_cast<int>(kinematics().incomingBeams().mode());

  // incoming beams information
  genparams_.inp1 = kinematics().incomingBeams().positive().momentum().pz();  // positive-z incoming beam
  if (HeavyIon::isHI(
          kinematics().incomingBeams().positive().integerPdgId())) {  // check if first incoming beam is a heavy ion
    const auto in1 = HeavyIon::fromPdgId(kinematics().incomingBeams().positive().integerPdgId());
    genparams_.a_nuc1 = in1.A;
    genparams_.z_nuc1 = static_cast<unsigned short>(in1.Z);
    if (genparams_.z_nuc1 > 1) {
      event().oneWithRole(Particle::Role::IncomingBeam1).setPdgId(in1);
      event().oneWithRole(Particle::Role::OutgoingBeam1).setPdgId(in1);
    }
  } else
    genparams_.a_nuc1 = genparams_.z_nuc1 = 1;
  genparams_.inp2 = kinematics().incomingBeams().negative().momentum().pz();  // negative-z incoming beam
  if (HeavyIon::isHI(
          kinematics().incomingBeams().negative().integerPdgId())) {  // check if second incoming beam is a heavy ion
    const auto in2 = HeavyIon::fromPdgId(kinematics().incomingBeams().negative().integerPdgId());
    genparams_.a_nuc2 = in2.A;
    genparams_.z_nuc2 = static_cast<unsigned short>(in2.Z);
    if (genparams_.z_nuc2 > 1) {
      event().oneWithRole(Particle::Role::IncomingBeam2).setPdgId(in2);
      event().oneWithRole(Particle::Role::OutgoingBeam2).setPdgId(in2);
    }
  } else
    genparams_.a_nuc2 = genparams_.z_nuc2 = 1;

  // intermediate partons information
  genparams_.iflux1 = static_cast<int>(phase_space_generator_->partons().at(0));
  genparams_.iflux2 = static_cast<int>(phase_space_generator_->partons().at(1));
}

double FortranFactorisedProcess::computeFactorisedMatrixElement() {
  //===========================================================================================
  // set all kinematics variables generated by CepGen for this phase space point
  //===========================================================================================

  // intermediate partons
  ktkin_.q1t = q1().p();
  ktkin_.q2t = q2().p();
  ktkin_.phiq1t = q1().phi();
  ktkin_.phiq2t = q2().phi();
  // central system
  ktkin_.y1 = m_y1_;
  ktkin_.y2 = m_y2_;
  ktkin_.ptdiff = m_pt_diff_;
  ktkin_.phiptdiff = m_phi_pt_diff_;
  // outgoing beam remnants
  ktkin_.m_x = mX();
  ktkin_.m_y = mY();
  const auto weight = func_();  // compute the event weight
  if (!utils::positive(weight))
    return 0.;

  //===========================================================================================
  // after weight was computed (and validated), start retrieving the kinematics information
  // from function-filled Fortran common blocks content
  //===========================================================================================

  // outgoing beam remnants (expressed per nucleon)
  pX() = Momentum(evtkin_.px) / genparams_.a_nuc1;
  pY() = Momentum(evtkin_.py) / genparams_.a_nuc2;
  // intermediate partons (including collinear/transverse components)
  q1() = pA() - pX();
  q2() = pB() - pY();
  event().oneWithRole(Particle::Role::Intermediate).setMomentum(q1() + q2());
  // central system
  auto oc = event()[Particle::Role::CentralSystem];  // retrieve all references
                                                     // to central system particles
  for (int i = 0; i < evtkin_.nout; ++i) {
    auto& p = oc[i].get();  // retrieve a reference to the specific particle
    p.setIntegerPdgId(evtkin_.pdg[i]);
    p.setStatus(Particle::Status::FinalState);
    p.setMomentum(Momentum(evtkin_.pc[i]));
  }

  return weight;
}
