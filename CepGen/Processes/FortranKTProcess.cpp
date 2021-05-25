#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/Processes/Fortran/KTStructures.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"

extern "C" {
extern cepgen::ktblock::Constants constants_;
extern cepgen::ktblock::Parameters genparams_;
extern cepgen::ktblock::KTKinematics ktkin_;
extern cepgen::ktblock::Cuts kincuts_;
extern cepgen::ktblock::Event evtkin_;

void cepgen_list_params_() { CG_LOG("cepgen_list_params") << "\t" << cepgen::proc::FortranKTProcess::kProcParameters; }

int cepgen_param_int_(char* pname, int& def) {
  //--- first check if the "integer" is a particle id
  if (cepgen::proc::FortranKTProcess::kProcParameters.has<cepgen::ParticleProperties>(pname))
    return cepgen::proc::FortranKTProcess::kProcParameters.get<cepgen::ParticleProperties>(pname).pdgid;
  //--- if not, proceed with retrieving the integer value
  return cepgen::proc::FortranKTProcess::kProcParameters.get<int>(pname, def);
}

double cepgen_param_real_(char* pname, double& def) {
  return cepgen::proc::FortranKTProcess::kProcParameters.get<double>(pname, def);
}
}

namespace cepgen {
  namespace proc {
    ParametersList FortranKTProcess::kProcParameters;  ///< List of parameters to steer the process

    FortranKTProcess::FortranKTProcess(const ParametersList& params, std::function<double(void)> func)
        : KTProcess(params, {{PDG::photon, PDG::photon}}, {PDG::muon, PDG::muon}), func_(func) {
      constants_.m_p = Process::mp_;
      constants_.units = constants::GEVM2_TO_PB;
      constants_.pi = M_PI;
      constants_.alpha_em = constants::ALPHA_EM;
    }

    void FortranKTProcess::preparePhaseSpace() {
      mom_ip1_ = event_->oneWithRole(Particle::IncomingBeam1).momentum();
      mom_ip2_ = event_->oneWithRole(Particle::IncomingBeam2).momentum();

      defineVariable(
          y1_, Mapping::linear, kin_.cuts.central.rapidity_single(), {-6., 6.}, "First central particle rapidity");
      defineVariable(
          y2_, Mapping::linear, kin_.cuts.central.rapidity_single(), {-6., 6.}, "Second central particle rapidity");
      defineVariable(pt_diff_,
                     Mapping::linear,
                     kin_.cuts.central.pt_diff(),
                     {0., 50.},
                     "Transverse momentum difference between central particles");
      defineVariable(phi_pt_diff_,
                     Mapping::linear,
                     kin_.cuts.central.phi_diff(),
                     {0., 2. * M_PI},
                     "Central particles azimuthal angle difference");

      //===========================================================================================
      // feed phase space cuts to the common block
      //===========================================================================================

      kin_.cuts.central.pt_single().save(kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max);
      kin_.cuts.central.energy_single().save(kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max);
      kin_.cuts.central.eta_single().save(kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max);
      kin_.cuts.central.mass_sum().save(kincuts_.iinvm, kincuts_.invm_min, kincuts_.invm_max);
      kin_.cuts.central.pt_sum().save(kincuts_.iptsum, kincuts_.ptsum_min, kincuts_.ptsum_max);
      kin_.cuts.central.rapidity_diff().save(kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max);

      //===========================================================================================
      // feed run parameters to the common block
      //===========================================================================================

      genparams_.icontri = (int)kin_.incoming_beams.mode();
      if (kin_.incoming_beams.structureFunctions())
        genparams_.sfmod = kin_.incoming_beams.structureFunctions()->name();

      //-------------------------------------------------------------------------------------------
      // incoming beams information
      //-------------------------------------------------------------------------------------------

      //--- positive-z incoming beam
      genparams_.inp1 = kin_.incoming_beams.positive().momentum.pz();
      //--- check if first incoming beam is a heavy ion
      const HeavyIon in1 = (HeavyIon)kin_.incoming_beams.positive().pdg;
      if (in1) {
        genparams_.a_nuc1 = in1.A;
        genparams_.z_nuc1 = (unsigned short)in1.Z;
        if (genparams_.z_nuc1 > 1) {
          event_->oneWithRole(Particle::IncomingBeam1).setPdgId((pdgid_t)in1);
          event_->oneWithRole(Particle::OutgoingBeam1).setPdgId((pdgid_t)in1);
        }
      } else
        genparams_.a_nuc1 = genparams_.z_nuc1 = 1;

      //--- negative-z incoming beam
      genparams_.inp2 = kin_.incoming_beams.negative().momentum.pz();
      //--- check if second incoming beam is a heavy ion
      const HeavyIon in2 = (HeavyIon)kin_.incoming_beams.negative().pdg;
      if (in2) {
        genparams_.a_nuc2 = in2.A;
        genparams_.z_nuc2 = (unsigned short)in2.Z;
        if (genparams_.z_nuc2 > 1) {
          event_->oneWithRole(Particle::IncomingBeam2).setPdgId((pdgid_t)in2);
          event_->oneWithRole(Particle::OutgoingBeam2).setPdgId((pdgid_t)in2);
        }
      } else
        genparams_.a_nuc2 = genparams_.z_nuc2 = 1;

      //-------------------------------------------------------------------------------------------
      // intermediate partons information
      //-------------------------------------------------------------------------------------------

      genparams_.iflux1 = (int)kin_.incoming_beams.positive().kt_flux;
      genparams_.iflux2 = (int)kin_.incoming_beams.negative().kt_flux;
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
      ktkin_.m_x = sqrt(mX2_);
      ktkin_.m_y = sqrt(mY2_);

      //--- compute the event weight
      return func_();
    }

    void FortranKTProcess::fillCentralParticlesKinematics() {
      //===========================================================================================
      // outgoing beam remnants
      //===========================================================================================

      pX_ = Momentum(evtkin_.px);
      pY_ = Momentum(evtkin_.py);
      // express these momenta per nucleon
      pX_ *= 1. / genparams_.a_nuc1;
      pY_ *= 1. / genparams_.a_nuc2;

      //===========================================================================================
      // intermediate partons
      //===========================================================================================

      const Momentum mom_par1 = mom_ip1_ - pX_, mom_par2 = mom_ip2_ - pY_;
      event_->oneWithRole(Particle::Parton1).setMomentum(mom_par1);
      event_->oneWithRole(Particle::Parton2).setMomentum(mom_par2);
      event_->oneWithRole(Particle::Intermediate).setMomentum(mom_par1 + mom_par2);

      //===========================================================================================
      // central system
      //===========================================================================================

      auto& oc = (*event_)[Particle::CentralSystem];  // retrieve all references
                                                      // to central system particles
      for (int i = 0; i < evtkin_.nout; ++i) {
        auto& p = oc[i];  // retrieve a reference to the specific particle
        p.setPdgId((long)evtkin_.pdg[i]);
        p.setStatus(Particle::Status::FinalState);
        p.setMomentum(Momentum(evtkin_.pc[i]));
      }
    }
  }  // namespace proc
}  // namespace cepgen
