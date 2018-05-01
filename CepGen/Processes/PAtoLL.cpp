#include "CepGen/Processes/PAtoLL.h"
#include "CepGen/Physics/PDG.h"

using namespace CepGen::Process;

extern "C"
{
  extern void pa_ll_( double& );
  extern struct {
    int mode, pdg_l, a_nuc, z_nuc;
    double m_l;
    double inp1, inp2;
  } params_;
  extern struct {
    double m_p, units, pi, alpha_em;
  } constants_;
  extern struct {
   double q1t, q2t, phiq1t, phiq2t, y1, y2, ptdiff, phiptdiff, m_x, m_y;
  } ktkin_;
}

PAtoLL::PAtoLL() : GenericKTProcess( "patoll", "pA ↝ ɣɣ → l⁺l¯", { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } )
{
  constants_.m_p = ParticleProperties::mass( PDG::Proton );
  constants_.units = Constants::GeV2toBarn;
  constants_.pi = M_PI;
  constants_.alpha_em = Constants::alphaEM;
  params_.mode = 1;
  params_.pdg_l = (int)PDG::Muon;
  params_.m_l = ParticleProperties::mass( (PDG)params_.pdg_l );
  params_.a_nuc = 82;
  params_.z_nuc = 208;
  params_.inp1 = params_.inp2 = 6500.;
  ktkin_.q1t = ktkin_.q2t = ktkin_.phiq1t = ktkin_.phiq2t = ktkin_.y1 = ktkin_.y2 = ktkin_.ptdiff = ktkin_.phiptdiff = 0.;
  ktkin_.m_x = ktkin_.m_y = ParticleProperties::mass( PDG::Proton );
}

void
PAtoLL::preparePhaseSpace()
{
}

double
PAtoLL::computeKTFactorisedMatrixElement()
{
  double weight = 0.;
  pa_ll_( weight );
  return weight;
}

void
PAtoLL::fillCentralParticlesKinematics()
{
}
