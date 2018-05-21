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
}

void
PAtoLL::preparePhaseSpace()
{
  registerVariable( y1_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "First outgoing lepton rapidity" );
  registerVariable( y2_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "Second outgoing lepton rapidity" );
  registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::pt_diff], { 0., 50. }, "Leptons transverse momentum difference" );
  registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::phi_pt_diff], { 0., 2.*M_PI }, "Leptons azimuthal angle difference" );
}

double
PAtoLL::computeKTFactorisedMatrixElement()
{
  ktkin_.q1t = qt1_;
  ktkin_.q2t = qt2_;
  ktkin_.phiq1t = phi_qt1_;
  ktkin_.phiq2t = phi_qt2_;
  ktkin_.y1 = y1_;
  ktkin_.y2 = y2_;
  ktkin_.ptdiff = pt_diff_;
  ktkin_.phiptdiff = phi_pt_diff_;
  ktkin_.m_x = MX_;
  ktkin_.m_y = MY_;
  double weight = 0.;
  pa_ll_( weight );
  return weight;
}

void
PAtoLL::fillCentralParticlesKinematics()
{
}
