#include "CepGen/Processes/PAtoLL.h"
#include "CepGen/Physics/PDG.h"

using namespace CepGen::Process;

extern "C"
{
  extern void pa_ll_( double& );
  extern struct {
    double m_p, units, pi, alpha_em;
  } constants_;
  extern struct {
    int icontri, imode, sfmod, pdg_l, a_nuc, z_nuc;
    double m_l;
    double inp1, inp2;
  } params_;
  extern struct {
   double q1t, q2t, phiq1t, phiq2t, y1, y2, ptdiff, phiptdiff, m_x, m_y;
  } ktkin_;
  extern struct {
    int ipt, iene, ieta, idely;
    double pt_min, pt_max, ene_min, ene_max, eta_min, eta_max, dely_min, dely_max;
  } kincuts_;
  extern struct {
    double p10, p1x, p1y, p1z, p20, p2x, p2y, p2z;
    double px0, pxx, pxy, pxz, py0, pyx, pyy, pyz;
  } evtkin_;
}

PAtoLL::PAtoLL() : GenericKTProcess( "patoll", "pA ↝ ɣɣ → l⁺l¯", { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } )
{
  constants_.m_p = ParticleProperties::mass( PDG::Proton );
  constants_.units = Constants::GeV2toBarn;
  constants_.pi = M_PI;
  constants_.alpha_em = Constants::alphaEM;
}

void
PAtoLL::preparePhaseSpace()
{
  registerVariable( y1_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "First outgoing lepton rapidity" );
  registerVariable( y2_, Mapping::linear, cuts_.cuts.central[Cuts::rapidity_single], { -6., 6. }, "Second outgoing lepton rapidity" );
  registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::pt_diff], { 0., 50. }, "Leptons transverse momentum difference" );
  registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central[Cuts::phi_pt_diff], { 0., 2.*M_PI }, "Leptons azimuthal angle difference" );

  // feed phase space cuts to the common block
  cuts_.cuts.central[Cuts::pt_single].save( (bool&)kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max );
  cuts_.cuts.central[Cuts::energy_single].save( (bool&)kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max );
  cuts_.cuts.central[Cuts::eta_single].save( (bool&)kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max );
  cuts_.cuts.central[Cuts::rapidity_diff].save( (bool&)kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max );

  // feed run parameters to the common block
  params_.icontri = (int)cuts_.mode;
  params_.imode = 2; // Budnev flux
  params_.sfmod = (int)cuts_.structure_functions;
  params_.pdg_l = (int)PDG::Muon;
  params_.m_l = ParticleProperties::mass( (PDG)params_.pdg_l );
  params_.a_nuc = 208;
  params_.z_nuc = 82;
  params_.inp1 = cuts_.inp.first;
  params_.inp2 = cuts_.inp.second;
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
  Particle& ol1 = event_->getByRole( Particle::CentralSystem )[0];
  ol1.setPdgId( (PDG)params_.pdg_l, +1 );
  ol1.setStatus( Particle::FinalState );
  ol1.setMomentum( Particle::Momentum( evtkin_.p1x, evtkin_.p1y, evtkin_.p1z, evtkin_.p10 ) );

  Particle& ol2 = event_->getByRole( Particle::CentralSystem )[1];
  ol2.setPdgId( (PDG)params_.pdg_l, -1 );
  ol2.setStatus( Particle::FinalState );
  ol2.setMomentum( Particle::Momentum( evtkin_.p2x, evtkin_.p2y, evtkin_.p2z, evtkin_.p20 ) );

  PX_ = Particle::Momentum( evtkin_.pxx, evtkin_.pxy, evtkin_.pxz, evtkin_.px0 );
  PY_ = Particle::Momentum( evtkin_.pyx, evtkin_.pyy, evtkin_.pyz, evtkin_.py0 );
}
