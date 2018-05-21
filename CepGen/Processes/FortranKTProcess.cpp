#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/Physics/PDG.h"

using namespace CepGen::Process;

extern "C"
{
  struct Constants {
    double m_p, units, pi, alpha_em;
  };
  struct Parameters {
    int icontri, iflux1, iflux2, sfmod, pdg_l, a_nuc1, z_nuc1, a_nuc2, z_nuc2, idum;
    double inp1, inp2;
  };
  struct KtKinematics {
   double q1t, q2t, phiq1t, phiq2t, y1, y2, ptdiff, phiptdiff, m_x, m_y;
  };
  struct KinematicsCuts {
    int ipt, iene, ieta, idely;
    double pt_min, pt_max, ene_min, ene_max, eta_min, eta_max, dely_min, dely_max;
  };
  struct EventKinematics {
    double p10, p1x, p1y, p1z, p20, p2x, p2y, p2z;
    double px0, pxx, pxy, pxz, py0, pyx, pyy, pyz;
  };

  extern Constants constants_;
  extern Parameters params_;
  extern KtKinematics ktkin_;
  extern KinematicsCuts kincuts_;
  extern EventKinematics evtkin_;
}

FortranKTProcess::FortranKTProcess( const char* name, const char* descr, std::function<void( double& )> func ) :
  GenericKTProcess( name, descr, { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } ),
  func_( func )
{
  constants_.m_p = ParticleProperties::mass( PDG::Proton );
  constants_.units = Constants::GeV2toBarn;
  constants_.pi = M_PI;
  constants_.alpha_em = Constants::alphaEM;
}

void
FortranKTProcess::preparePhaseSpace()
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
  params_.sfmod = (int)cuts_.structure_functions;
  params_.pdg_l = (int)cuts_.central_system[0];
  params_.iflux1 = cuts_.incoming_beams.first.kt_flux;
  params_.iflux2 = cuts_.incoming_beams.second.kt_flux;
  params_.a_nuc1 = cuts_.incoming_beams.first.hi.A;
  params_.z_nuc1 = cuts_.incoming_beams.first.hi.Z;
  params_.a_nuc2 = cuts_.incoming_beams.second.hi.A;
  params_.z_nuc2 = cuts_.incoming_beams.second.hi.Z;
  params_.inp1 = cuts_.incoming_beams.first.pz;
  params_.inp2 = cuts_.incoming_beams.second.pz;
}

double
FortranKTProcess::computeKTFactorisedMatrixElement()
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
  func_( weight );
  return weight;
}

void
FortranKTProcess::fillCentralParticlesKinematics()
{
  const Particle::Momentum& mom_ip1 = event_->getOneByRole( Particle::IncomingBeam1 ).momentum();
  const Particle::Momentum& mom_ip2 = event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

  //===============================================================================================
  // outgoing beam remnants
  //===============================================================================================

  PX_ = Particle::Momentum( evtkin_.pxx, evtkin_.pxy, evtkin_.pxz, evtkin_.px0 );
  PX_ *= 1./params_.z_nuc1;
  PY_ = Particle::Momentum( evtkin_.pyx, evtkin_.pyy, evtkin_.pyz, evtkin_.py0 );
  PY_ *= 1./params_.z_nuc2;

  //===============================================================================================
  // intermediate partons
  //===============================================================================================

  const Particle::Momentum mom_par1 = mom_ip1-PX_, mom_par2 = mom_ip2-PY_;
  Particle& par1 = event_->getOneByRole( Particle::Parton1 );
  par1.setMomentum( mom_par1 );
  Particle& par2 = event_->getOneByRole( Particle::Parton2 );
  par2.setMomentum( mom_par2 );
  event_->getOneByRole( Particle::Intermediate ).setMomentum( mom_par1+mom_par2 );

  //===============================================================================================
  // central system
  //===============================================================================================

  Particle& ol1 = event_->getByRole( Particle::CentralSystem )[0];
  ol1.setPdgId( (PDG)params_.pdg_l, +1 );
  ol1.setStatus( Particle::FinalState );
  ol1.setMomentum( Particle::Momentum( evtkin_.p1x, evtkin_.p1y, evtkin_.p1z, evtkin_.p10 ) );

  Particle& ol2 = event_->getByRole( Particle::CentralSystem )[1];
  ol2.setPdgId( (PDG)params_.pdg_l, -1 );
  ol2.setStatus( Particle::FinalState );
  ol2.setMomentum( Particle::Momentum( evtkin_.p2x, evtkin_.p2y, evtkin_.p2z, evtkin_.p20 ) );
}

