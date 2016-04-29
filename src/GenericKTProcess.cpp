#include "GenericKTProcess.h"

GenericKTProcess::GenericKTProcess(std::string name_="<generic process>",
                                   Particle::ParticleCode p1_=Particle::Photon,
                                   Particle::ParticleCode p2_=Particle::Photon) :
  GenericProcess(name_+" (kT-factorization approach)"),
  kIntermediatePart1(p1_), kIntermediatePart2(p2_)
{}

GenericKTProcess::~GenericKTProcess()
{}

void
GenericKTProcess::AddEventContent()
{
  IncomingState is; OutgoingState os;
  is.insert(ParticleWithRole(Particle::IncomingBeam1,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::IncomingBeam2,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::Parton1,          kIntermediatePart1));
  is.insert(ParticleWithRole(Particle::Parton2,          kIntermediatePart2));
  os.insert(ParticleWithRole(Particle::OutgoingBeam1,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::OutgoingBeam2,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::CentralParticle1, Particle::Muon));
  os.insert(ParticleWithRole(Particle::CentralParticle2, Particle::Muon));
  GenericProcess::SetEventContent(is, os);
}

int
GenericKTProcess::GetNdim(ProcessMode process_mode_) const
{
  switch (process_mode_) {
    default:
    case ElasticElastic:     return 8;
    case ElasticInelastic:
    case InelasticElastic:   return 9;
    case InelasticInelastic: return 10;
  }
}

double
GenericKTProcess::ComputeWeight()
{
  PrepareKTKinematics();
  ComputeOutgoingPrimaryParticlesMasses();
  
  const double jac = ComputeJacobian(),
               integrand = ComputeKTFactorisedMatrixElement(),
               weight = jac*integrand;
  DebugInsideLoop(Form("Jacobian = %f\n\tIntegrand = %f\n\tdW = %f", jac, integrand, weight));
  
  return weight;
}

void
GenericKTProcess::ComputeOutgoingPrimaryParticlesMasses()
{
  switch (fCuts.kinematics) {
    case 0: default: { Error("This kT factorisation process is intended for p-on-p collisions! Aborting!"); exit(0); break; }
    case 1: 
      fMX = GetParticle(Particle::IncomingBeam1)->M();
      fMY = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 2:
      fMX = GetParticle(Particle::IncomingBeam1)->M();
      fMY = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      break;
    case 3:
      fMX = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      fMY = GetParticle(Particle::IncomingBeam2)->M();
      break;
    case 4:
      fMX = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(8);
      fMY = fCuts.mxmin+(fCuts.mxmax-fCuts.mxmin)*x(9);
      break;
  }
  DebugInsideLoop(Form("outgoing remnants invariant mass: %f / %f (%.2f < M(X/Y) < %.2f)", fMX, fMY, fCuts.mxmin, fCuts.mxmax));
}

void
GenericKTProcess::FillPrimaryParticlesKinematics()
{
  //=================================================================
  //     outgoing protons
  //=================================================================
  Particle *op1 = GetParticle(Particle::OutgoingBeam1),
           *op2 = GetParticle(Particle::OutgoingBeam2);
  switch (fCuts.kinematics) {
    case 1:
      op1->status = Particle::FinalState;
      op2->status = Particle::FinalState;
      break;
    case 2:
      op1->status = Particle::Undecayed; op1->SetM();
      op2->status = Particle::FinalState;
      break;
    case 3:
      op1->status = Particle::FinalState;
      op2->status = Particle::Undecayed; op2->SetM();
      break;
    case 4:
      op1->status = Particle::Undecayed; op1->SetM();
      op2->status = Particle::Undecayed; op2->SetM();
      break;    
  }
  
  if (!op1->SetMomentum(fPX)) { Error(Form("Invalid outgoing proton 1: energy: %.2f", fPX.E())); }
  if (!op2->SetMomentum(fPY)) { Error(Form("Invalid outgoing proton 2: energy: %.2f", fPY.E())); }
}

double
GenericKTProcess::ElasticFlux(double x_, double kt2_) const
{
  double f_ela;

  const double mp2 = pow(Particle::GetMassFromPDGId(Particle::Proton), 2);

  const double Q2_ela = (kt2_+pow(x_, 2)*mp2)/(1.-x_);
  const double G_dip = 1./pow(1.+Q2_ela/0.71, 2);
  const double G_E = G_dip;
  const double G_M = 2.79*G_dip;

  const double ela1 = pow(kt2_/(kt2_+pow(x_, 2)*mp2), 2);
  const double ela2 = (4.*mp2*pow(G_E, 2)+Q2_ela*pow(G_M, 2))/(4.*mp2+Q2_ela);
  //const double ela3 = 1.-(Q2_ela-kt2_)/Q2_ela;
  //const double ela3 = 1.-pow(x_, 2)*mp2/Q2_ela/(1.-x_);
  //f_ela = alpha_em/Constants::Pi*(1.-x_+pow(x_, 2)/4.)*ela1*ela2*ela3/kt2_;
  f_ela = Constants::AlphaEM/Constants::Pi*ela1*ela2/Q2_ela;
  //f_ela = Constants::AlphaEM/Constants::Pi*((1.-x_)*ela1*ela2*ela3+pow(x_, 2)/2.*pow(G_M, 2))/kt2_;

  return f_ela;
}

double
GenericKTProcess::InelasticFlux(double x_, double kt2_, double mx_) const
{
  double f_ine;

  const double mx2 = mx_*mx_,
               mp2 = pow(Particle::GetMassFromPDGId(Particle::Proton), 2);
  //const double mpi = pow(Particle::GetMassFromPDGId(Particle::PiZero), 2);

  const double Q02 = 0.8; // introduced to shift the Q2 scale
  double term1, term2;
  double f_aux;

  // F2 structure function
  const double Q2min = 1./(1.-x_)*(x_*(mx2-mp2)+pow(x_, 2)*mp2), Q2 = kt2_/(1.-x_)+Q2min;
  float x_Bjorken = Q2/(Q2+mx2-mp2);

  float mu2 = Q2+Q02; // scale is shifted

  float xuv, xdv, xus, xds, xss, xg;
  grv95lo_(x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg);
  DebugInsideLoop(Form("Form factor content at xB = %e (scale = %f GeV^2):\n\t"
                       "  valence quarks: u / d     = %e / %e\n\t"
                       "  sea quarks:     u / d / s = %e / %e / %e\n\t"
                       "  gluons:                   = %e",
                       x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg));

  const double F2_aux = 4./9.*(xuv + 2.*xus)
                      + 1./9.*(xdv + 2.*xds)
                      + 1./9.*2.*xss;

  /*F2_aux = 4./9.*(xuv + 2.*xus)
         + 1./9.*(0. + 2.*xds)
         + 1./9.*2.*xss;*/

  // F2 corrected for low Q^2 behaviour
  const double F2_corr = Q2/(Q2+Q02)*F2_aux;

  ///////term1 = pow(1.- x_/2.*(mx2-mp2+Q2)/Q2, 2);
  //term1 = (1.-x_*(mx2-mp2+Q2)/Q2);
  term1 = (1.-(Q2-kt2_)/Q2);
  //term1 = (1.-Q2min/Q2);
  //term1 = 1.;
  term2 = pow(kt2_/(kt2_+x_*(mx2-mp2)+pow(x_, 2)*mp2), 2);

  f_aux = F2_corr/(mx2+Q2-mp2)*term1*term2;

  f_ine = Constants::AlphaEM/Constants::Pi*(1.-x_)*f_aux/kt2_;

  return f_ine;
}
