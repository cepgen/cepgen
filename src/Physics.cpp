#ifndef __CINT__
#include "Physics.h"

PhysicsBoundaries::PhysicsBoundaries() :
  wmin(20.), wmax(0.),
  q2min(4.), q2max(100.),
  zmin(0.), zmax(1.)
{}

PhysicsBoundaries::~PhysicsBoundaries()
{}

/*double
GetBRFromProcessId(Particle::ParticleCode vmId_)
{
  switch ((VMDecay)vmId_) {
  case RHO_TO_PIPI:         return 1.0;    // rho0->pi+ pi-
  case OMEGA_TO_PIPI:       return 0.0221; // omega->pi+ pi-
  case PHI_TO_KK:           return 0.491;  // phi->K+ K-
  case PHI_TO_KLKS:         return 0.344;  // phi->KL0 KS0 //FIXME FIXME FIXME
  case JPSI_TO_LL:          return 0.0598; // J/psi->l+ l-
  case PSIP_TO_LLX:         return 0.0425; // psi'->l+ l- X
  case UPS1S_TO_LL:         return 0.0250; // Upsilon(1s)->l+ l-
  case UPS2S_TO_LLX:        return 0.0200; // Upsilon(2s)->l+ l- X
  case UPS3S_TO_LLX:        return 0.0217; // Upsilon(3s)->l+ l- X
  case RHO1450_TO_PIPIRHO0: // rho(1450)->pi+ pi- rho0
  case PHI1680_TO_KKBAR: // phi(1680)->K Kbar
///  case RHO_770_0:         return 1.0;    // rho0->pi+ pi-
  case OMEGA_782:       return 0.0221; // omega->pi+ pi-
  case PHI_TO_KK:           return 0.491;  // phi->K+ K-
  case PHI_TO_KLKS:         return 0.344;  // phi->KL0 KS0 //FIXME FIXME FIXME
  case JPSI_TO_LL:          return 0.0598; // J/psi->l+ l-
  case PSIP_TO_LLX:         return 0.0425; // psi'->l+ l- X
  case UPS1S_TO_LL:         return 0.0250; // Upsilon(1s)->l+ l-
  case UPS2S_TO_LLX:        return 0.0200; // Upsilon(2s)->l+ l- X
  case UPS3S_TO_LLX:        return 0.0217; // Upsilon(3s)->l+ l- X
  case RHO1450_TO_PIPIRHO0: // rho(1450)->pi+ pi- rho0
  case PHI1680_TO_KKBAR: // phi(1680)->K Kbar
  default: return -1;
  }
}*/

/*Particles
VMDecayer(Particle part_, GenericHadroniser *had_)
{
  if (part_.status!=1) {
    error.str(""); error << __PRETTY_FUNCTION__ << " ERROR: Particle has status" << part_.status;
    throw std::runtime_error(error.str());
  }
  if (part_.M()<2.5) {
    error.str(""); error << __PRETTY_FUNCTION__ << " ERROR: Particle has a too small mass (" << part_.M() << " GeV)";
    throw std::runtime_error(error.str());
  }
  //if (iret<=-2)...
  if (!had_->Hadronise(&part_)) {

  }
}*/

void
Lorenb(double u_, const Particle::Momentum& ps_, double pi_[4], double pf_[4])
{
  double fn;

  if (ps_.E()!=u_) {
    pf_[3] = (pi_[3]*ps_.E()+pi_[2]*ps_.Pz()+pi_[1]*ps_.Py()+pi_[0]*ps_.Px())/u_;
    fn = (pf_[3]+pi_[3])/(ps_.E()+u_);
    for (unsigned int i=0; i<3; i++) {
      pf_[i] = pi_[i]+fn*ps_.P(i);
    }
  }
  else {
    std::copy(pi_, pi_+4, pf_);
  }
}
#endif

// values of a, b, c provided from the fits on ep data and retrieved from
// http://dx.doi.org/10.1016/0550-3213(76)90231-5 with 1.110 <= w2 <=1.990

double abrass[56] = {5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178,6.125,5.999,
                     5.769,5.622,5.431,5.288,5.175,5.131,5.003,5.065,5.045,5.078,
                     5.145,5.156,5.234,5.298,5.371,5.457,5.543,5.519,5.465,5.384,
                     5.341,5.320,5.275,5.290,5.330,5.375,5.428,5.478,5.443,5.390,
                     5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158,5.159,5.178,
                     5.182,5.195,5.160,5.195,5.163,5.172};
double bbrass[56] = {0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,
                     2.041,2.089,2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,
                     2.609,2.678,2.771,2.890,2.982,3.157,3.183,3.315,3.375,3.450,
                     3.477,3.471,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519,
                     4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322,5.545,5.623,
                     5.775,5.894,6.138,6.151,6.301,6.542};
double cbrass[56] = { 0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080,-0.065,-0.056,
                     -0.065,-0.056,-0.043,-0.034,-0.054,-0.018,-0.046,-0.015,-0.029,-0.048,
                     -0.032,-0.045,-0.084,-0.115,-0.105,-0.159,-0.164,-0.181,-0.203,-0.223,
                     -0.245,-0.254,-0.239,-0.302,-0.299,-0.318,-0.383,-0.393,-0.466,-0.588,
                     -0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798,-1.048,-0.980,
                     -1.021,-1.092,-1.313,-1.341,-1.266,-1.473};

bool
PSF(double q2_, double mX2_, double* sigT_, double* w1_, double* w2_)
{
  int nBin;
  double xBin, dx, nu2, logqq0, gd2;
  double sigLow, sigHigh;
  double mX = std::sqrt(mX2_);
  //const double m_min = Particle::GetMassFromPDGId(Particle::Proton)+0.135;
  const double m_min = 1.07, mP = 0.938;

  if (mX>=m_min && mX<1.99) {
    if (mX<1.11) {
      nBin = 0;
      xBin = mX-m_min;
      dx = 1.11-m_min; // Delta w bin sizes
    }
    else if (mX<1.77) { // w in [1.11, 1.77[
      dx = 0.015; // Delta w bin sizes
      nBin = (mX-1.11)/dx+1;
      xBin = fmod(mX-1.11, dx);
    }
    else { // w in [1.77, 1.99[
      dx = 0.02; // Delta w bin sizes
      nBin = (mX-1.77)/dx+45;
      xBin = fmod(mX-1.77, dx);
    }
  }
  else {
    *sigT_ = 0.;
    *w1_ = 0.;
    *w2_ = 0.;
    return false;
  }
  nu2 = std::pow((mX2_-q2_-std::pow(mP, 2))/(2.*mP), 2);
  logqq0 = log((nu2-q2_)/std::pow((mX2_-std::pow(mP, 2))/(2.*mP), 2))/2.;
  gd2 = std::pow(1./(1-q2_/.71), 4); // dipole form factor of the proton

  sigLow = (nBin!=0) ?
    exp(abrass[nBin-1]+bbrass[nBin-1]*logqq0+cbrass[nBin-1]*std::pow(fabs(logqq0), 3))*gd2 :
    0.;
  sigHigh =
    exp(abrass[nBin]  +bbrass[nBin]  *logqq0+cbrass[nBin]  *std::pow(fabs(logqq0), 3))*gd2;

  *sigT_ = sigLow+xBin*(sigHigh-sigLow)/dx;
  *w1_ = (mX2_-std::pow(mP, 2))/(8.*std::pow(Constants::Pi, 2)*mP*Constants::AlphaEM)/Constants::GeV2toBarn*1.e6*(*sigT_);
  *w2_ = (*w1_)*q2_/(q2_-nu2);

  return true;
}

double
ElasticFlux(double x_, double kt2_)
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
InelasticFlux(double x_, double kt2_, double mx_)
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

FormFactors
TrivialFormFactors()
{
  FormFactors ff;
  ff.FE = 1.;
  ff.FM = 1.;
  return ff;
}

FormFactors
ElasticFormFactors(double q2, double mi2)
{
  const double GE = pow(1.+q2/0.71, -2.), GM = 2.79*GE;
  FormFactors ff;
  ff.FE = (4.*mi2*GE*GE+q2*GM*GM)/(4.*mi2+q2);
  ff.FM = GM*GM;
  return ff;
}

FormFactors
SuriYennieFormFactors(double q2, double mi2, double mf2)
{
  // values extracted from experimental fits
  const double cc1 = 0.86926, // 0.6303
               cc2 = 2.23422, // 2.2049
               dd1 = 0.12549, // 0.0468
               cp = 0.96, // 1.23
               bp = 0.63, // 0.61
               rho = 0.585; // 1.05
  const double x = q2/(q2+mf2),
               dm2 = mf2-mi2,
               en = dm2+q2,
               tau = -q2/4./mi2,
               rhot = rho+q2;
  FormFactors ff;
  ff.FM = -(-cc1*std::pow(rho/rhot, 2)*dm2-cc2*mi2*std::pow(1.-x, 4)/(x*(x*cp-2*bp)+1.))/q2;
  ff.FE = (-tau*ff.FM+dd1*dm2*q2*(rho/rhot)*std::pow(dm2/en, 2)/(rhot*mi2))/(1.+en*en/(4.*mi2*q2));
  return ff;
}

FormFactors
FioreBrasseFormFactors(double q2, double mi2, double mf2)
{
  const double k = 2.*sqrt(mi2);
  double dummy, psfw1, psfw2; PSF(-q2, mf2, &dummy, &psfw1, &psfw2);
  FormFactors ff;
  ff.FM = psfw1*k/q2;
  ff.FE = psfw2/k;
  return ff;
}
