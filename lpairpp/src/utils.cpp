#include "utils.h"

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

double GetMassFromPDGId(int pdgId_)
{
  switch(abs(pdgId_)) {
  case 1:    return 0.33;           // d (from PYTHIA6.4)
  case 2:    return 0.33;           // u (from PYTHIA6.4)
  case 11:   return 0.510998928e-3; // electron
  case 13:   return 0.1056583715;   // muon
  case 15:   return 1.77682;        // tau
  case 21:   return 0.;             // gluon
  case 22:   return 0.;             // photon
  case 211:  return 0.13957018;     // pi+
  case 111:  return 0.1349766;      // pi0
  case 2101: return 0.57933;        // (ud)0 (from PYTHIA6.4)
  case 2103: return 0.77133;        // (ud)1 (from PYTHIA6.4)
  case 2203: return 0.77133;        // (uu)1 (from PYTHIA6.4)
  case 2212: return 0.938272046;    // proton
  default:   return -1.;
  }
}

bool PSF(double q2_, double mX2_, double* sigT_, double* w1_, double* w2_)
{
  int nBin;
  double xBin, dx, nu2, logqq0, gd2;
  double sigLow, sigHigh;
  double mX = std::sqrt(mX2_);
  double mP = GetMassFromPDGId(2212);
  double mPI = 0.135; //FIXME pi0 mass ???

  if (mX>=mP+mPI && mX<1.99) {
    if (mX<1.11) {
      nBin = 0;
      xBin = mX-mP-mPI;
      dx = 1.11-mP-mPI; // Delta w bin sizes
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

  sigLow = (nBin==0) ? 0. : exp(abrass[nBin]+bbrass[nBin]*logqq0+cbrass[nBin]*std::pow(fabs(logqq0), 3))*gd2;
  sigHigh = exp(abrass[nBin+1]+bbrass[nBin+1]*logqq0+cbrass[nBin+1]*std::pow(fabs(logqq0), 3))*gd2;

  *sigT_ = sigLow+xBin*(sigHigh-sigLow)/dx;
  *w1_ = (mX2_-std::pow(mP, 2))/(8.*std::pow(pi, 2)*mP*alphaF)*muBarn*(*sigT_);
  *w2_ = (*w1_)*q2_/(q2_-nu2);

  return true;
}

void Map(double expo_, double xmin_, double xmax_, double* out_, double* dout_)
{
  double y, out;
  y = xmax_/xmin_;
  out = xmin_*std::pow(y, expo_);
  *out_ = out;
  *dout_ = out*log(y);
#ifdef DEBUG
  std::cout << "=====================================" << std::endl;
  std::cout << "[Map] [DEBUG]"
            << "\n  min = " << xmin_
            << "\n  max = " << xmax_
            << "\n  max/min = " << y
            << "\n  exponent = " << expo_
            << "\n  output = " << *out_
            << "\n  d(output) = "<< *dout_
            << std::endl;
  std::cout << "=====================================" << std::endl;
#endif
}

void Mapla(double y_, double z_, int u_, double xm_, double xp_, double* x_, double* d_)
{
  double xmb, xpb, c, yy, zz, alp, alm, am, ap, ax;

  xmb = xm_-y_-z_;
  xpb = xp_-y_-z_;
  c = -4.*y_*z_;
  alp = std::sqrt(std::pow(xpb, 2)+c);
  alm = std::sqrt(std::pow(xmb, 2)+c);
  am = xmb+alm;
  ap = xpb+alp;
  yy = ap/am;
  zz = std::pow(yy, u_);

  *x_ = y_+z_+(am*zz-c/(am*zz))/2.;
  ax = std::sqrt(std::pow(*x_-y_-z_, 2)+c);
  *d_ = ax*log(yy);
}

void Lorenb(double u_, double ps_[4], double pi_[4], double pf_[4])
{
  double fn;

  if (ps_[3]!=u_) {
    pf_[3] = (pi_[3]*ps_[3]+pi_[2]*ps_[2]+pi_[1]*ps_[1]+pi_[0]*ps_[0])/u_;
    fn = (pf_[3]+pi_[3])/(ps_[3]+u_);
    pf_[0] = pi_[0]+fn*ps_[0];
    pf_[1] = pi_[1]+fn*ps_[1];
    pf_[2] = pi_[2]+fn*ps_[2];
  }
  else {
    std::copy(pi_, pi_+4, pf_);
  }
}
