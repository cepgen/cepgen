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
  double m;
  switch(abs(pdgId_)) {
  case 2212: // proton
    //m = .938272046;
    m = .93830001354217529; // FROM LPAIR
    break;
  case 11: // electron
    m = .510998928e-3;
    break;
  case 13: // muon
    //m = .1056583715;
    m = .10570000112056732; // FROM LPAIR
    break;
  case 15: // tau
    m = 1.77682;
    break;
  case 22: // photon
    m = 0.;
    break;
  default:
    m = 0.;
    break;
  }
  return m;
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
  //std::cout << "===> map : a = " << y << ", b = " << expo_ << ", a^b = " << std::pow(y, expo_) << std::endl;
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

InputParameters::InputParameters() :
  p1mod(2), p2mod(2),
  pair(13),
  mcut(0),
  minpt(0.5), maxpt(-1.),
  minenergy(1.), maxenergy(-1.),
  mintheta(5.), maxtheta(175.),
  minmx(1.07), maxmx(320.),
  itmx(10),
  generation(true), store(false), debug(false)
{
  /*for (int i=0; i<MAX_HISTOS; i++) {
    this->plot[i] = new Gnuplot();
  }*/
}

InputParameters::~InputParameters()
{
  /*for (int i=0; i<MAX_HISTOS; i++) {
    delete[] this->plot[i];
  }//FIXME ???*/
  //delete[] this->plot;
}

void InputParameters::Dump()
{
  std::cout << "[InputParameters] Input Parameters DUMP =============" << std::endl;
  std::cout << "===== Kinematics =====" << std::endl;
  std::cout << "  Incoming protons-like particles" << std::endl;
  std::cout << "  -> Mode : " << p1mod << " and " << p2mod << std::endl;
  std::cout << "  -> Momenta : " << in1p << " and " << in2p << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "  Leptons pair : " << pair << std::endl;
  std::cout << "  Cuts mode : " << mcut << std::endl;
  std::cout << "  -> pT in range [" << minpt << ", " << maxpt << "]" << std::endl;
  std::cout << "  -> Energy in range [" << minenergy << ", " << maxenergy << "]" << std::endl;
  std::cout << "  -> Azimuthal angle theta in range [" << mintheta << ", " << maxtheta << "]" << std::endl;
  std::cout << "  Outgoing remnants' mass in range [" << minmx << ", " << maxmx << "]" << std::endl;
  std::cout << "===== VEGAS =====" << std::endl;
  std::cout << "  Maximum number of iterations : " << itmx << std::endl;
  std::cout << "===== General informaton =====" << std::endl;
  std::cout << "  Events generation ? " << generation << std::endl;
  std::cout << "  Events storage ? " << store << std::endl;
  std::cout << "  Debugging mode ? " << debug << std::endl;
  std::cout << "[InputParameters] End of DUMP =======================" << std::endl;
}

bool InputParameters::ReadConfigFile(std::string inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_.c_str(), std::fstream::in);
  if (!f.is_open()) {
    return false;
  }
#ifdef DEBUG
  std::cout << "[InputParameters::ReadConfigFile] [DEBUG] File " << inFile_ << " succesfully opened !" << std::endl;
  std::cout << "======================================================" << std::endl;
  std::cout << "Configuration file content : " << std::endl;
  std::cout << "======================================================" << std::endl;
#endif
  while (f >> key >> value) {
    if (key=="IEND") {
      int iend = (int)atoi(value.c_str());
      if (iend>1) {
        this->generation = true;
      }
    }
    else if (key=="NCVG") {
      this->ncvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Number of Vegas calls                               : " << this->ncvg << std::endl;
#endif
    }
    else if (key=="ITVG") {
      this->itvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Number of Vegas iterations                          : " << this->itvg << std::endl;
#endif
    }
    else if (key=="INPP") {
      this->in1p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << " * First incoming particles' momentum                  : " << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="PMOD") {
      this->p1mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * First incoming particles' mode                      : " << this->p1mod << " --> ";
      switch (this->p1mod) {
        case 1:
          std::cout << "electron";
          break;
        case 2:
        default:
          std::cout << "elastic proton [EPA]";
          break;
        case 11:
          std::cout << "dissociating proton [structure functions]";
          break;
        case 12:
          std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q**2 < 5 GeV**2]";
          break;
        case 101:
          std::cout << "dissociating proton [parton model, only valence quarks]";
          break;
        case 102:
          std::cout << "dissociating proton [parton model, only sea quarks]";
          break;
        case 103:
          std::cout << "dissociating proton [parton model, valence and sea quarks]";
          break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="INPE") {
      this->in2p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << " * Second incoming particles' momentum                 : " << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="EMOD") {
      this->p2mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Second incoming particles' mode                     : " << this->p2mod << " --> ";
      switch (this->p2mod) {
        case 1:
          std::cout << "electron";
          break;
        case 2:
        default:
          std::cout << "elastic proton [EPA]";
          break;
        case 11:
          std::cout << "dissociating proton [structure functions]";
          break;
        case 12:
          std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q**2 < 5 GeV**2]";
          break;
        case 101:
          std::cout << "dissociating proton [parton model, only valence quarks]";
          break;
        case 102:
          std::cout << "dissociating proton [parton model, only sea quarks]";
          break;
        case 103:
          std::cout << "dissociating proton [parton model, valence and sea quarks]";
          break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PAIR") {
      this->pair = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Outgoing leptons pair PDG id                        : " << this->pair << " --> ";
      switch (this->pair) {
        case 11:
        default:
          std::cout << "electrons";
          break;
        case 13:
          std::cout << "muons";
          break;
        case 15:
          std::cout << "taus";
          break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="MCUT") {
      this->mcut = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Set of cuts to apply on the total generation        : " << this->mcut << " --> ";
      switch (this->mcut) {
        case 3:
          std::cout << "cuts on at least one outgoing lepton";
          break;
        case 2:
          std::cout << "cuts on both the outgoing leptons";
          break;
        case 1:
          std::cout << "Vermaseren's hypothetical detector cuts";
          break;
        case 0:
        default:
          std::cout << "no cuts";
          break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PTCT") {
      this->minpt = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << " * Outgoing lepton pairs' minimal transverse momentum  : " << this->minpt << " GeV/c" << std::endl;
#endif
    }
    else if (key=="ECUT") {
      this->minpt = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << " * Outgoing lepton pairs' minimal energy               : " << this->minpt << " GeV/c" << std::endl;
#endif
    }
    else if (key=="ITMX") {
      this->itmx = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << " * Maximal number of VEGAS iterations                  : " << this->itmx << std::endl;
#endif
    }
    else {
      std::cout << "[InputParameters::ReadConfigFile] <WARNING> Unrecognized argument : [" << key << "] = " << value << std::endl;
    }
  }
  f.close();
  std::cout << "======================================================" << std::endl;
  return true;
}

bool InputParameters::StoreConfigFile(std::string outFile_)
{
  std::ofstream f;
  f.open(outFile_.c_str(), std::fstream::out | std::fstream::trunc);
  if (!f.is_open()) {
    return false;
  }
  // ...
  f.close();
  return true;
}

Cuts::Cuts() :
  ptmin(3.), ptmax(-1.), emin(0.), emax(-1.),
  thetamin(0.), thetamax(180.)
{}

Cuts::~Cuts() {}
