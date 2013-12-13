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
    m = .938272046;
    //m = .93830001354217529; // FROM LPAIR
    break;
  case 11: // electron
    m = .510998928e-3;
    break;
  case 13: // muon
    m = .1056583715;
    //m = .10570000112056732; // FROM LPAIR
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
  //mintheta(0.), maxtheta(180.),
  minq2(0.), maxq2(1.e5),
  minmx(1.07), maxmx(320.),
  ncvg(14000), itvg(10),
  //ncvg(100000), itvg(10),
  ntreat(1), npoints(100),
  generation(true), store(false), debug(false),
  maxgen(1e5),
  gpdf(5), spdf(4), qpdf(12),
  symmetrise(true)
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

void InputParameters::SetEtaRange(double etamin_, double etamax_)
{
  this->mintheta = 2.*atan(exp(-etamax_))/pi*180.;
  this->maxtheta = 2.*atan(exp(-etamin_))/pi*180.;
#ifdef DEBUG
  std::cout << "[InputParameters::SetEtaRange] [DEBUG]"
	    << "\n\teta(min) = " << std::setw(5) << etamin_ << " -> theta(min) = " << this->mintheta
	    << "\n\teta(max) = " << std::setw(5) << etamax_ << " -> theta(max) = " << this->maxtheta
	    << std::endl;
#endif
}

void InputParameters::Dump()
{
  std::string cutsmode, particles;

  switch(mcut) {
    case 1:
      cutsmode = "Vermaseren"; break;
    case 2:
      cutsmode = "both"; break;
    case 3:
      cutsmode = "single"; break;
    case 0:
    default:
      cutsmode = "none"; break;
  }
  switch(pair) {
    case 11:
      particles = "electrons"; break;
    case 13:
    default: 
      particles = "muons"; break;
    case 15:
      particles = "taus"; break;
  }
  std::cout 
    << std::left
    << "[InputParameters::Dump] BEGINNING dump ===============" << std::endl << std::endl
    << " _" << std::setfill('_') << std::setw(52) << "_/¯ INCOMING- AND OUTGOING KINEMATICS ¯\\_" << std::setfill(' ') << "_ " << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Incoming protons-like particles " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Mode" << std::setw(4) << p1mod << ", " << std::setw(4) << p2mod << " |" << std::endl
    << "| " << std::setw(40) << "Momenta [GeV/c]" << std::setw(4) << in1p << ", " << std::setw(4) << in2p << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Outgoing leptons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Pair" << std::setw(2) << pair << " " << std::setw(7) << particles << " |" << std::endl
    << "| " << std::setw(40) << "Cuts mode" << std::setw(1) << mcut << " (" << std::setw(6) << cutsmode << ")" << " |" << std::endl
    << "| " << std::setw(40) << "pT [GeV/c]" << "[" << std::setw(3) << minpt << ", " << std::setw(3) << maxpt << "]" << " |" << std::endl
    << "| " << std::setw(40) << "Energy [GeV]" << "[" << std::setw(3) << minenergy << ", " << std::setw(3) << maxenergy << "]" << " |" << std::endl
    << "| " << std::setw(40) << "Polar angle theta [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw(3) << maxtheta << "]" << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Outgoing remnants " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Minimal mass [GeV/c**2]" << std::setw(10) << minmx << " |" << std::endl
    << "| " << std::setw(40) << "Maximal mass [GeV/c**2]" << std::setw(10) << maxmx << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(52) << "_/¯ VEGAS INTEGRATION PARAMETERS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Maximum number of iterations" << std::setw(10) << itvg << " |" << std::endl
    << "| " << std::setw(40) << "Number of function calls" << std::setw(10) << ncvg << " |" << std::endl
    << "| " << std::setw(40) << "Number of points to try per bin" << std::setw(10) << npoints << " |" << std::endl
    << "| " << std::setw(40) << "Is the integration smoothed (TREAT) ? " << std::setw(10) << ntreat << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(52) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Events generation ? " << std::setw(10) << generation << " |" << std::endl
    << "| " << std::setw(40) << "Number of events to generate" << std::setw(10) << maxgen << " |" << std::endl
    << "| " << std::setw(40) << "Events storage ? " << std::setw(10) << store << " |" << std::endl
    << "| " << std::setw(40) << "Debugging mode ? " << std::setw(10) << debug << " |" << std::endl
    << "| " << std::setw(40) << "Is Output file opened ? " << std::setw(10) << file->is_open() << " |" << std::endl
    //<< "| " << std::setw(40) << "Is Debug file opened ? " << std::setw(10) << file_debug->is_open() << " |" << std::endl
    << "|_" << std::right << std::setfill('_') << std::setw(52) << "_|" << std::left << std::endl
    //<< " -" << std::right << std::setfill('-') << std::setw(52) << "- " << std::left << std::endl
    << std::endl
    << "[InputParameters::Dump] END of dump ==================" << std::endl;
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
  std::cout << "[InputParameters::ReadConfigFile] [DEBUG] File " << inFile_ << " succesfully opened !" << std::endl
            << "======================================================" << std::endl
            << "Configuration file content : " << std::endl
            << "======================================================" << std::endl
            << std::left;
#endif
  while (f >> key >> value) {
    //std::cout << std::setw(60) << "[" << key << "] = " << value << std::endl;
    //if (strncmp(key.c_str(), "#")==0) continue; // FIXME need to ensure there is no extra space before !
    if (key=="IEND") {
      int iend = (int)atoi(value.c_str());
      if (iend>1) {
        this->generation = true;
      }
    }
    else if (key=="NCVG") {
      this->ncvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of function calls" << this->ncvg << std::endl;
#endif
    }
    else if (key=="NCSG") {
      this->npoints = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of points to probe" << this->npoints << std::endl;
#endif
    }
    else if (key=="ITVG") {
      this->itvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of Vegas iterations" << this->itvg << std::endl;
#endif
    }
    else if (key=="INPP") {
      this->in1p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * First incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="PMOD") {
      this->p1mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * First incoming particles' mode" << this->p1mod << " --> ";
      switch (this->p1mod) {
        case 1:          std::cout << "electron"; break;
        case 2: default: std::cout << "elastic proton"; break;
        case 11:         std::cout << "dissociating proton [structure functions]"; break;
        case 12:         std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
        case 101:        std::cout << "dissociating proton [parton model, only valence quarks]"; break;
        case 102:        std::cout << "dissociating proton [parton model, only sea quarks]"; break;
        case 103:        std::cout << "dissociating proton [parton model, valence and sea quarks]"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="INPE") {
      this->in2p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Second incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="EMOD") {
      this->p2mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Second incoming particles' mode" << this->p2mod << " --> ";
      switch (this->p2mod) {
        case 1:           std::cout << "electron"; break;
        case 2: default:  std::cout << "elastic proton [EPA]"; break;
        case 11:          std::cout << "dissociating proton [structure functions]"; break;
        case 12:          std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
        case 101:         std::cout << "dissociating proton [parton model, only valence quarks]"; break;
        case 102:         std::cout << "dissociating proton [parton model, only sea quarks]"; break;
        case 103:         std::cout << "dissociating proton [parton model, valence and sea quarks]"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PAIR") {
      this->pair = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Outgoing leptons' PDG id   " << this->pair << " --> ";
      switch (this->pair) {
        case 11: default: std::cout << "electrons"; break;
        case 13:          std::cout << "muons"; break;
        case 15:          std::cout << "taus"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="MCUT") {
      this->mcut = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Set of cuts to apply on the total generation  " << this->mcut << " --> ";
      switch (this->mcut) {
        case 0: default: std::cout << "no cuts"; break;
        case 3:          std::cout << "cuts on at least one outgoing lepton"; break;
        case 2:          std::cout << "cuts on both the outgoing leptons"; break;
        case 1:          std::cout << "Vermaseren's hypothetical detector cuts"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PTCT") {
      this->minpt = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Single outgoing lepton's minimal transverse momentum" << this->minpt << " GeV/c" << std::endl;
#endif
    }
    else if (key=="ECUT") {
      this->minenergy = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Single outgoing lepton's minimal energy" << this->minenergy << " GeV" << std::endl;
#endif
    }
    else if (key=="NTRT") {
      this->ntreat = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of TREAT calls" << this->ntreat << std::endl;
#endif
    }
    else if (key=="NGEN") {
      this->maxgen = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of events to generate" << this->maxgen << std::endl;
#endif
    }
    else if (key=="THMN") {
      this->mintheta = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal polar production angle for the leptons" << this->mintheta << std::endl;
#endif
    }
    else if (key=="THMX") {
      this->maxtheta = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal polar production angle for the leptons" << this->maxtheta << std::endl;
#endif
    }
    else if (key=="Q2MN") {
      this->minq2 = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal Q^2 for the incoming photons" << this->minq2 << " GeV^2" << std::endl;
#endif
    }
    else if (key=="Q2MX") {
      this->maxq2 = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal Q^2 for the incoming photons" << this->maxq2 << " GeV^2" << std::endl;
#endif
    }
    else if (key=="MXMN") {
      this->minmx = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal invariant mass of proton remnants" << this->minmx << " GeV/c^2" << std::endl;
#endif
    }
    else if (key=="MXMX") {
      this->maxmx = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal invariant mass of proton remnants" << this->maxmx << " GeV/c^2" << std::endl;
#endif
    }
    else if (key=="GPDF") {
      this->gpdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * GPDF" << this->gpdf << std::endl;
#endif
    }
    else if (key=="SPDF") {
      this->spdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * SPDF" << this->spdf << std::endl;
#endif
    }
    else if (key=="QPDF") {
      this->qpdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * QPDF" << this->qpdf << std::endl;
#endif
    }
    else {
      std::cout << std::setw(60) << "[InputParameters::ReadConfigFile] <WARNING> Unrecognized argument : [" << key << "] = " << value << std::endl;
    }
  }
  f.close();
  std::cout << std::right << "======================================================" << std::endl;
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
  if (this->itvg!=-1) f << "ITVG  " << this->itvg << std::endl;
  if (this->minenergy!=-1) f << "ECUT  " << this->minenergy << std::endl;
  if (this->minenergy!=-1) f << "PTCT  " << this->minpt << std::endl;
  if (this->ntreat!=-1) f << "NTRT  " << this->ntreat << std::endl;
  // ...
  f.close();
  return true;
}

