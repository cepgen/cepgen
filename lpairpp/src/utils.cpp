#include "utils.h"

double GetMassFromPDGId(int pdgId_)
{
  switch (abs(pdgId_)) {
  case 1:    return 0.33;           // d (from PYTHIA6.4)
  case 2:    return 0.33;           // u (from PYTHIA6.4)
  case 11:   return 0.510998928e-3; // electron
  case 13:   return 0.1056583715;   // muon
  case 15:   return 1.77682;        // tau
  case 21:   return 0.;             // gluon
  case 22:   return 0.;             // photon
  case 211:  return 0.13957018;     // pi+
  case 111:  return 0.1349766;      // pi0
  case 553:  return 20.;            // J/psi //FIXME FIXME FIXME
  case 2101: return 0.57933;        // (ud)0 (from PYTHIA6.4)
  case 2103: return 0.77133;        // (ud)1 (from PYTHIA6.4)
  case 2203: return 0.77133;        // (uu)1 (from PYTHIA6.4)
  case 2212: return 0.938272046;    // proton
  default:   return -1.;
  }
}

double GetWidthFromPDGId(int pdgId_)
{
  switch (abs(pdgId_)) {
  case 553:  return 10.; //FIXME
  default:   return -1.;
  }
}

double GetBRFromPDGId(int pdgId_)
{
  switch (abs(pdgId_)) {
  case 113:   return 1.0;    // rho0->pi+ pi-
  case 223:   return 0.0221; // omega->pi+ pi-
  case 333:   return 0.491;  // phi->K+ K-
  case 3332:  return 0.344;  // phi->KL0 KS0 //FIXME FIXME FIXME
  case 444:   return 0.0598; // J/psi->l+ l-
  case 20443: return 0.0425; // psi'->l+ l- X
  case 553:   return 0.0250; // Upsilon(1s)->l+ l-
  case 20553: return 0.0200; // Upsilon(2s)->l+ l- X
  case 30553: return 0.0217; // Upsilon(3s)->l+ l- X
    //case 40113: // rho(1450)->pi+ pi- rho0
    //case 10333: // phi(1680)->K Kbar
  default: return -1;
  }
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

double RanBW(double er_, double gamma_, double emin_, double emax_)
{
  double a, b, e;

  if (gamma_<1.e-3*er_) {
    return er_;
  }
  a = atan(2.*(emax_-er_)/gamma_);
  b = atan(2.*(emin_-er_)/gamma_);
  e = er_+gamma_*tan(drand()*(a-b)+b)/2.;
  if (e<emax_) {
    return e;
  }
  return emax_;
}

double GenerT(double tmin_, double tmax_, double b_, double anexp_)
{
  double c0, c1, bloc, z, t;
  int iter;

  // Generate spectrum by method of R. Lausen

  bloc = b_;
  if (b_<.1) {
    std::cerr << "[GenerT] ERROR: B=" << b_ << std::endl;
    //CALL ERRLOG (20, 'W: GENERT: B < 0.1')
    bloc = .1;
  }
  if (tmin_>=tmax_) {
    std::cerr << "[GenerT] ERROR: TMIN=" << tmin_ << ", TMAX=" << tmax_ << " => return TMIN=" << tmin_ << std::endl;
    //CALL ERRLOG (21, 'S: GENERT: TMIN >= TMAX')
    return tmin_;
  }

  iter = 0;
  do {
    if (anexp_<=1.) {
      // power law exponent is 0 or illegal                                                                                                                      
      //  => generate pure exp(bt) spectrum 
      if (bloc*(tmax_-tmin_)>=25.) {
	t = tmin_-log(drand())/bloc;
#ifdef DEBUG
	std::cout << "[GenerT] DEBUG: Method 1: T=" << t << std::endl;
#endif
      }
      else {
	t = tmin_-log(1.-drand()*(1.-exp(bloc*(tmin_-tmax_))))/bloc;
#ifdef DEBUG
	std::cout << "[GenerT] DEBUG: Method 2: T=" << t << std::endl;
#endif
      }
    }
    else {
      // New 16.5.07 BL:
      // Generate mixed exp(bt)/power law spectrum
      // d sigma/d t = exp (-n*ln(-bt/n+1)) = (-bt/n+1)^-n
      // Limit for small bt: exp (bt + c t^2) with c=b^2/2n
      // Limit for large bt>>n: t^-n
      c1 = std::pow(anexp_+bloc*tmin_, 1.-anexp_);
      c0 = std::pow(anexp_+bloc*tmax_, 1.-anexp_);
      z = drand();
      t = -(anexp_-std::pow(z*(c1-c0)+c0, 1./(1.-anexp_)))/bloc;
    }
    iter++;
  } while ((t<tmin_ or t>tmax_) and iter<=100);
  if (iter>100) {
    //CALL ERRLOG (22, 'W: GENT: More than 100 iterations!')
    std::cout << "[GenerT] WARNING: more than 100 iterations!" << std::endl
	      << "TMIN: " << tmin_ << ", TMAX: " << tmax_ << " BLOC: " << bloc << ", T: " << t
	      << std::endl;
  }
  return t;
}

double GenTDL(double tmin_, double tmax_, double b_, int n_)
{
  int iter;
  double t, w;

  if (tmin_>tmax_) {
    std::cerr << "[GenTDL] ERROR: TMIN=" << tmin_ << ", TMAX=" << tmax_ << " => return TMIN=" << tmin_ << std::endl;
    return tmin_;
  }

  iter = 0;
  do {
    if (b_*(tmax_-tmin_)>=25.) {
      t = tmin_-log(drand())/b_;
#ifdef DEBUG
      std::cout << "[GenTDL] DEBUG: Method 1: T=" << t << std::endl;
#endif
    }
    else {
      t = tmin_-log(1.-drand()*(1.-exp(b_*(tmin_-tmax_))))/b_;
#ifdef DEBUG
      std::cout << "[GenTDL] DEBUG: Method 2: T=" << t << std::endl;
#endif
    }
    w = std::pow((1.+1.41*tmin_)/(1.+1.41*t), n_);
    iter += 1;
  } while ((t<tmin_ or t>tmax_ or w<drand()) and iter<=100);
  if (iter>100) {
    //CALL ERRLOG (22, 'W: GENTDL: More than 100 iterations!')
    std::cout << "[GenTDL] WARNING: more than 100 iterations!" << std::endl
	      << "TMIN: " << tmin_ << ", TMAX: " << tmax_ << ", T: " << t
	      << std::endl;
  }
  return t;
}

int Heli(double longFr_)
{
  if (drand()<longFr_) return 0; // longitudinal photon
  else if (drand()<.5) return 1; // transverse photon
  else return -1;
}

