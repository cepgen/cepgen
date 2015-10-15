#include "utils.h"

void Map(double expo_, double xmin_, double xmax_, double* out_, double* dout_)
{
  double y, out;
  y = xmax_/xmin_;
  out = xmin_*std::pow(y, expo_);
  *out_ = out;
  *dout_ = out*log(y);
  DebugInsideLoop(Form("min = %f\n\t"
                            "max = %f\n\t"
                            "max/min = %f\n\t"
                            "exponent = %f\n\t"
                            "output = %f\n\t"
                            "d(output) = %f",
                            xmin_, xmax_, y, expo_, *out_, *dout_));
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
    Info(Form("ERROR: B=%f < 1", b_));
    bloc = .1;
  }
  if (tmin_>=tmax_) {
    Info(Form("ERROR: TMIN=%f >= TMAX=%f => return TMIN=%f", tmin_, tmax_, tmin_))
    return tmin_;
  }

  iter = 0;
  do {
    if (anexp_<=1.) {
      // power law exponent is 0 or illegal
      //  => generate pure exp(bt) spectrum 
      if (bloc*(tmax_-tmin_)>=25.) {
	      t = tmin_-log(drand())/bloc;
        Debug(Form("Method 1: T=%f", t));
      }
      else {
	      t = tmin_-log(1.-drand()*(1.-exp(bloc*(tmin_-tmax_))))/bloc;
        Debug(Form("Method 2: T=%f", t));
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
    Info(Form("WARNING: more than 100 iterations!\n\t"
                   "TMIN: %f, TMAX: %f, BLOC: %f, T: %f", tmin_, tmax_, bloc, t));
  }
  return t;
}

double GenTDL(double tmin_, double tmax_, double b_, int n_)
{
  int iter;
  double t, w;

  if (tmin_>tmax_) {
    Info(Form("ERROR: TMIN=%f, TMAX=%f => return TMIN=%f", tmin_, tmax_, tmin_));
    return tmin_;
  }

  iter = 0;
  do {
    if (b_*(tmax_-tmin_)>=25.) {
      t = tmin_-log(drand())/b_;
      Debug(Form("Method 1: T=%f", t));
    }
    else {
      t = tmin_-log(1.-drand()*(1.-exp(b_*(tmin_-tmax_))))/b_;
      Debug(Form("Method 2: T=%f", t));
    }
    w = std::pow((1.+1.41*tmin_)/(1.+1.41*t), n_);
    iter += 1;
  } while ((t<tmin_ or t>tmax_ or w<drand()) and iter<=100);
  if (iter>100) {
    Info(Form("WARNING: more than 100 iterations!\n\t"
                   "TMIN: %f, TMAX: %f, T: %f", tmin_, tmax_, t));
  }
  return t;
}

int Heli(double longFr_)
{
  if (drand()<longFr_) return 0; // longitudinal photon
  else if (drand()<.5) return 1; // transverse photon
  else return -1;
}

double ThetaToEta(double theta_)
{
  return -log(tan(theta_/180.*pi/2.));
}

double EtaToTheta(double eta_)
{
  return 2.*atan(exp(-eta_))*180./pi;
}

double EtaToY(double eta_, double m_, double pt_)
{
  const double mt = pow(m_, 2)+pow(pt_, 2);
  return asinh(sqrt((((pow(mt, 2)-pow(m_, 2))*cosh(2*eta_)+pow(m_, 2))/pow(mt, 2)-1.)/2.));
}

std::string
Form(const std::string fmt, ...)
{  
  int size = ((int)fmt.size()) * 2 + 50;   // Use a rubric appropriate for your code
  std::string str;
  va_list ap;
  while (true) {     // Maximum two passes on a POSIX system...
    str.resize(size);
    va_start(ap, fmt);
    int n = vsnprintf((char*)str.data(), size, fmt.c_str(), ap);
    va_end(ap);
    if (n>-1 and n<size) {  // Everything worked
      str.resize(n);
      return str;
    }
    if (n>-1)  // Needed size returned
         size = n + 1;   // For null char
    else size *= 2;      // Guess at a larger size (OS specific)
  }
  return str;
}  

