#include "utils.h"

void Map(double expo_, double xmin_, double xmax_, double* out_, double* dout_, const std::string& var_name_)
{
  double y, out;
  y = xmax_/xmin_;
  out = xmin_*std::pow(y, expo_);
  *out_ = out;
  *dout_ = out*log(y);
  DebuggingInsideLoop(Form("Mapping variable \"%s\"\n\t"
                           "min = %f\n\tmax = %f\n\tmax/min = %f\n\t"
                           "exponent = %f\n\t"
                           "output = %f\n\td(output) = %f",
                           var_name_.c_str(), xmin_, xmax_, y, expo_, *out_, *dout_));
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

double BreitWigner(double er, double gamma, double emin, double emax, double x)
{
  if (x==-1.) x = drand();
  if (gamma<1.e-3*er) { return er; }

  const double a = atan(2.*(emax-er)/gamma),
               b = atan(2.*(emin-er)/gamma),
               e = er+gamma*tan(x*(a-b)+b)/2.;

  if (e>emax) { return emax; }
  return e;
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

