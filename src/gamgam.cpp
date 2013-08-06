#include "gamgam.h"

GamGam::GamGam(const unsigned int ndim_, double q2min_, double q2max_, int nOpt_, double x_[]) :
  _ep1(-1), _w1(-1),
  _ep2(-1), _w2(-1.),
  _w3(-1.), _w4(-1.), _w5(-1.), _w6(-1.), _w7(-1.),
  _sqs(-1.),
  _p12(0.), _p13(0.), _p14(0.), _p15(0.), _p23(0.),
  _p24(0.), _p25(0.), _p34(0.), _p35(0.), _p45(0.),
  _p1k2(0.), _p2k1(0.),
  setp1(false), setp2(false), setp3(false), setp5(false),
  setll(false), setin(false), setout(false), setkin(false),
  _wmin(0.), _wmax(-1.),
  _cotth1(-99999.), _cotth2(99999.)
{
  _ndim = ndim_;
  _q2min = q2min_;
  _q2max = q2max_;
  _qp2min = -q2min_;
  _qp2max = -q2max_;
  _nOpt = nOpt_;
  std::cout << std::setprecision(16);
  
  //srand(time(0));
  
  //FIXME ???
  _w12 = 0.;
  _w31 = 0.;
  _w52 = 0.;

  _x = new double[ndim_];
  std::copy(x_, x_+ndim_, _x);
    
  _part = new std::map<int,Particle>;

#ifdef DEBUG
  std::cout << "[GamGam::GamGam] [DEBUG] number of integration parameters : " << ndim_ << std::endl;
  for (int i=0; i<ndim_; i++) {
    std::cout << "  _x[" << i << "] = " << _x[i] << std::endl;
  }
#endif
}

GamGam::~GamGam()
{
  delete[] _x;
  delete _part;
}

bool GamGam::SetOutgoingParticles(int part_, int pdgId_)
{
  double mass_;

  mass_ = GetMassFromPDGId(pdgId_);
  if (mass_<0) {
    return false;
  }

  switch(part_) {
  case 3:
    // First outgoing proton (or remnant)
    _mp3 = mass_;
    _w3 = std::pow(_mp3, 2);
    _pdg3 = pdgId_;
    setp3 = true;
    break;
  case 5:
    // Second outgoing proton (or remnant)
    _mp5 = mass_;
    _w5 = std::pow(_mp5, 2);
    _pdg5 = pdgId_;
    setp5 = true;
    break;
  case 6:
  case 7:
    // First outgoing lepton
    _ml6 = mass_;
    _w6 = std::pow(_ml6, 2);
    _pdg6 = pdgId_;
    // Second outgoing lepton
    _ml7 = mass_;
    _w7 = std::pow(_ml7, 2);
    _pdg7 = pdgId_;
    setll = true;
    break;
  default:
    return false;
  }
  setout = setp3 && setp5 && setll;
  setkin = setin && setout;
#ifdef DEBUG
  std::cout << "[GamGam::SetOutgoingParticles] [DEBUG] Particle \"" << part_ << "\" has PDG id " << pdgId_ << std::endl;
  if (setout) {
    std::cout << "  --> Outgoing state is fully set" << std::endl;
  }
  if (setkin) {
    std::cout << "  --> Kinematics is fully set" << std::endl;
  }
#endif
  return true;
}

bool GamGam::SetIncomingKinematics(Particle ip1_, Particle ip2_)
{
  _p3_p1[0] = ip1_.px;
  _p3_p1[1] = ip1_.py;
  _p3_p1[2] = ip1_.pz;
  _ep1 = ip1_.e;
  _mp1 = ip1_.m;
  _pp1 = ip1_.p;
  _pdg1 = ip1_.pdgId;
  _w1 = std::pow(_mp1, 2);

  _p3_p2[0] = ip2_.px;
  _p3_p2[1] = ip2_.py;
  _p3_p2[2] = ip2_.pz;
  _ep2 = ip2_.e;
  _mp2 = ip2_.m;
  _pp2 = ip2_.p;
  _pdg2 = ip2_.pdgId;
  _w2 = std::pow(_mp2, 2);

  _etot = _ep1+_ep2;
  _ptot = std::sqrt(std::pow(_p3_p1[0]-_p3_p2[0], 2)
                   +std::pow(_p3_p1[1]-_p3_p2[1], 2)
                   +std::pow(_p3_p1[2]-_p3_p2[2], 2));
  
  setp1 = setp2 = setin = true;
  setkin = setin && setout;
  return setp1;
}

bool GamGam::SetIncomingKinematics(int part_, double momentum_[3], int pdgId_)
{
  std::string name;
  double mass, energy, mom2tmp;

  mass = GetMassFromPDGId(pdgId_);

  // Compute the particle energy
  // using its mass and momentum
  mom2tmp = 0.;
  for (int i=0; i<3; i++) {
    mom2tmp += std::pow(momentum_[i],2);
  }
  energy = std::sqrt(mom2tmp+std::pow(mass,2));

  switch(part_) {
  case 1:
    name = "incoming particle 1";
    std::copy(momentum_, momentum_+3, _p3_p1);
    _pp1 = std::sqrt(mom2tmp);
    _ep1 = energy;
    _mp1 = mass;
    _w1 = std::pow(_mp1, 2);
    _pdg1 = pdgId_;
    setp1 = true;
    break;
  case 2:
    name = "incoming particle 2";
    std::copy(momentum_, momentum_+3, _p3_p2);
    _pp2 = std::sqrt(mom2tmp);
    _ep2 = energy;
    _mp2 = mass;
    _w2 = std::pow(_mp2, 2);
    _pdg2 = pdgId_;
    setp2 = true;
    break;
  default:
    return false;
  }
  setin = setp1 && setp2;

  if (setin) { // if the incoming kinematics is fully specified
    _etot = _ep1+_ep2;
    _ptot = std::sqrt(std::pow(_p3_p1[0]-_p3_p2[0], 2)
                     +std::pow(_p3_p1[1]-_p3_p2[1], 2)
                     +std::pow(_p3_p1[2]-_p3_p2[2], 2));
  }

  setkin = setin && setout;
#ifdef DEBUG
  std::cout << "[GamGam::SetIncomingKinematics] [DEBUG] \"" << name_ << "\" (PDG id " << pdgId_ << ") set to\n"
            << "  E = " << energy_ << " GeV,\n"
            << "  p = (" << momentum_[0] << ", " << momentum_[1] << ", " << momentum_[2] << ") GeV/c,\n"
            << "  M = " << mass_ << " GeV/cc"
            << std::endl;
#endif
  return true;
}

bool GamGam::Pickin()
{
  double sig, sig1, sig2; // sig1 = sigma and sig2 = sigma' in [1]
  double sp, ss, st;
  double sb, sd, se;
  double splus, s2x;
  double s2min, s2max;
  double smax, ds2;
  double s1p, s1m, s1pp, s1pm, s2p;
  double sl2, sl3, sl4, sl5, sl6, sl7;
  double t1min, t1max, dt1;
  double t2min, t2max, dt2;
  double t13, t25;
  double rl1, rl2, rl4;
  double r1, r2, r3, r4;
  double b, c;
  double ap;
  double yy4;
  double dd, delb;
  double sbb, sdd, see, ssb, ssd, sse;

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] _nOpt = " << _nOpt << std::endl;
#endif
  _dj = 0.;

  _w4 = std::pow(_mc4, 2);

  sig = _mc4+_mp5;
  sig1 = std::pow(sig, 2);
  sig2 = std::pow(sig, 2);

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] mc4 = " << _mc4 << std::endl;
  std::cout << "[GamGam::Pickin] [DEBUG] sig1 = " << sig1 << std::endl;
  std::cout << "[GamGam::Pickin] [DEBUG] sig2 = " << sig2 << std::endl;
#endif

  // Mass difference between the first outgoing particle and the first incoming
  // particle
  _d1 = _w3-_w1;
  // Mass difference between the second outgoing particle and the second
  // incoming particle
  _d2 = _w5-_w2;
  // Mass difference between the two incoming particles
  _d5 = _w1-_w2;
  // Mass difference between the central two-photons system and the second
  // outgoing particle
  _d6 = _w4-_w5;

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG]"
            << "\n\tw1 = " << _w1
            << "\n\tw2 = " << _w2
            << "\n\tw3 = " << _w3
            << "\n\tw4 = " << _w4
            << "\n\tw5 = " << _w5
            << std::endl;
#endif

  ss = _s+_d5;
  rl1 = std::pow(ss, 2)-4.*_w1*_s; // lambda(s, m1**2, m2**2)

  if (rl1<=0.) {
    std::cerr << "[GamGam::Pickin] [ERROR]"
              << "\n  rl1 = " << rl1 << " <= 0" << std::endl;
    return false;
  }
  _sl1 = std::sqrt(rl1);

  if (_nOpt==0) {
    smax = _s+_w3-2.*_mp3*_sqs;
    //std::cout << "smax = " << smax << std::endl;
    Map(_x[2], sig1, smax, &_s2, &ds2);
    //std::cout << "After mapping : s2 = " << _s2 << ", ds2 = " << ds2 << std::endl;
    sig1 = _s2; //FIXME!!!!!!!!!!!!!!!!!!!!
  }
#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] _s2 = " << _s2 << std::endl;
#endif

  sp = _s+_w3-sig1;

  _d3 = sig1-_w2;
  //std::cout << "sig1 = " << sig1 << ", w2 = " << _w2 << ", d3 = " << _d3 << std::endl;

  rl2 = std::pow(sp, 2)-4.*_s*_w3; // lambda(s, m3**2, sigma)
  if (rl2<=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]"
              << "\n  rl2 = " << rl2 << " <= 0" << std::endl;
#endif
    return false;
  }
  sl2 = std::sqrt(rl2);

  //std::cout << " ==> sl2 = " << sl2 << std::endl;

  t1max = _w1+_w3-(ss*sp+_sl1*sl2)/(2.*_s); // definition from eq. (A.4) in [1]
  t1min = (_d1*_d3+(_d3-_d1)*(_d3*_w1-_d1*_w2)/_s)/t1max; // definition from eq. (A.5) in [1]

  /*std::cout << "first definition of t1min = " << t1min << std::endl;
  std::cout << "d1, d3 = " << _d1 << ", " << _d3 << std::endl;
  std::cout << "ss, sp = " << ss << ", " << sp << std::endl;
  std::cout << "sl1, sl2 = " << _sl1 << ", " << sl2 << std::endl;
  std::cout << "t1max, s, w2 = " << t1max << ", " << _s << ", " << _w2 << std::endl;*/
  //std::cout << "[before the thresholds] t1min = " << t1min << ", t1max = " << t1max << std::endl;

  // FIXME dropped in CDF version
  if (t1max>_qp2min or t1min<_qp2max) {
    std::cerr << "[GamGam::Pickin] [ERROR]"
              << "\n  t1max>_qp2min or t1min<_qp2max :"
              << "\n  t1max = " << t1max
              << "\n  t1min = " << t1min
              << std::endl;
    return false;
  }
  if (t1max<_qp2max) {
    t1max = _qp2max;
  }
  if (t1min>_qp2min) {
    t1min = _qp2min;
  }
  /////
  //std::cout << "[after the thresholds] t1min = " << t1min << ", t1max = " << t1max << std::endl;

  // t1, the first photon propagator, is defined here
  Map(_x[0], t1min, t1max, &_t1, &dt1);
  // changes wrt mapt1 : dx->-dx
  dt1 = -dt1;
#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] definition of t1 according to"
            << " (t1min, t1max) = (" << t1min << ", " << t1max << ")"
            << "\n  _t1 = " << _t1
            << std::endl;
#endif

  _d4 = _w4-_t1;
  _d8 = _t1-_w2;

  //std::cout << "d1 = " << _d1 << ", d5 = " << _d5 << ", d8 = " << _d8 << std::endl;

  t13 = _t1-_w1-_w3;

  //std::cout << "===> t13 = " << t13 << std::endl;

  _sa1 =-std::pow(_t1-_d1, 2)/4.+_w1*_t1;
  if (_sa1>=0.) {
    std::cerr << "[GamGam::Pickin] [ERROR]"
              << "\n  _sa1>=0 : " << _sa1
              << std::endl;
    return false;
  }

  sl3 = std::sqrt(-_sa1);

  // one computes splus and (s2x=s2max)
  if (_w1!=0.) {
    sb =(_s*(_t1-_d1)+_d5*t13)/(2.*_w1)+_w3;
    sd = _sl1*sl3/_w1;
    se =(_s*(_t1*(_s+t13-_w2)-_w2*_d1)+_w3*(_d5*_d8+_w2*_w3))/_w1;
    if (fabs((sb-sd)/sd)>=1.) {
      splus = sb-sd;
      s2max = se/splus;
    }
    else {
      s2max = sb+sd;
      splus = se/s2max;
    }
  }
  else {
    // 3
    s2max = (_s*(_t1*(_s+_d8-_w3)-_w2*_w3)+_w2*_w3*(_w2+_w3-_t1))/(ss*t13);
    splus = sig2;
  }
  // 4
  s2x = s2max;
#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] s2x = s2max = " << s2x << std::endl;
#endif

  if (_nOpt<0) { // 5
    if (splus>sig2) {
      sig2 = splus;
#ifdef DEBUG
      std::cout << "[GamGam::Pickin] [DEBUG] sig2 truncated to splus = " << splus << std::endl;
#endif
    }
    if (_nOpt<-1) {
      Map(_x[2], sig2, s2max, &_s2, &ds2);
    }
    else { // _nOpt==-1
      Mapla(_t1, _w2, _x[2], sig2, s2max, &_s2, &ds2);
    }
    s2x = _s2;
  }
  else if (_nOpt==0) { // 6
    s2x = _s2;
  }

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] s2x = " << s2x << std::endl;
#endif
  // 7
  r1 = s2x-_d8;
  r2 = s2x-_d6;
  //std::cout << (std::pow(r1, 2)-4.*_w2*s2x) << ", s2x = " << std::pow(r2, 2) << " <-> " << -4.*_w5*s2x << ", " << (4.*_w2*s2x) << std::endl;
  rl4 = (std::pow(r1, 2)-4.*_w2*s2x)*(std::pow(r2, 2)-4.*_w5*s2x);
  if (rl4<=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]\n  rl4<=0 : " << rl4 << std::endl;
#endif
    return false;
  }
  sl4 = std::sqrt(rl4);
  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  t2max = _w2+_w5-(r1*r2+sl4)/(2.*s2x);
  t2min = (_d2*_d4+(_d4-_d2)*(_d4*_w2-_d2*_t1)/s2x)/t2max;

  // t2, the second photon propagator, is defined here
  Map(_x[1], t2min, t2max, &_t2, &dt2);
  // changes wrt mapt2 : dx->-dx
  //std::cout << "_t2 = " << _t2 << std::endl;

  dt2 = -dt2;

  _d7 = _t1-_t2;
  r3 = _d4-_t2;
  r4 = _d2-_t2;

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG]\n"
            << "  r1 = " << r1 << "\n"
            << "  r2 = " << r2 << "\n"
            << "  r3 = " << r3 << "\n"
            << "  r4 = " << r4 << std::endl;
#endif

  b = r3*r4-2.*(_t1+_w2)*_t2;
  c = _t2*_d6*_d8+(_d6-_d8)*(_d6*_w2-_d8*_w5);
  t25 = _t2-_w2-_w5;

  _sa2 = -std::pow(r4, 2)/4.+_w2*_t2;
  if (_sa2>=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]\n  _sa2 = " << _sa2 << " >= 0" << std::endl;
#endif
    return false;
  }
  sl6 = 2.*std::sqrt(-_sa2);
  _g4 = -std::pow(r3, 2)/4.+_t1*_t2;
  if (_g4>=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]\n  _g4 = " << _g4 << " >= 0" << std::endl;
#endif
    return false;
  }
  sl7 = std::sqrt(-_g4)*2.;
  sl5 = sl6*sl7;
  if (fabs((sl5-b)/sl5)>=1.) {
    s2p = (sl5-b)/(2.*_t2);
    s2min = c/(_t2*s2p);
  }
  else { // 8
    s2min = (-sl5-b)/(2.*_t2);
    s2p = c/(_t2*s2min);
  }
  // 9
  if (_nOpt>1) {
    Map(_x[2], s2min, s2max, &_s2, &ds2);
  }
  else if (_nOpt==1) { // _nOpt=1
    Mapla(_t1, _w2, _x[2], s2min, s2max, &_s2, &ds2);
  }
  ap = -std::pow(_s2+_d8, 2)/4.+_s2*_t1;
  if (_w1!=0.) {
    _dd1 = -_w1*(_s2-s2max)*(_s2-splus)/4.;
  }
  else { // 10
    _dd1 = ss*t13*(_s2-s2max)/4.;
  }
  // 11
  _dd2 = -_t2*(_s2-s2p)*(_s2-s2min)/4.;

  // FIXME there should be a more beautiful way to check for nan!
  // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
  // FIXME dropped in CDF version
  if (!(_dd2>=0.) && !(_dd2<0.)) { // NaN
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR] : dd2 == NaN" << std::endl;
    std::cerr << "  dd2 = " << _dd2 << "\n"
              << "  s2 = " << _s2 << "\n"
              << "  s2p = " << s2p << "\n"
              << "  s2min = " << s2min << "\n"
              << "  t2min = " << t2min << "\n"
              << "  t2max = " << t2max
              << std::endl;
#endif
  }
  /////
  if (_x[3]>1. or _x[3]<-1.) {
    std::cerr << "[GamGam::Pickin] [ERROR] x[3] = " << _x[3] << std::endl;
    return false;
  }
  yy4 = cos(pi*_x[3]);
  dd = _dd1*_dd2;
  _p12 = (_s-_w1-_w2)/2.;
  st = _s2-_t1-_w2;
  delb = (2.*_w2*r3+r4*st)*(4.*_p12*_t1-(_t1-_d1)*st)/(16.*ap);

  if (dd<=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]\n  dd = " << dd << " <= 0" << std::endl;
#endif
    return false;
  }

  _delta = delb-yy4*st*std::sqrt(dd)/(2.*ap);
  _s1 = _t2+_w1+(2.*_p12*r3-4.*_delta)/st;

  //std::cout << "pickin's _s1 = " << _s1 << std::endl;

  if (ap>=0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Pickin] [ERROR]\n  ap = " << ap << " >= 0" << std::endl;
#endif
    return false;
  }

  _dj = ds2*dt1*dt2*std::pow(pi, 2)/(8.*_sl1*std::sqrt(-ap));

#ifdef DEBUG
  std::cout << "[GamGam::Pickin] [DEBUG] _dj = " << _dj << std::endl;
#endif

  _gram = (1.-std::pow(yy4, 2))*dd/ap;

  _p13 = -t13/2.;
  _p14 = (_d7+_s1-_w3)/2.;
  _p15 = (_s+_t2-_s1-_w2)/2.;
  _p23 = (_s+_t1-_s2-_w1)/2.;
  _p24 = (_s2-_d7-_w5)/2.;
  _p25 = -t25/2.;
  _p34 = (_s1-_w3-_w4)/2.;
  _p35 = (_s+_w4-_s1-_s2)/2.;
  _p45 = (_s2-_w4-_w5)/2.;

  _p1k2 = (_s1-_t2-_w1)/2.;
  _p2k1 = st/2.;


  if (_w2!=0.) {
    sbb = (_s*(_t2-_d2)-_d5*t25)/(2.*_w2)+_w5;
    sdd = _sl1*sl6/(2.*_w2);
    see = (_s*(_t2*(_s+t25-_w1)-_w1*_d2)+_w5*(_w1*_w5-_d5*(_t2-_w1)))/_w2;
    if (sbb/sdd>=0.) {
      s1p = sbb+sdd;
      s1m = see/s1p;
      // FIXME there should be a more beautiful way to check for nan!
      // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
      // FIXME dropped in CDF version
      if (!(_dd2>=0.) && !(_dd2<0.)) { // NaN
#ifdef ERROR
        std::cout << "[GamGam::Pickin] [ERROR] : dd2 == NaN" << std::endl;
        std::cout << "  dd2 = " << _dd2 << "\n  s1 = " << _s1 << "\n"
                  << "  s1p = " << s1p << "\n  s1m = " << s1m << "\n"
                  << "  w2 = " << _w2
                  << std::endl;
#endif
      }
      /////
    }
    else { // 12
      s1m = sbb-sdd;
      s1p = see/s1m;
    }
    _dd3 = -_w2*(s1p-_s1)*(s1m-_s1)/4.; // 13
    //std::cout << "--> (13) : s1m=" << s1m << ", s1p=" << s1p << ", dd3=" << _dd3 << std::endl;
  }
  else { // 14
    s1p = (_s*(_t2*(_s-_w5+_t2-_w1)-_w1*_w5)+_w1*_w5*(_w1+_w5-_t2))/(t25*(_s-_d5));
    _dd3 = -t25*(_s-_d5)*(s1p-_s1)/4.;
  }
  // 15
  _acc3 = (s1p-_s1)/(s1p+_s1);

  ssb = _t2+_w1-r3*(_d1-_t1)/(2.*_t1);
  ssd = sl3*sl7/_t1;
  sse = (_t2-_w1)*(_w4-_w3)+(_t2-_w4+_d1)*((_t2-_w1)*_w3-(_w4-_w3)*_w1)/_t1;

  if (ssb/ssd>=0.) {
    s1pp = ssb+ssd;
    s1pm = sse/s1pp;
  }
  else { // 16
    s1pm = ssb-ssd;
    s1pp = sse/s1pm;
  }
  // 17
  _dd4 = -_t1*(_s1-s1pp)*(_s1-s1pm)/4.;
  _acc4 = (_s1-s1pm)/(_s1+s1pm);
  _dd5 = _dd1+_dd3+((_p12*(_t1-_d1)/2.-_w1*_p2k1)*(_p2k1*(_t2-_d2)-_w2*r3)-_delta*(2.*_p12*_p2k1-_w2*(_t1-_d1)))/_p2k1;

  return true;
}

bool GamGam::Orient()
{
  double re, rr, a1;

  bool pck = Pickin();
  if (!pck || _dj==0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR] Pickin() failed : (bool)" << pck << ", _dj = " << _dj << std::endl;
#endif
    return false;
  }
  re = 1./(2.*_sqs);
  _ep1 = re*(_s+_w12);
  _ep2 = re*(_s-_w12);

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG]\n  re = " << re << "\n  _w12 = " << _w12 << std::endl;
  std::cout << "[GamGam::Orient] [DEBUG] incoming particles' energy = " << _ep1 << ", " << _ep2 << std::endl;
#endif

  _p = re*_sl1;

  _de3 = re*(_s2-_w3+_w12);
  _de5 = re*(_s1-_w5-_w12);

  // Final state energies
  _ep3 = _ep1-_de3;
  _ec4 = _de3+_de5;
  _ep5 = _ep2-_de5;

  if (_ec4<_mc4) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  _ec4<_mc4 : _ec4 = " << _ec4 << ", _mc4 = " << _mc4 << std::endl;
#endif
    return false;
  }
  // What if the protons' momenta are not along the z-axis?
  _pp3 = std::sqrt(std::pow(_ep3, 2)-_w3);
  _pc4 = std::sqrt((std::pow(_ec4, 2)-std::pow(_mc4, 2)));

  if (_pc4==0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  _pzc4==0" << std::endl;
#endif
    return false;
  }
  _pp5 = std::sqrt(std::pow(_ep5, 2)-_w5);
  _p_p3 = std::sqrt(_dd1/_s)/_p;

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG] central system's energy : E4 = " << _ec4 << std::endl;
  std::cout << "[GamGam::Orient] [DEBUG] central system's momentum : P4 = " << _pc4 << std::endl;
  std::cout << "[GamGam::Orient] [DEBUG] central system's invariant mass : M4 = " << _mc4 << std::endl;
  std::cout << "[GamGam::Orient] [DEBUG] outgoing particles' energy : E3 = " << _ep3 << ", E5 = " << _ep5 << std::endl;
#endif
  //std::cout << "dd3 = " << _dd3 << std::endl;

  _p_p5 = std::sqrt(_dd3/_s)/_p;
  _st3 = _p_p3/_pp3;
  _st5 = _p_p5/_pp5;

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG] _st3 = " << _st3 << ", _st5 = " << _st5 << std::endl;
#endif

  // FIXME there should be a more beautiful way to check for nan!
  // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
  // FIXME dropped in CDF version
  if (!(_dd3>=0.) && !(_dd3<0.)) { // NaN
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  dd3 == NaN" << std::endl;
#endif
  }
  if (!(_dd1>=0.) && !(_dd1<0.)) { // NaN
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  dd1 == NaN" << std::endl;
#endif
  }
  /////

  if (_st3>1. or _st5>1.) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  _st3>1 or _st5>1 : _st3 = " << _st5 << ", _st5 = " << _st5 << std::endl;
#endif
    return false;
  }
  _ct3 = std::sqrt(1.-std::pow(_st3, 2));
  _ct5 = std::sqrt(1.-std::pow(_st5, 2));

  if (_ep1*_ep3<_p13) {
    _ct3 = -_ct3;
  }

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG] _ct3 = " << _ct3 << ", _ct5 = " << _ct5 << std::endl;
#endif

  if (_ep2*_ep5>_p25) {
    _ct5 = -_ct5;
  }
  _al3 = std::pow(_st3, 2)/(1.+_ct3);
  _be5 = std::pow(_st5, 2)/(1.-_ct5);

  if (_dd5<0.) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  _dd5<0 : " << _dd5 << std::endl;
#endif
    return false;
  }

  // Centre of mass system kinematics (theta4 and phi4)
  _p_p4 = std::sqrt(_dd5/_s)/_p;
  _st4 = _p_p4/_pc4;

  if (_st4>1.) {
#ifdef ERROR
    std::cerr << "[GamGam::Orient] [ERROR]\n  _st4>1 : " << _st4 << std::endl;
#endif
    return false;
  }
  _ct4 = std::sqrt(1.-std::pow(_st4, 2));
  if (_ep1*_ec4<_p14) {
    _ct4 =-_ct4;
  }

  _al4 = 1.-_ct4;
  _be4 = 1.+_ct4;


  if (_ct4<0.) {
    _be4 = std::pow(_st4, 2)/_al4;
  }
  else {
    _al4 = std::pow(_st4, 2)/_be4;
  }

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG] _ct4 = " << _ct4 << ", _al4 = " << _al4 << ", _be4 = " << _be4 << std::endl;
#endif

  //std::cout << "pp3 = " << _p_p3 << ", pp5 = " << _p_p5 << std::endl;

  rr  = std::sqrt(-_gram/_s)/(_p*_p_p4);
  _sp3 = rr/_p_p3;
  _sp5 = -rr/_p_p5;
  //std::cout << "rr = " << rr << ", sp3 = " << fabs(_sp3) << ", sp5 = " << fabs(_sp5) << std::endl;

  if (fabs(_sp3)>1. || fabs(_sp5)>1.) {
    std::cerr << "[GamGam::Orient] [ERROR]\n  |_sp3|>1 or |_sp5|>1 :\n"
              << "   |_sp3| = " << fabs(_sp3) << "\n"
              << "  _sp3**2 = " << std::pow(_sp3, 2) << "\n"
              << "   |_sp5| = " << fabs(_sp5) << "\n"
              << "  _sp5**2 = " << std::pow(_sp5, 2) << std::endl;
    return false;
  }

  _cp3 = -std::sqrt(1.-std::pow(_sp3, 2));
  _cp5 = -std::sqrt(1.-std::pow(_sp5, 2));

  a1 = _p_p3*_cp3-_p_p5*_cp5; //OK!!!

#ifdef DEBUG
  std::cout << "[GamGam::Orient] [DEBUG] Kinematic quantities"
            << "\n  cos(theta3) = " << _ct3
            << "\n  sin(theta3) = " << _st3
            << "\n  cos(phi3) = " << _cp3
            << "\n  sin(phi3) = " << _sp3
            << "\n  cos(theta4) = " << _ct4
            << "\n  sin(theta4) = " << _st4
            << "\n  cos(phi4) = " << _ct4
            << "\n  sin(phi4) = " << _st4
            << "\n  cos(theta5) = " << _ct5
            << "\n  sin(theta5) = " << _ct5
            << "\n  cos(phi5) = " << _cp5
            << "\n  sin(phi5) = " << _cp5
            << std::endl;
#endif

  if (fabs(_p_p4+_p_p3*_cp3+_cp5*_p_p5)<fabs(fabs(a1)-_p_p4)) {
#ifdef DEBUG
    std::cout << "[GamGam::Orient] [DEBUG] fabs(_p_p4+_p_p3*_cp3+_cp5*_p_p5)<fabs(fabs(a1)-_p_p4)"
              << "\n  pp4 = " << _p_p4
              << "\n  pp5 = " << _p_p5
              << "\n  cp3 = " << _cp3
              << "\n  cp5 = " << _cp5
              << "\n  a1 = " << a1
              << std::endl;
#endif
    return true;
  }
  if (a1<0.) {
    _cp5 = -_cp5;
  }
  else {
    _cp3 = -_cp3;
  }

  return true;
}

void GamGam::ComputeSqS()
{
  double k;

  k = 0.;
  for (int i=0; i<3; i++) {
    k += _p3_p1[i]*_p3_p2[i];
  }
  _s = std::pow(_mp1,2)+std::pow(_mp2,2)+2.*(_ep1*_ep2+k);
  _sqs = sqrt(_s);

#ifdef DEBUG
  std::cout << "[GamGam::ComputeSqS] [DEBUG] Centre of mass energy : " << _sqs << " GeV" << std::endl;
#endif

}

double GamGam::ComputeXsec(int nm_)
{

  // WARNING ====> PP5-->_p_p5
  //               P5--->_pp5
  double xsec;
  double dw4;
  double wmin, wmax;
  double e1mp1, e3mp3;
  double eg, pg;
  double pgx, pgy, pgz, pgp, pgg;
  double stg, cpg, spg, ctg;
  double xx6;
  double amap, bmap, ymap, beta;
  double ddd;

  double pp6, pp7;
  double p6x, p6y, p6z;
  double pz6, pz7;

  double qcx, qcz;
  double pc6x, pc6z;
  double pcm6x, pcm6y, pcm6z, pcm6, ecm6;

  double phicm6, spcm6, cpcm6;

  double b1, b2, b3;
  double c1, c2, c3;
  double h1, h2, hq;
  double r12, r13, r22, r23;

  double cott6, cott7, cost6, cost7;
  bool lcut, lmu1, lmu2;

  xsec = 0.;

  if (!setout) {
    std::cout << "[GamGam::ComputeXsec] [ERROR]\n  : output state not set !" << std::endl;
    return 0.;
  }
  ComputeSqS();
  //std::cout << "wmax = " << _wmax << std::endl;
  if (_wmax<0) _wmax = _s;

  // the minimal energy for the central system is its outgoing leptons' mass
  // energy (or wmin_ if specified)
  wmin = std::pow(_ml6+_ml7,2); // POW ok

  if (wmin<_wmin) wmin = _wmin; // lower limit on wmin if specified

  // the maximal energy for the central system is its CM energy with the
  // outgoing particles' mass energy substracted (or _wmax if specified)
  wmax = std::pow(_sqs-_mp3-_mp5,2); // POW ok

  if (wmax>_wmax) wmax = _wmax; // upper limit on wmax if specified

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG]"
            << "\n  wmin = " << wmin
            << "\n  wmax = " << wmax
            << "\n  wmax/wmin = " << wmax/wmin
            << std::endl;
#endif
  Map(_x[4],wmin,wmax,&_w4,&dw4);
  _mc4 = std::sqrt(_w4);

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG] Computed value for w4 = " << _w4 << " -> mc4 = " << _mc4 << std::endl;
#endif

  if (!Orient()) {
    return 0.;
  }

  if (_t1>0. or _t2>0.) {
    _dj = 0.;
  }
  if (_dj==0.) {
#ifdef ERROR
    std::cout << "[GamGam::ComputeXsec] [ERROR]\n  _dj = " << _dj << std::endl;
#endif
    return 0.;
  }
  ecm6 = (_w4+_w6-_w7)/(2.*_mc4);
  pcm6 = std::sqrt(std::pow(ecm6, 2)-_w6);
  
  //std::cout << "pcm6 = " << pcm6 << ", ecm6 = " << ecm6 << ", w6 = " << _w6 << std::endl; // OKKKKK!!!!

  _dj *= dw4*pcm6/(_mc4*sconstb*_s);
  
  // Let the most obscure part of this code begin...

  e1mp1 = _w1/(_ep1+_p);
  e3mp3 = _w3/(_ep3+_pp3);
  eg = (_w4+_t1-_t2)/(2.*_mc4);
  pg = std::sqrt(std::pow(eg, 2)-_t1);
  
  pgx = -_p_p3*_cp3*_ct4-_st4*(_de3-e1mp1+e3mp3+_pp3*_al3);
  pgy = -_p_p3*_sp3;
  pgz = _mc4*_de3/(_ec4+_pc4)-_ec4*_de3*_al4/_mc4-_p_p3*_cp3*_ec4*_st4/_mc4+_ec4*_ct4/_mc4*(_pp3*_al3+e3mp3-e1mp1);

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG] pg3 = ("
            << pg3[0] << ", "
            << pg3[1] << ", "
            << pg3[2] << "), pg3**2 = "
            << std::sqrt(std::pow(pg3[0], 2)+std::pow(pg3[1], 2)+std::pow(pg3[2], 2))
            << std::endl;
#endif

  pgp = std::sqrt(std::pow(pgx, 2)+std::pow(pgy, 2)); // outgoing proton (3)'s transverse momentum
  pgg = std::sqrt(std::pow(pgp, 2)+std::pow(pgz, 2)); // outgoing proton (3)'s momentum
  if (pgg>pgp*.9 && pgg>pg) { //FIXME ???
    pg = pgg;
  }
  
  // phi angle
  cpg = pgx/pgp;
  spg = pgy/pgp;
  
  // theta angle
  stg = pgp/pg;
  ctg = std::sqrt(1.-std::pow(stg, 2));
  if (pgz<0.) {
    ctg = -ctg;
  }
  
  xx6 = _x[5];
  
  if (nm_!=0) {
    amap = (_w4-_t1-_t2)/2.;
    bmap = std::sqrt((std::pow(_w4-_t1-_t2, 2)-4.*_t1*_t2)*(1.-4.*_w6/_w4))/2.;
    ymap = (amap+bmap)/(amap-bmap);
    beta = std::pow(ymap, (double)(2.*xx6-1.));
    xx6 = (amap/bmap*(beta-1.)/(beta+1.)+1.)/2.;
    if (xx6>1.) {
      xx6 = 1.;
    }
    if (xx6<0.) {
      xx6 = 0.;
    }
    _ctcm6 = 1.-2.*xx6;
    /*ddd = (std::pow(amap, 2)-std::pow(bmap*_ctcm6, 2))/(amap*bmap)*log(ymap);*/
    ddd = (amap+bmap*_ctcm6)*(amap-bmap*_ctcm6)/amap/bmap*log(ymap);
    _dj *= ddd/2.;
  }
  // 3D rotation of the first outgoing lepton wrt the CM system
  _ctcm6 = 1.-2.*xx6; // cos(theta_cm,6) is between -1 and 1
  _stcm6 = 2.*std::sqrt(xx6*(1.-xx6)); // definition is OK (according to _ctcm6 def)
#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG]\n\tctcm6 = " << _ctcm6 << "\n\tstcm6 = " << _stcm6 << std::endl;
#endif

  phicm6 = 2.*pi*_x[6];

  cpcm6 = cos(phicm6);
  spcm6 = sin(phicm6);

  // First outgoing lepton's 3-momentum in the centre of mass system
  pcm6x = pcm6*_stcm6*cpcm6;
  pcm6y = pcm6*_stcm6*spcm6;
  pcm6z = pcm6*_ctcm6;

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG] p3cm6 = ("
            << pcm6x << ", "
            << pcm6y << ", "
            << pcm6z << ")"
            << std::endl;
#endif
  pc6z = ctg*pcm6z-stg*pcm6x;

  h1 = stg*pcm6z+ctg*pcm6x;
  //std::cout << "h1=" << h1 << std::endl;
  
  pc6x = cpg*h1-spg*pcm6y;
  
  qcx = 2.*pc6x;
  qcz = 2.*pc6z;
  // qcy == QCY is never defined
  //std::cout << "qcx = " << qcx << std::endl;
  //std::cout << "qcz = " << qcz << std::endl; // QCZ ok!!
  
  // First outgoing lepton's 3-momentum
  p6y = cpg*pcm6y+spg*h1;
  _el6 = (_ec4*ecm6+_pc4*pc6z)/_mc4;
  h2 = (_ec4*pc6z+_pc4*ecm6)/_mc4;
  p6x = _ct4*pc6x+_st4*h2;
  p6z = _ct4*h2-_st4*pc6x;
  
  /*std::cout << "p6x = " << p6x << std::endl; // wrong (-0.97171198099796430)
  std::cout << "p6y = " << p6y << std::endl;
  std::cout << "p6z = " << p6z << std::endl;*/

  _qve[0] = _pc4*qcz/_mc4;
  _qve[2] = 2.*p6y;
  hq = _ec4*qcz/_mc4;
  _qve[1] = _ct4*qcx+_st4*hq;
  _qve[3] = _ct4*hq-_st4*qcx;

  _pl6 = std::sqrt(std::pow(_el6, 2)-_w6); // first outgoing lepton's |p|
  /*double test = std::sqrt(std::pow(p6x,2)+std::pow(p6y,2)+std::pow(p6z,2));
  std::cout << (test-_pl6) << "\t_pl6 = " << _pl6 << " <==> " << test << " = norm" << std::endl;*/

  // available energy for the second lepton is the two-photon system's energy
  // with the first lepton's energy removed
  _el7 = _ec4-_el6;
  _pl7 = std::sqrt(std::pow(_el7, 2)-_w7); // second outgoing lepton's |p|

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG] (outgoing kinematics)"
            << "\n\tfirst outgoing lepton : p, E = " << _pl6 << ", " << _el6
            << "\n\tsecond outgoing lepton : p, E = " << _pl7 << ", " << _el7
            << std::endl;
#endif
  double p7x, p7y, p7z;
  // Second outgoing lepton's 3-momentum
  p7x = _p_p4-p6x;
  p7y = -p6y;
  p7z = _pc4*_ct4-p6z;
  // should be P7X=   1.0886407692489968      P7Y= -0.17366271681007375      P7Z=  -32.904157577959587     
  
  /*std::cout << "p7x=" << p7x << std::endl;
  std::cout << "p7y=" << p7y << std::endl;
  std::cout << "p7z=" << p7z << std::endl;*/
  
  pp6 = std::sqrt(std::pow(p6x, 2)+std::pow(p6y, 2)); // pT of first outgoing lepton
  pp7 = std::sqrt(std::pow(p7x, 2)+std::pow(p7y, 2)); // pT of second outgoing lepton

  // First outgoing lepton's kinematics (sin/cos theta/phi)
  _ct6 = p6z/_pl6;
  _st6 = pp6/_pl6;
  _cp6 = p6x/pp6;
  _sp6 = p6y/pp6;
  if (_st6<0) {
    std::cout << "st6<0 : " << _st6 << std::endl;
  }

  /*std::cout << "pp6 = " << pp6 << std::endl;
  std::cout << "pp7 = " << pp7 << std::endl;*/
  // Second outgoing lepton's kinematics (sin/cos theta/phi)
  _ct7 = p7z/_pl7;
  _st7 = pp7/_pl7;
  _cp7 = p7x/pp7;
  _sp7 = p7y/pp7;

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG] (outgoing trajectories)"
            << "\n\tfirst outgoing lepton : cos(theta) = " << _ct6 << ", sin(theta) = " << _st6
            << "\n\tfirst outgoing lepton : cos(phi) = " << _cp6 << ", sin(phi) = " << _sp6
            << "\n\tsecond outgoing lepton : cos(theta) = " << _ct7 << ", sin(theta) = " << _st7
            << "\n\tsecond outgoing lepton : cos(phi) = " << _cp7 << ", sin(phi) = " << _sp7
            << std::endl;
#endif

  _q1dq = eg*(2.*ecm6-_mc4)-2.*pg*pcm6*_ctcm6;
  _q1dq2 = (_w4-_t1-_t2)/2.;

  _bb = _t1*_t2+(_w4*std::pow(_stcm6, 2)+4.*_w6*std::pow(_ctcm6, 2))*std::pow(pg, 2);
  // Q0=QVE[0], QX=QVE[1], QY=QVE[2], QZ=QVE[3]
  c1 = (_qve[1]*_sp3-_qve[2]*_cp3)*_p_p3;
  c2 = (_qve[3]*_ep1-_qve[0]*_p)*_p_p3;
  c3 = (_w31*std::pow(_ep1, 2)+2.*_w1*_de3*_ep1-_w1*std::pow(_de3, 2)+std::pow(_p_p3, 2)*std::pow(_ep1, 2))/(_ep3*_p+_pp3*_ct3*_ep1);
  
  b1 = (_qve[1]*_sp5-_qve[2]*_cp5)*_p_p5;
  b2 = (_qve[3]*_ep2+_qve[0]*_p)*_p_p5;
  b3 = (_w52*std::pow(_ep2, 2)+2.*_w2*_de5*_ep2-_w2*std::pow(_de5, 2)+std::pow(_p_p5*_ep2, 2))/(_ep2*_pp5*_ct5-_ep5*_p);
  
  //std::cout << "b1, b2, b3 = " << b1 << ", " << b2 << ", " << b3 << std::endl;

  r12 = c2*_sp3+_qve[2]*c3;
  r13 = -c2*_cp3-_qve[1]*c3;

#ifdef DEBUG
  std::cout << "[GamGam::ComputeXsec] [DEBUG]" << std::endl;
  for (unsigned int i=0; i<sizeof(_qve)/sizeof(double); i++) {
    std::cout << "  _qve[" << i << "] = " << _qve[i] << std::endl;
  }
#endif

  r22 = b2*_sp5+_qve[2]*b3;
  r23 = -b2*_cp5-_qve[1]*b3;
  
  //std::cout << "r22 = " << r22 << ", r23 = " << r23 << std::endl;

  _epsi = _p12*c1*b1+r12*r22+r13*r23;

  _g5 = _w1*std::pow(c1, 2)+std::pow(r12, 2)+std::pow(r13, 2);
  _g6 = _w2*std::pow(b1, 2)+std::pow(r22, 2)+std::pow(r23, 2);

  _a5 = -(_qve[1]*_cp3+_qve[2]*_sp3)*_p_p3*_p1k2-(_ep1*_qve[0]-_p*_qve[3])*(_cp3*_cp5+_sp3*_sp5)*_p_p3*_p_p5+(_de5*_qve[3]+_qve[0]*(_p+_pp5*_ct5))*c3;
  _a6 = -(_qve[1]*_cp5+_qve[2]*_sp5)*_p_p5*_p2k1-(_ep2*_qve[0]+_p*_qve[3])*(_cp3*_cp5+_sp3*_sp5)*_p_p3*_p_p5+(_de3*_qve[3]-_qve[0]*(_p-_pp3*_ct3))*b3;
  
  ////////////////////////////////////////////////////////////////
  // END of GAMGAM subroutine in the FORTRAN version
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // INFO from f.f
  ////////////////////////////////////////////////////////////////

  _gamma = _etot/_sqs;
  _betgam = _ptot/_sqs;

  if (_cuts.mode==0) {
#ifdef DEBUG
    std::cout << "[GamGam::ComputeXsec] [DEBUG] No cuts applied on the outgoing leptons kinematics !" << std::endl;
#endif
  }
  // Kinematics computation for both muons
  _pt_l6 = _pl6*_st6;
  pz6 = _betgam*_el6+_gamma*_pl6*_ct6;
  _e6lab = _gamma*_el6+_betgam*_pl6*_ct6; //OK!

  _pt_l7 = _pl7*_st7;
  pz7 = _betgam*_el7+_gamma*_pl7*_ct7;  
  _e7lab = _gamma*_el7+_betgam*_pl7*_ct7; //OK!
  
  /*std::cout << "st6 = " << _st6 << std::endl;
  std::cout << "st7 = " << _st7 << std::endl;
  
  std::cout << "pt6 = " << _pt_l6 << ", pz6 = " << pz6 << std::endl;
  std::cout << "pt7 = " << _pt_l7 << ", pz7 = " << pz7 << std::endl;

  std::cout << "e6lab = " << _e6lab << std::endl;
  std::cout << "e7lab = " << _e7lab << std::endl;*/
  /*std::cout << "ct6 = " << _ct6 << ", st6 = " << _st6 << std::endl;
  std::cout << "cp6 = " << _cp6 << ", sp6 = " << _sp6 << std::endl;
  std::cout << "ct7 = " << _ct7 << ", st7 = " << _st7 << std::endl;
  std::cout << "cp7 = " << _cp7 << ", sp7 = " << _sp7 << std::endl;*/
  //FIXME Introduce here the kinematic cuts!!!!!
  // Still to be implemented :
  // - cut on mx for the proton remnants (inelastic case)
  // - ...

  lcut = false; // Event discarded by default
  cott6 = pz6/_pt_l6;
  cott7 = pz7/_pt_l7;
  lmu1 = cott6>=_cotth1
      && cott6<=_cotth2
      && (_pt_l6>=_cuts.ptmin || _cuts.ptmin<=0.)
      && (_pt_l6<=_cuts.ptmax || _cuts.ptmax<=0.)
      && ( _e6lab>=_cuts.emin || _cuts.emin <=0.)
      && ( _e6lab<=_cuts.emax || _cuts.emax <=0.);
  lmu2 = cott7>=_cotth1
      && cott7<=_cotth2
      && (_pt_l7>=_cuts.ptmin || _cuts.ptmin<=0.)
      && (_pt_l7<=_cuts.ptmax || _cuts.ptmax<=0.)
      && ( _e7lab>=_cuts.emin || _cuts.emin <=0.)
      && ( _e7lab<=_cuts.emax || _cuts.emax <=0.);
  switch (_cuts.mode) {
    case 0:
      lcut = true;
      break;
    case 1: // Vermaseren's hypothetical detector cuts
      cost6 = pz6/std::sqrt(std::pow(pz6,2)+std::pow(_pt_l6,2));
      cost7 = pz7/std::sqrt(std::pow(pz7,2)+std::pow(_pt_l7,2));
      lcut = ((fabs(cost6)<=0.75 && _pt_l6>=1.) || (fabs(cost6)<=0.95 && fabs(cost6)>0.75 && fabs(_p3_l6[2])>1.)) &&
             ((fabs(cost7)<=0.75 && _pt_l7>=1.) || (fabs(cost7)<=0.95 && fabs(cost7)>0.75 && fabs(_p3_l7[2])>1.));
      break;
    case 2:
      lcut = lmu1 and lmu2;
      break;
    case 3:
    default:
      lcut = lmu1 or lmu2;
      break;
  }

  // Cut on mass of final hadronic system (MX)
  /*if (_ndim>7 && (std::pow(mx,2)<std::pow(_cuts.mxmin,2) || std::pow(mx,2)>std::pow(_cuts.mxmax,2))) {
    lcut = false;
  }*/
  // Cut on the proton's Q2 (first photon propagator T1)
  if (_t1<_qp2max || _t1>_qp2min) {
    lcut = false;
  }

  if (!lcut) { // Dismiss the cuts-failing events in the xsection computation
    return 0.;
  }
  this->FillKinematics();
  
  int intgp, intge;
  
  switch (_ndim) {
    case 7: // elastic case
      intgp = 2;
      intge = 2;
      xsec = sconst*_dj*this->PeriPP(intgp, intge);
      break;
    case 8: // single-dissociative case
      intgp = 3;
      intge = 2;
      break;
    case 9: // double-dissociative case
      intgp = 3;
      intge = 3;
      break;
    default:
      intgp = 1;
      intge = 1;
      break;
  }
  //std::cout << "dj = " << _dj << std::endl;
  //std::cout << "Elastic x-sec : " << xsec << std::endl;
  return xsec;
}

void GamGam::FillKinematics()
{
  //double gmux, gmuy, gmuw, gmunu;
  double px, py, pz;
  double ranphi, cp, sp;
  int rany, ransign;
  Particle part;
  
  // Needed to parametrise a random rotation around z-axis
  //std::cout << "rand = " << ((double)rand()/RAND_MAX)*2.*pi << std::endl;
  rany = ((double)rand()>=.5*(double)RAND_MAX) ? 1 : -1;
  ransign = ((double)rand()>=.5*(double)RAND_MAX) ? 1 : -1;
  ranphi = ((double)rand()/(double)RAND_MAX)*2.*pi;
  cp = cos(ranphi);
  sp = sin(ranphi);
  /*cp = ((double)rand()/RAND_MAX)*2.-1.;
  sp = sqrt(1-pow(cp, 2));*/
  
  /*std::ofstream of;
  of.open("haha", std::fstream::app);
  of << ranphi << "\t" << rany << "\t" << ransign << std::endl;
  of.close();*/
  /*cp = 1.;
  sp = -1.;*/
  
  // First incoming proton
  part.SetP(0., 0., _gamma*_pp1+_betgam*_ep1, _gamma*_ep1+_betgam*_pp1);
  part.role = 1;
  part.pdgId = _pdg1;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // Second incoming proton
  Particle ip2;
  part.SetP(0., 0., -_gamma*_pp1+_betgam*_ep2, _gamma*_ep2-_betgam*_pp1);
  part.role = 2;
  part.pdgId = _pdg2;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // First outgoing proton
  Particle op1;
  px = _pp3*_st3*_cp3;
  py = _pp3*_st3*_sp3;
  pz = _pp3*_ct3;
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, _gamma*pz  +_betgam*_ep3, _gamma*_ep3+_betgam*pz);
  part.role = 3;
  part.pdgId = _pdg3;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // Second outgoing proton
  Particle op2;
  px = _pp5*_st5*_cp5;
  py = _pp5*_st5*_sp5;
  pz = _pp5*_ct5;
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, _gamma*pz  +_betgam*_ep5, _gamma*_ep5+_betgam*pz);
  part.role = 5;
  part.pdgId = _pdg5;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // First outgoing lepton
  Particle ol1;
  px = _pl6*_st6*_cp6;
  py = _pl6*_st6*_sp6;
  pz = _pl6*_ct6;
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, _gamma*pz+_betgam*_el6, _gamma*_el6+_betgam*(_pl6*_ct6));
  part.role = 6;
  part.pdgId = ransign*abs(_pdg6);
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // Second outgoing lepton
  Particle ol2;
  px = _pl7*_st7*_cp7;
  py = _pl7*_st7*_sp7;
  pz = _pl7*_ct7;
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, _gamma*pz+_betgam*_el7, _gamma*_el7+_betgam*pz);
  part.role = 7;
  part.pdgId = -ransign*abs(_pdg7);
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // First incoming photon
  // Equivalent in LPAIR : PLAB(x, 3)
  Particle ph1;
  px = _p3_p1[0]-_pp3*_st3*_cp3;
  py = _p3_p1[1]-_pp3*_st3*_sp3;
  pz = _p3_p1[2]-_p3_p3[2];
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, pz, _ep1-_ep5);
  part.role = 41;
  part.pdgId = 22;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  // Second incoming photon
  // Equivalent in LPAIR : PLAB(x, 4)
  Particle ph2;
  px = _p3_p2[0]-_pp5*_st5*_cp5;
  py = _p3_p2[1]-_pp5*_st5*_sp5;
  pz = _p3_p2[2]-_p3_p5[2];
  part.SetP(px*cp+rany*py*sp, -px*sp+rany*py*cp, pz, _ep2-_ep5);
  part.role = 42;
  part.pdgId = 22;
  this->_part->insert(std::pair<int,Particle>(part.role, part));
  
  /*gmux = -_t2/(_ep1*_eg2-_pp1*_p3_g2[2])/2.;
  gmuy = (_ep1*_eg2-_pp1*_p3_g2[2])/(_ep2*_eg2+_pp2*_p3_g2[2]);
  gmuw = std::pow(_ep1+_eg2, 2)-std::pow(_pp1+_p3_g2[2], 2);
  if (gmuw>=0.) {
    gmuw = std::sqrt(gmuw);
  }
  else {
    std::cerr << "[GamGam::FillKinematics] [ERROR] W**2 = " << gmuw << " < 0" << std::endl;
    gmuw = 0.;
  }
  gmunu = gmuy*2.*GetMassFromPDGId(2212)/_ep1/_ep2;
#ifdef DEBUG
  std::cout << "[GamGam::FillKinematics] [DEBUG]"
            << "\n\tgmux = " << gmux
            << "\n\tgmuy = " << gmuy
            << "\n\tgmuw = " << gmuw
            //<< "\n\tgmunu = " << gmunu
            << std::endl;
#endif
  */
}

void GamGam::SetCuts(Cuts cuts_)
{
  _cuts = cuts_;
  _cotth1 = 1./tan(cuts_.thetamax*pi/180.);
  _cotth2 = 1./tan(cuts_.thetamin*pi/180.);  
}

void GamGam::SetWRange(double wmin_, double wmax_)
{
  _wmin = wmin_;
  _wmax = wmax_;
  if (setin) {

  }
}

double GamGam::PeriPP(int nup_, int ndown_)
{
  double rho = .585; // 1.05
  double cc1 = .86926; // .6303
  double cc2 = 2.23422; // 2.2049
  double dd1 = .12549; // .0468
  double cp = .96; // 1.23
  double bp = .63; // .61

  double dummy, psfw1, psfw2;
  double u1, u2, v1, v2;
  double en, x, xt;
  double rhot, qqq, qdq;
  double t11, t12, t21, t22;
  double peripp;

#ifdef DEBUG
  std::cout << "[GamGam::PeriPP] [DEBUG]\n  Nup = " << nup_ << "\n  Ndown = " << ndown_ << std::endl;
#endif

  switch(nup_) {
  case 1:
    u1 = 1.;
    u2 = 1.;
    break;
  case 2:
    xt = std::pow(1.-_t1/.71, 2);
    _tau = _t1/(4*_w1);
    u1 = std::pow(2.79/xt, 2);
    u2 = (1./std::pow(xt, 2)-u1*_tau)/(1.-_tau);
    break;
  case 4:
    std::cout << "[GamGam::PeriPP] [DEBUG] Result of PSF : " << PSF(_t1, _w3, &dummy, &psfw1, &psfw2) << std::endl;
    std::cout << "after PSF : " << psfw1 << "\t" << psfw2 << std::endl;
    u1 = -psfw1*(2.*_mp1)/_t1;
    u2 = psfw2/(2.*_mp1);
    break;
  default:
    x = _t1/(_t1-_w3);
    en = _w31-_t1;
    _tau = _t1/(4.*_w1);
    rhot = rho-_t1;
    u1 = (-cc1*std::pow(rho/rhot, 2)*_w31-cc2*_w1*std::pow(1.-x, 4)/(x*(x*cp-2*bp)+1.))/_t1;
    u2 = (-_tau*u1-dd1*_w31*_t1*(rho/rhot)*std::pow(_w31/en, 2)/(rhot*_w1))/(1.-std::pow(en, 2)/(4.*_w1*_t1));
    break;
  }

  switch(ndown_) {
  case 1:
    v1 = 1.;
    v2 = 1.;
    break;
  case 2:
    xt = std::pow(1.-_t2/.71, 2);
    _tau = _t2/(4.*_w2);
    v1 = std::pow(2.79/xt, 2);
    v2 = (1./std::pow(xt, 2)-v1*_tau)/(1.-_tau);
    break;
  default:
    x = _t2/(_t2-_w5);
    en = _w52-_t2;
    _tau = _t2/(4.*_w2);
    rhot = rho-_t2;
    v1 = (-cc1*std::pow(rho/rhot, 2)*_w52-cc2*_w2*std::pow(1.-x, 4)/(x*(x*cp-2.*bp)+1))/_t2;
    v2 = (-_tau*_w1-dd1*_w52*_t2*(rho/rhot)*std::pow(_w52/en, 2)/(rhot*_w2))/(1.-std::pow(en, 2)/(4.*_w2*_t2));
    break;
  }
#ifdef DEBUG
  std::cout << "[GamGam::PeriPP] [DEBUG]"
            << "\n  u1 = " << u1
            << "\n  u2 = " << u2
            << "\n  v1 = " << v1
            << "\n  v2 = " << v2
            << std::endl;
#endif

  qqq = std::pow(_q1dq, 2);
  qdq = 4.*_w6-_w4;
  t11 = 64.*(_bb*(qqq-_g4-qdq*(_t1+_t2+2.*_w6))-2.*(_t1+2.*_w6)*(_t2+2.*_w6)*qqq)*_t1*_t2;
  t12 = 128.*(-_bb*(_d2+_g6)-2.*(_t1+2.*_w6)*(_a2*qqq+std::pow(_a6, 2)))*_t1;
  t21 = 128.*(-_bb*(_d4+_g5)-2.*(_t2+2.*_w6)*(_a1*qqq+std::pow(_a5, 2)))*_t2;
  t22 = 512.*(_bb*(std::pow(_delta, 2)-_gram)-std::pow(_epsi-_delta*(qdq+_q1dq2), 2)-_a1*std::pow(_a6, 2)-_a2*std::pow(_a5, 2)-_a1*_a2*qqq);

  peripp = (((u1*v1*t11+u2*v1*t21+u1*v2*t12+u2*v2*t22)/(_t1*_t2*_bb))/(_t1*_t2*_bb))/4.;
#ifdef DEBUG
  std::cout << "[GamGam::PeriPP] [DEBUG]"
            << "\n  t11 = " << t11
            << "\n  t12 = " << t12
            << "\n  t21 = " << t21
            << "\n  t22 = " << t22
            << "\n  tau = " << _tau
            << "\n  --> PeriPP = " << peripp
            << std::endl;
#endif
  return peripp;
}

Particle* GamGam::GetParticle(int role_)
{
  std::map<int,Particle>::iterator out = this->_part->find(role_);
  if (out!=this->_part->end()) {
    return &(out->second);
  }
  else {
    return new Particle();
  }
}

void GamGam::StoreEvent(std::ofstream *file_, double weight_=1.)
{
  Particle *ol1 = GetParticle(6);
  Particle *ol2 = GetParticle(7);
  *file_ << ol1->e << "\t"
         << ol1->px << "\t"
         << ol1->py << "\t"
         << ol1->pz << "\t"
         << ol1->pt << "\t"
         << ol1->m << "\t"
         << ol1->pdgId << "\t"
         << weight_
         << std::endl;
  *file_ << ol2->e << "\t"
         << ol2->px << "\t"
         << ol2->py << "\t"
         << ol2->pz << "\t"
         << ol2->pt << "\t"
         << ol2->m << "\t"
         << ol2->pdgId << "\t"
         << weight_
         << std::endl;
}

