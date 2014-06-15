#include "gamgamll.h"

GamGamLL::GamGamLL(int nOpt_) :
  _nOpt(nOpt_),
  _ep1(-1), _w1(-1), _ep2(-1), _w2(-1.),
  _w3(-1.), _w4(-1.), _w5(-1.), _w6(-1.), _w7(-1.),
  _p12(0.), _p13(0.), _p14(0.), _p15(0.), _p23(0.), _p24(0.), _p25(0.), _p34(0.), _p35(0.), _p45(0.),
  _p1k2(0.), _p2k1(0.),
  setp1(false), setp2(false), setp3(false), setp5(false), setll(false),
  _cotth1(-99999.), _cotth2(99999.)
{
  _name = PROCESS_NAME;
}

bool
GamGamLL::SetOutgoingParticles(int part_, int pdgId_, int mothRole_)
{
  double mass_, outm, dm;

  if (!_point_set) return false;

  mass_ = GetMassFromPDGId(pdgId_);

  if (mass_<0 or pdgId_==2) { //FIXME!!!
    switch (_cuts.kinematics) {
    case 1: // elastic
    default:
      return false;
    case 2: // single-dissociative
      outm = _mp1;
      mass_ = this->ComputeMX(x(7), outm, &dm);
      break;
    case 3: // double-dissociative
      int ind;
      if (part_==3) { // First outgoing proton remnant
	      outm = _mp1;
	      ind = 7;
      }
      else if (part_==5 && _mp3>0.) { // Second outgoing proton remnant (if first is defined)
	      outm = _mp3;
	      ind = 8;
      }
      else return false;
      mass_ = this->ComputeMX(x(ind), outm, &dm);
    }
  }

  switch(part_) {
  case 3:
    // First outgoing proton (or remnant)
    _mp3 = mass_;
    _w3 = std::pow(_mp3, 2);
    _pdg3 = pdgId_;
    _dw31 = dm;
    setp3 = true;
    break;
  case 5:
    // Second outgoing proton (or remnant)
    _mp5 = mass_;
    _w5 = std::pow(_mp5, 2);
    _pdg5 = pdgId_;
    _dw52 = dm;
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
  _setout = setp3 and setp5 and setll;
  _setkin = _setin and _setout;
#ifdef DEBUG
  std::cout << "[GamGamLL::SetOutgoingParticles] [DEBUG] Particle \"" << part_ << "\" has PDG id " << pdgId_ << std::endl;
  if (_setout) {
    std::cout << "  --> Outgoing state is fully set" << std::endl;
  }
  if (_setkin) {
    std::cout << "  --> Kinematics is fully set" << std::endl;
  }
#endif
  return true;
}

bool
GamGamLL::SetIncomingParticles(Particle ip1_, Particle ip2_)
{
  Particle *p1, *p2;
  double k, *p31, *p32;
  int role1, role2;

  role1 = (ip1_.Pz()>0.) ? 1:2;
  role2 = (ip2_.Pz()>0.) ? 1:2;
  if (role1==role2) return false;
  ip1_.role = role1;
  ip2_.role = role2;

  k = 0.;
  p31 = ip1_.P4();
  p32 = ip2_.P4();
  for (int i=0; i<3; i++) k += p31[i]*p32[i];
  _s = ip1_.M2()+ip2_.M2()+2.*(ip1_.E()*ip2_.E()-k);
  _ecm = sqrt(_s);

  this->_ev->AddParticle(&ip1_);
  this->_ev->AddParticle(&ip2_);

  p1 = this->_ev->GetOneByRole(1);
  p2 = this->_ev->GetOneByRole(2);

  _ep1 = p1->E();
  _mp1 = p1->M();
  _w1 = p1->M2();
  _pp1 = p1->P();
  _pdg1 = p1->pdgId;
  setp1 = true;

  _ep2 = p2->E();
  _mp2 = p2->M();
  _w2 = p2->M2();
  _pp2 = p2->P();
  _pdg2 = p2->pdgId;

  _etot = p1->E()+p2->E();
  _ptot = std::sqrt(std::pow(p1->Px()+p2->Px(), 2)
                   +std::pow(p1->Py()+p2->Py(), 2)
                   +std::pow(p1->Pz()+p2->Pz(), 2));

  _setin = p1->Valid() and p2->Valid();
  _setkin = _setin && _setout;
  return _setkin;
}

bool
GamGamLL::Pickin()
{
  double sig, sig1, sig2; // sig1 = sigma and sig2 = sigma' in [1]
  double sp, ss, st;
  double sb, sd, se;
  double smax;
  double splus, s2x;
  double s2min, s2max;
  double ds2;
  double s1p, s1m, s1pp, s1pm, s2p;
  double sl2, sl3, sl4, sl5, sl6, sl7;
  double dt1, dt2;
  double t13, t25;
  double rl1, rl2, rl4;
  double r1, r2, r3, r4;
  double b, c;
  double ap;
  double yy4;
  double dd, delb;
  double sbb, sdd, see, ssb, ssd, sse;
  double d6, d8;

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] _nOpt = " << _nOpt << std::endl;
#endif
  _dj = 0.;

  _w4 = std::pow(_mc4, 2);

  sig = _mc4+_mp5;
  sig1 = std::pow(sig, 2);
  sig2 = std::pow(sig, 2);

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] mc4 = " << _mc4 << std::endl;
  std::cout << "[GamGamLL::Pickin] [DEBUG] sig1 = " << sig1 << std::endl;
  std::cout << "[GamGamLL::Pickin] [DEBUG] sig2 = " << sig2 << std::endl;
#endif

  // Mass difference between the first outgoing particle and the first incoming
  // particle
  _w31 = _w3-_w1;
  // Mass difference between the second outgoing particle and the second
  // incoming particle
  _w52 = _w5-_w2;
  // Mass difference between the two incoming particles
  _w12 = _w1-_w2;
  // Mass difference between the central two-photons system and the second
  // outgoing particle
  d6 = _w4-_w5;

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG]" << std::endl
            << "  w1 = " << _w1 << std::endl
            << "  w2 = " << _w2 << std::endl
            << "  w3 = " << _w3 << std::endl
            << "  w4 = " << _w4 << std::endl
            << "  w5 = " << _w5 << std::endl;
#endif

  ss = _s+_w12;
  rl1 = std::pow(ss, 2)-4.*_w1*_s; // lambda(s, m1**2, m2**2)

  if (rl1<=0.) {
    return false;
  }
  _sl1 = std::sqrt(rl1);

  _s2 = ds2 = 0.;
  if (_nOpt==0) {
    smax = _s+_w3-2.*_mp3*_ecm;
    Map(x(2), sig1, smax, &_s2, &ds2);
    sig1 = _s2; //FIXME!!!!!!!!!!!!!!!!!!!!
  }
#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] _s2 = " << _s2 << std::endl;
#endif

  //std::cout << "s=" << _s << ", w3=" << _w3 << ", sig1=" << sig1 << std::endl;
  sp = _s+_w3-sig1;

  _d3 = sig1-_w2;

  rl2 = std::pow(sp, 2)-4.*_s*_w3; // lambda(s, m3**2, sigma)
  if (rl2<=0.) {
    return false;
  }
  sl2 = std::sqrt(rl2);

  //std::cout << "ss=" << ss << ", sp=" << sp << ", sl1=" << _sl1 << ", sl2=" << sl2 << std::endl;
  _t1max = _w1+_w3-(ss*sp+_sl1*sl2)/(2.*_s); // definition from eq. (A.4) in [1]
  _t1min = (_w31*_d3+(_d3-_w31)*(_d3*_w1-_w31*_w2)/_s)/_t1max; // definition from eq. (A.5) in [1]

  // FIXME dropped in CDF version
  if (_t1max>-_cuts.q2min or (_cuts.q2max!=-1. and _t1min<-_cuts.q2max)) {
    /*std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
              << "     t1max = " << std::setw(8) << _t1max << " > -q2min = " << std::setw(8) << -_cuts.q2min << std::endl
              << "  or t1min = " << std::setw(8) << _t1min << " < -q2max = " << std::setw(8) << -_cuts.q2max << std::endl;*/
    return false;
  }
  if (_cuts.q2max!=-1. and _t1max<-_cuts.q2max) {
    _t1max = -_cuts.q2max;
  }
  if (_t1min>-_cuts.q2min) {
    _t1min = -_cuts.q2min;
  }
  /////

  // t1, the first photon propagator, is defined here
  Map(x(0), _t1min, _t1max, &_t1, &dt1);
  // changes wrt mapt1 : dx->-dx
  dt1 = -dt1;
#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] definition of t1 according to"
            << " (t1min, t1max) = (" << _t1min << ", " << _t1max << ")" << std::endl
            << "  _t1 = " << _t1 << std::endl;
#endif

  _dd4 = _w4-_t1;
  d8 = _t1-_w2;

  t13 = _t1-_w1-_w3;

  _sa1 =-std::pow(_t1-_w31, 2)/4.+_w1*_t1;
  if (_sa1>=0.) {
    std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
              << "  _sa1>=0 : " << _sa1 << std::endl;
    return false;
  }

  sl3 = std::sqrt(-_sa1);

  // one computes splus and (s2x=s2max)
  if (_w1!=0.) {
    sb =(_s*(_t1-_w31)+_w12*t13)/(2.*_w1)+_w3;
    sd = _sl1*sl3/_w1;
    se =(_s*(_t1*(_s+t13-_w2)-_w2*_w31)+_w3*(_w12*d8+_w2*_w3))/_w1;
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
    std::cout << 3 << std::endl;
    // 3
    s2max = (_s*(_t1*(_s+d8-_w3)-_w2*_w3)+_w2*_w3*(_w2+_w3-_t1))/(ss*t13);
    splus = sig2;
  }
  // 4
  s2x = s2max;
#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] s2x = s2max = " << s2x << std::endl;
#endif

  if (_nOpt<0) { // 5
    if (splus>sig2) {
      sig2 = splus;
#ifdef DEBUG
      std::cout << "[GamGamLL::Pickin] [DEBUG] sig2 truncated to splus = " << splus << std::endl;
#endif
    }
    if (_nOpt<-1) {
      Map(x(2), sig2, s2max, &_s2, &ds2);
    }
    else { // _nOpt==-1
      Mapla(_t1, _w2, x(2), sig2, s2max, &_s2, &ds2);
    }
    s2x = _s2;
  }
  else if (_nOpt==0) { // 6
    s2x = _s2;
  }

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] s2x = " << s2x << std::endl;
#endif
  // 7
  r1 = s2x-d8;
  r2 = s2x-d6;
  rl4 = (std::pow(r1, 2)-4.*_w2*s2x)*(std::pow(r2, 2)-4.*_w5*s2x);
  if (rl4<=0.) {
    /*std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
      << "  rl4<=0 : " << rl4 << std::endl;*/
    return false;
  }
  sl4 = std::sqrt(rl4);
  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  _t2max = _w2+_w5-(r1*r2+sl4)/(2.*s2x);
  _t2min = (_w52*_dd4+(_dd4-_w52)*(_dd4*_w2-_w52*_t1)/s2x)/_t2max;

  // t2, the second photon propagator, is defined here
  Map(x(1), _t2min, _t2max, &_t2, &dt2);
  // changes wrt mapt2 : dx->-dx

  dt2 = -dt2;

  _tau = _t1-_t2;
  r3 = _dd4-_t2;
  r4 = _w52-_t2;

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG]" << std::endl
            << "  r1 = " << r1 << std::endl
            << "  r2 = " << r2 << std::endl
            << "  r3 = " << r3 << std::endl
            << "  r4 = " << r4 << std::endl;
#endif

  b = r3*r4-2.*(_t1+_w2)*_t2;
  c = _t2*d6*d8+(d6-d8)*(d6*_w2-d8*_w5);
  t25 = _t2-_w2-_w5;

  _sa2 = -std::pow(r4, 2)/4.+_w2*_t2;
  if (_sa2>=0.) {
    std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
	      << "  _sa2 = " << _sa2 << " >= 0" << std::endl;
    return false;
  }
  sl6 = 2.*std::sqrt(-_sa2);
  _g4 = -std::pow(r3, 2)/4.+_t1*_t2;
  if (_g4>=0.) {
    std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
	      << "  _g4 = " << _g4 << " >= 0" << std::endl;
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
    Map(x(2), s2min, s2max, &_s2, &ds2);
  }
  else if (_nOpt==1) { // _nOpt=1
    Mapla(_t1, _w2, x(2), s2min, s2max, &_s2, &ds2);
  }
  ap = -std::pow(_s2+d8, 2)/4.+_s2*_t1;
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
    std::cerr << "[GamGamLL::Pickin] [ERROR] : dd2 == NaN" << std::endl;
    std::cerr << "  dd2 = " << _dd2 << std::endl
              << "  s2 = " << _s2 << std::endl
              << "  s2p = " << s2p << std::endl
              << "  s2min = " << s2min << std::endl
              << "  t2min = " << _t2min << std::endl
              << "  t2max = " << _t2max << std::endl;
#endif
  }
  /////
  if (x(3)>1. or x(3)<-1.) {
    std::cerr << "[GamGamLL::Pickin] [FATAL] x[3] = " << x(3) << std::endl;
    return false;
  }
  yy4 = cos(pi*x(3));
  dd = _dd1*_dd2;
  _p12 = (_s-_w1-_w2)/2.;
  st = _s2-_t1-_w2;
  delb = (2.*_w2*r3+r4*st)*(4.*_p12*_t1-(_t1-_w31)*st)/(16.*ap);

  if (dd<=0.) {
    /*std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
      << "  dd = " << dd << " <= 0" << std::endl;*/
    return false;
  }

  _delta = delb-yy4*st*std::sqrt(dd)/(2.*ap);
  _s1 = _t2+_w1+(2.*_p12*r3-4.*_delta)/st;

  if (ap>=0.) {
    std::cerr << "[GamGamLL::Pickin] [FATAL]" << std::endl
	      << "  ap = " << ap << " >= 0" << std::endl;
    return false;
  }

  _dj = ds2*dt1*dt2*std::pow(pi, 2)/(8.*_sl1*std::sqrt(-ap));

#ifdef DEBUG
  std::cout << "[GamGamLL::Pickin] [DEBUG] _dj = " << _dj << std::endl;
#endif

  _gram = (1.-std::pow(yy4, 2))*dd/ap;

  _p13 = -t13/2.;
  _p14 = (_tau+_s1-_w3)/2.;
  _p15 = (_s+_t2-_s1-_w2)/2.;
  _p23 = (_s+_t1-_s2-_w1)/2.;
  _p24 = (_s2-_tau-_w5)/2.;
  _p25 = -t25/2.;
  _p34 = (_s1-_w3-_w4)/2.;
  _p35 = (_s+_w4-_s1-_s2)/2.;
  _p45 = (_s2-_w4-_w5)/2.;

  _p1k2 = (_s1-_t2-_w1)/2.;
  _p2k1 = st/2.;


  if (_w2!=0.) {
    sbb = (_s*(_t2-_w52)-_w12*t25)/(2.*_w2)+_w5;
    sdd = _sl1*sl6/(2.*_w2);
    see = (_s*(_t2*(_s+t25-_w1)-_w1*_w52)+_w5*(_w1*_w5-_w12*(_t2-_w1)))/_w2;
    if (sbb/sdd>=0.) {
      s1p = sbb+sdd;
      s1m = see/s1p;
      // FIXME there should be a more beautiful way to check for nan!
      // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
      // FIXME dropped in CDF version
      if (!(_dd2>=0.) && !(_dd2<0.)) { // NaN
#ifdef ERROR
        std::cout << "[GamGamLL::Pickin] [ERROR] : dd2 == NaN" << std::endl;
        std::cout << "  dd2 = " << _dd2 << std::endl
		  << "   s1 = " << _s1 << std::endl
                  << "  s1p = " << s1p << std::endl
		  << "  s1m = " << s1m << std::endl
                  << "   w2 = " << _w2 << std::endl;
#endif
      }
      /////
    }
    else { // 12
      s1m = sbb-sdd;
      s1p = see/s1m;
    }
    _dd3 = -_w2*(s1p-_s1)*(s1m-_s1)/4.; // 13
  }
  else { // 14
    s1p = (_s*(_t2*(_s-_w5+_t2-_w1)-_w1*_w5)+_w1*_w5*(_w1+_w5-_t2))/(t25*(_s-_w12));
    _dd3 = -t25*(_s-_w12)*(s1p-_s1)/4.;
  }
  // 15
  _acc3 = (s1p-_s1)/(s1p+_s1);

  ssb = _t2+_w1-r3*(_w31-_t1)/(2.*_t1);
  ssd = sl3*sl7/_t1;
  sse = (_t2-_w1)*(_w4-_w3)+(_t2-_w4+_w31)*((_t2-_w1)*_w3-(_w4-_w3)*_w1)/_t1;

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
  _dd5 = _dd1+_dd3+((_p12*(_t1-_w31)/2.-_w1*_p2k1)*(_p2k1*(_t2-_w52)-_w2*r3)-_delta*(2.*_p12*_p2k1-_w2*(_t1-_w31)))/_p2k1;

  return true;
}

bool
GamGamLL::Orient()
{
  double re, rr, a1;

  bool pck = Pickin();
  if (!pck || _dj==0.) {
    //std::cerr << "[GamGamLL::Orient] [FATAL] Pickin() failed : (bool)" << pck << ", _dj = " << _dj << std::endl;
    return false;
  }
  re = 1./(2.*_ecm);
  _ep1 = re*(_s+_w12);
  _ep2 = re*(_s-_w12);

#ifdef DEBUG
  std::cout << "[GamGamLL::Orient] [DEBUG]" << std::endl
	    << "  re = " << re << std::endl
	    << "  _w12 = " << _w12 << std::endl;
  std::cout << "[GamGamLL::Orient] [DEBUG] incoming particles' energy = " << _ep1 << ", " << _ep2 << std::endl;
#endif

  _p = re*_sl1;

  _de3 = re*(_s2-_w3+_w12);
  _de5 = re*(_s1-_w5-_w12);

  // Final state energies
  _ep3 = _ep1-_de3;
  _ec4 = _de3+_de5;
  _ep5 = _ep2-_de5;

  /*if (_ep3>_ep1 or _ep5>_ep2) {
    // FIXME workaround to avoid very unphysical events...
    return false;
  }*/

  if (_ec4<_mc4) {
    /*std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
      << "  _ec4<_mc4 : _ec4 = " << _ec4 << ", _mc4 = " << _mc4 << std::endl;*/
    return false;
  }
  // What if the protons' momenta are not along the z-axis?
  _pp3 = std::sqrt(std::pow(_ep3, 2)-_w3);
  _pc4 = std::sqrt((std::pow(_ec4, 2)-std::pow(_mc4, 2)));

  if (_pc4==0.) {
    /*std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
      << "  _pzc4==0" << std::endl;*/
    return false;
  }
  _pp5 = std::sqrt(std::pow(_ep5, 2)-_w5);
  _p_p3 = std::sqrt(_dd1/_s)/_p;

#ifdef DEBUG
  std::cout << "[GamGamLL::Orient] [DEBUG] central system's energy : E4 = " << _ec4 << std::endl;
  std::cout << "[GamGamLL::Orient] [DEBUG] central system's momentum : P4 = " << _pc4 << std::endl;
  std::cout << "[GamGamLL::Orient] [DEBUG] central system's invariant mass : M4 = " << _mc4 << std::endl;
  std::cout << "[GamGamLL::Orient] [DEBUG] outgoing particles' energy : E3 = " << _ep3 << ", E5 = " << _ep5 << std::endl;
#endif

  _p_p5 = std::sqrt(_dd3/_s)/_p;
  _st3 = _p_p3/_pp3;
  _st5 = _p_p5/_pp5;

#ifdef DEBUG
  std::cout << "[GamGamLL::Orient] [DEBUG] _st3 = " << _st3 << ", _st5 = " << _st5 << std::endl;
#endif

  // FIXME there should be a more beautiful way to check for nan!
  // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
  // FIXME dropped in CDF version
  if (!(_dd3>=0.) && !(_dd3<0.)) { // NaN
#ifdef ERROR
    std::cerr << "[GamGamLL::Orient] [ERROR]" << std::endl
	      << "  dd3 == NaN" << std::endl;
#endif
  }
  if (!(_dd1>=0.) && !(_dd1<0.)) { // NaN
#ifdef ERROR
    std::cerr << "[GamGamLL::Orient] [ERROR]" << std::endl
	      << "  dd1 == NaN" << std::endl;
#endif
  }
  /////

  if (_st3>1. or _st5>1.) {
    std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
	      << "  _st3>1 or _st5>1 : _st3 = " << _st5 << ", _st5 = " << _st5 << std::endl;
    return false;
  }
  _ct3 = std::sqrt(1.-std::pow(_st3, 2));
  _ct5 = std::sqrt(1.-std::pow(_st5, 2));

  if (_ep1*_ep3<_p13) {
    _ct3 = -_ct3;
  }

#ifdef DEBUG
  std::cout << "[GamGamLL::Orient] [DEBUG] _ct3 = " << _ct3 << ", _ct5 = " << _ct5 << std::endl;
#endif

  if (_ep2*_ep5>_p25) {
    _ct5 = -_ct5;
  }
  _al3 = std::pow(_st3, 2)/(1.+_ct3);
  _be5 = std::pow(_st5, 2)/(1.-_ct5);

  if (_dd5<0.) {
    /*std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
      << "  _dd5<0 : " << _dd5 << std::endl;*/
    return false;
  }

  // Centre of mass system kinematics (theta4 and phi4)
  _p_p4 = std::sqrt(_dd5/_s)/_p;
  _st4 = _p_p4/_pc4;

  if (_st4>1.) {
    std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
	      << "  _st4>1 : " << _st4 << std::endl;
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
  std::cout << "[GamGamLL::Orient] [DEBUG] _ct4 = " << _ct4 << ", _al4 = " << _al4 << ", _be4 = " << _be4 << std::endl;
#endif

  //std::cout << "pp3 = " << _p_p3 << ", pp5 = " << _p_p5 << std::endl;

  rr  = std::sqrt(-_gram/_s)/(_p*_p_p4);
  _sp3 = rr/_p_p3;
  _sp5 = -rr/_p_p5;
  //std::cout << "rr = " << rr << ", sp3 = " << fabs(_sp3) << ", sp5 = " << fabs(_sp5) << std::endl;

  if (fabs(_sp3)>1. || fabs(_sp5)>1.) {
    /*std::cerr << "[GamGamLL::Orient] [FATAL]" << std::endl
	      << "  |_sp3|>1 or |_sp5|>1 :" << std::endl
              << "   |_sp3| = " << fabs(_sp3) << std::endl
              << "  _sp3**2 = " << std::pow(_sp3, 2) << std::endl
              << "   |_sp5| = " << fabs(_sp5) << std::endl
              << "  _sp5**2 = " << std::pow(_sp5, 2) << std::endl;*/
    return false;
  }

  _cp3 = -std::sqrt(1.-std::pow(_sp3, 2));
  _cp5 = -std::sqrt(1.-std::pow(_sp5, 2));

  a1 = _p_p3*_cp3-_p_p5*_cp5; //OK!!!

#ifdef DEBUG
  std::cout << "[GamGamLL::Orient] [DEBUG] Kinematic quantities" << std::endl
            << "  cos(theta3) = " << _ct3 << "\t  sin(theta3) = " << _st3 << std::endl
            << "  cos( phi3 ) = " << _cp3 << "\t  sin( phi3 ) = " << _sp3 << std::endl
            << "  cos(theta4) = " << _ct4 << "\t  sin(theta4) = " << _st4 << std::endl
            << "  cos( phi4 ) = " << _ct4 << "\t  sin( phi4 ) = " << _st4 << std::endl
            << "  cos(theta5) = " << _ct5 << "\t  sin(theta5) = " << _ct5 << std::endl
            << "  cos( phi5 ) = " << _cp5 << "\t  sin( phi5 ) = " << _cp5 << std::endl;
#endif

  if (fabs(_p_p4+_p_p3*_cp3+_cp5*_p_p5)<fabs(fabs(a1)-_p_p4)) {
#ifdef DEBUG
    std::cout << "[GamGamLL::Orient] [DEBUG] fabs(_p_p4+_p_p3*_cp3+_cp5*_p_p5)<fabs(fabs(a1)-_p_p4)" << std::endl
              << "  pp4 = " << _p_p4 << std::endl
              << "  pp5 = " << _p_p5 << std::endl
              << "  cos(phi3) = cp3 = " << _cp3 << std::endl
              << "  cos(phi5) = cp5 = " << _cp5 << std::endl
              << "  a1 = " << a1 << std::endl;
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

double
GamGamLL::ComputeMX(double x_, double outmass_, double *dw_)
{
  double wx2min, wx2max;
  double mx2, dmx2;

  /*wx2min = std::pow(std::max(GetMassFromPDGId(2212)+GetMassFromPDGId(211), _cuts.mxmin), 2);
    wx2max = std::pow(std::min(_ecm-_mp2-2.*outmass_, _cuts.mxmax), 2);*/
  
  wx2min = std::pow(GetMassFromPDGId(2212)+GetMassFromPDGId(211), 2);
  wx2max = std::pow(_ecm-_mp2-2.*outmass_, 2);
  Map(x_, wx2min, wx2max, &mx2, &dmx2);

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeMX] [DEBUG]" << std::endl
	    << "\tMX**2 in range [" << wx2min << ", " << wx2max << "]" << std::endl
	    << "\tx = " << x_ << std::endl
	    << "\tMX**2 = " << mx2 << ", dMX**2 = " << dmx2 << std::endl
	    << "\tMX = " << sqrt(mx2) << ", dMX = " << sqrt(dmx2) << std::endl;
#endif

  *dw_ = sqrt(dmx2);
  return sqrt(mx2);
}

double
//GamGamLL::ComputeWeight(int nm_)
GamGamLL::ComputeWeight()
{
  int nm_ = 1; //FIXME...

  // WARNING ====> PP5-->_p_p5
  //                P5-->_pp5
  double weight;
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

  // COMMON /QVEC/   // 0 = E, 1-3 = p
  double qve[4];

  weight = 0.;

  if (!_setout) {
    std::cerr << "[GamGamLL::ComputeWeight] [FATAL]" << std::endl
	      << "  Output state not set !" << std::endl;
    return 0.;
  }

  if (_cuts.wmax<0) _cuts.wmax = _s;

  // The minimal energy for the central system is its outgoing leptons' mass energy (or wmin_ if specified)
  wmin = std::pow(_ml6+_ml7,2);
  if (fabs(wmin)<fabs(_cuts.wmin)) wmin = _cuts.wmin;

  // The maximal energy for the central system is its CM energy with the outgoing particles' mass energy substracted (or _wmax if specified)
  wmax = std::pow(_ecm-_mp3-_mp5,2);
  if (fabs(wmax)>fabs(_cuts.wmax)) wmax = _cuts.wmax;

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG]" << std::endl
            << "  wmin = " << wmin << std::endl
            << "  wmax = " << wmax << std::endl
            << "  wmax/wmin = " << wmax/wmin << std::endl;
#endif
  Map(x(4),wmin,wmax,&_w4,&dw4);
  _mc4 = std::sqrt(_w4);

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG] Computed value for w4 = " << _w4 << " -> mc4 = " << _mc4 << std::endl;
#endif

  if (!Orient()) {
    return 0.;
  }

  if (_t1>0. or _t2>0.) {
    _dj = 0.;
  }
  if (_dj==0.) {
    /*std::cerr << "[GamGamLL::ComputeWeight] [FATAL]" << std::endl
      << "  _dj = " << _dj << std::endl;*/
    return 0.;
  }
  ecm6 = (_w4+_w6-_w7)/(2.*_mc4);
  pcm6 = std::sqrt(std::pow(ecm6, 2)-_w6);
  
  _dj *= dw4*pcm6/(_mc4*sconstb*_s);
  
  // Let the most obscure part of this code begin...

  e1mp1 = _w1/(_ep1+_p);
  e3mp3 = _w3/(_ep3+_pp3);

  // 2-photon system kinematics ?!
  eg = (_w4+_t1-_t2)/(2.*_mc4);
  pg = std::sqrt(std::pow(eg, 2)-_t1);
  pgx = -_p_p3*_cp3*_ct4-_st4*(_de3-e1mp1+e3mp3+_pp3*_al3);
  pgy = -_p_p3*_sp3;
  pgz = _mc4*_de3/(_ec4+_pc4)-_ec4*_de3*_al4/_mc4-_p_p3*_cp3*_ec4*_st4/_mc4+_ec4*_ct4/_mc4*(_pp3*_al3+e3mp3-e1mp1);

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG] pg3 = ("
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
  
  // Phi angle for the 2-photon system ?!
  cpg = pgx/pgp;
  spg = pgy/pgp;
  
  // Theta angle for the 2-photon system ?!
  stg = pgp/pg;
  ctg = std::sqrt(1.-std::pow(stg, 2));
  if (pgz<0.) {
    ctg = -ctg;
  }
  
  xx6 = x(5);
  
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
    ddd = (amap+bmap*_ctcm6)*(amap-bmap*_ctcm6)/amap/bmap*log(ymap);
    _dj *= ddd/2.;
  }

  // 3D rotation of the first outgoing lepton wrt the CM system
  _ctcm6 = 1.-2.*xx6; // cos(theta_cm,6) is between -1 and 1
  _stcm6 = 2.*std::sqrt(xx6*(1.-xx6)); // definition is OK (according to _ctcm6 def)
#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG]" << std::endl
	    << "  ctcm6 = " << _ctcm6 << std::endl
	    << "  stcm6 = " << _stcm6 << std::endl;
#endif

  phicm6 = 2.*pi*x(6);

  cpcm6 = cos(phicm6);
  spcm6 = sin(phicm6);

  // First outgoing lepton's 3-momentum in the centre of mass system
  pcm6x = pcm6*_stcm6*cpcm6;
  pcm6y = pcm6*_stcm6*spcm6;
  pcm6z = pcm6*_ctcm6;

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG] p3cm6 = ("
            << pcm6x << ", "
            << pcm6y << ", "
            << pcm6z << ")"
            << std::endl;
#endif

  pc6z = ctg*pcm6z-stg*pcm6x;

  h1 = stg*pcm6z+ctg*pcm6x;
  
  pc6x = cpg*h1-spg*pcm6y;
  
  qcx = 2.*pc6x;
  qcz = 2.*pc6z;
  // qcy == QCY is never defined
  
  // First outgoing lepton's 3-momentum
  p6y = cpg*pcm6y+spg*h1;
  _el6 = (_ec4*ecm6+_pc4*pc6z)/_mc4;
  h2 = (_ec4*pc6z+_pc4*ecm6)/_mc4;
  p6x = _ct4*pc6x+_st4*h2;
  p6z = _ct4*h2-_st4*pc6x;
  
  qve[0] = _pc4*qcz/_mc4;
  qve[2] = 2.*p6y;
  hq = _ec4*qcz/_mc4;
  qve[1] = _ct4*qcx+_st4*hq;
  qve[3] = _ct4*hq-_st4*qcx;

  _pl6 = std::sqrt(std::pow(_el6, 2)-_w6); // first outgoing lepton's |p|

  // Available energy for the second lepton is the 2-photon system's energy with the first lepton's energy removed
  _el7 = _ec4-_el6;
  _pl7 = std::sqrt(std::pow(_el7, 2)-_w7); // second outgoing lepton's |p|

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG] (outgoing kinematics)" << std::endl
            << "   first outgoing lepton : p, E = " << _pl6 << ", " << _el6 << std::endl
            << "  second outgoing lepton : p, E = " << _pl7 << ", " << _el7 << std::endl;
#endif

  double p7x, p7y, p7z;
  // Second outgoing lepton's 3-momentum
  p7x = _p_p4-p6x;
  p7y = -p6y;
  p7z = _pc4*_ct4-p6z;
  
  pp6 = std::sqrt(std::pow(p6x, 2)+std::pow(p6y, 2)); // pT of first outgoing lepton
  pp7 = std::sqrt(std::pow(p7x, 2)+std::pow(p7y, 2)); // pT of second outgoing lepton

  // First outgoing lepton's kinematics (sin/cos theta/phi)
  _ct6 = p6z/_pl6;
  _st6 = pp6/_pl6;
  _cp6 = p6x/pp6;
  _sp6 = p6y/pp6;

  // Second outgoing lepton's kinematics (sin/cos theta/phi)
  _ct7 = p7z/_pl7;
  _st7 = pp7/_pl7;
  _cp7 = p7x/pp7;
  _sp7 = p7y/pp7;

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG] (outgoing trajectories)" << std::endl
            << "   first outgoing lepton : cos(theta) = " << _ct6 << ", sin(theta) = " << _st6 << std::endl
            << "   first outgoing lepton : cos( phi ) = " << _cp6 << ", sin( phi ) = " << _sp6 << std::endl
            << "  second outgoing lepton : cos(theta) = " << _ct7 << ", sin(theta) = " << _st7 << std::endl
            << "  second outgoing lepton : cos( phi ) = " << _cp7 << ", sin( phi ) = " << _sp7 << std::endl;
#endif

  _q1dq = eg*(2.*ecm6-_mc4)-2.*pg*pcm6*_ctcm6;
  _q1dq2 = (_w4-_t1-_t2)/2.;

  _bb = _t1*_t2+(_w4*std::pow(_stcm6, 2)+4.*_w6*std::pow(_ctcm6, 2))*std::pow(pg, 2);
  // Q0=QVE[0], QX=QVE[1], QY=QVE[2], QZ=QVE[3]
  c1 = (qve[1]*_sp3-qve[2]*_cp3)*_p_p3;
  c2 = (qve[3]*_ep1-qve[0]*_p)*_p_p3;
  c3 = (_w31*std::pow(_ep1, 2)+2.*_w1*_de3*_ep1-_w1*std::pow(_de3, 2)+std::pow(_p_p3, 2)*std::pow(_ep1, 2))/(_ep3*_p+_pp3*_ct3*_ep1);
  
  b1 = (qve[1]*_sp5-qve[2]*_cp5)*_p_p5; //OK
  b2 = (qve[3]*_ep2+qve[0]*_p)*_p_p5; //OK
  b3 = (_w52*std::pow(_ep2, 2)+2.*_w2*_de5*_ep2-_w2*std::pow(_de5, 2)+std::pow(_p_p5*_ep2, 2))/(_ep2*_pp5*_ct5-_ep5*_p); //OK
  
  r12 = c2*_sp3+qve[2]*c3;
  r13 = -c2*_cp3-qve[1]*c3;

#ifdef DEBUG
  std::cout << "[GamGamLL::ComputeWeight] [DEBUG]" << std::endl;
  for (unsigned int i=0; i<sizeof(qve)/sizeof(double); i++) {
    std::cout << "  qve[" << i << "] = " << qve[i] << std::endl;
  }
#endif

  r22 = b2*_sp5+qve[2]*b3;
  r23 = -b2*_cp5-qve[1]*b3;
  
  _epsi = _p12*c1*b1+r12*r22+r13*r23;

  _g5 = _w1*std::pow(c1, 2)+std::pow(r12, 2)+std::pow(r13, 2);
  _g6 = _w2*std::pow(b1, 2)+std::pow(r22, 2)+std::pow(r23, 2);

  _a5 = -(qve[1]*_cp3+qve[2]*_sp3)*_p_p3*_p1k2-(_ep1*qve[0]-_p*qve[3])*(_cp3*_cp5+_sp3*_sp5)*_p_p3*_p_p5+(_de5*qve[3]+qve[0]*(_p+_pp5*_ct5))*c3; //OK
  _a6 = -(qve[1]*_cp5+qve[2]*_sp5)*_p_p5*_p2k1-(_ep2*qve[0]+_p*qve[3])*(_cp3*_cp5+_sp3*_sp5)*_p_p3*_p_p5+(_de3*qve[3]-qve[0]*(_p-_pp3*_ct3))*b3; //OK

  ////////////////////////////////////////////////////////////////
  // END of GAMGAMLL subroutine in the FORTRAN version
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // INFO from f.f
  ////////////////////////////////////////////////////////////////

  _gamma = _etot/_ecm;
  _betgam = _ptot/_ecm;

  if (_cuts.mode==0) {
#ifdef DEBUG
    std::cout << "[GamGamLL::ComputeWeight] [DEBUG] No cuts applied on the outgoing leptons kinematics !" << std::endl;
#endif
  }
  // Kinematics computation for both leptons
  _pt_l6 = _pl6*_st6;
  pz6 = _betgam*_el6+_gamma*_pl6*_ct6;
  _e6lab = _gamma*_el6+_betgam*_pl6*_ct6;

  _pt_l7 = _pl7*_st7;
  pz7 = _betgam*_el7+_gamma*_pl7*_ct7;
  _e7lab = _gamma*_el7+_betgam*_pl7*_ct7;

  lcut = false; // Event discarded by default
  cott6 = pz6/_pt_l6;
  cott7 = pz7/_pt_l7;

  // Cuts on outgoing leptons' kinematics

  lmu1 = cott6>=_cotth1
     and cott6<=_cotth2
     and (_pt_l6>=_cuts.ptmin or _cuts.ptmin<=0.)
     and (_pt_l6<=_cuts.ptmax or _cuts.ptmax<=0.)
     and (_e6lab>=_cuts.emin  or _cuts.emin <=0.)
     and (_e6lab<=_cuts.emax  or _cuts.emax <=0.);
  lmu2 = cott7>=_cotth1
     and cott7<=_cotth2
     and (_pt_l7>=_cuts.ptmin or _cuts.ptmin<=0.)
     and (_pt_l7<=_cuts.ptmax or _cuts.ptmax<=0.)
     and (_e7lab>=_cuts.emin  or _cuts.emin <=0.)
     and (_e7lab<=_cuts.emax  or _cuts.emax <=0.);

  switch (_cuts.mode) {
    case 0:
    default:
      lcut = true;
      break;
    case 1: // Vermaseren's hypothetical detector cuts
      cost6 = pz6/std::sqrt(std::pow(pz6,2)+std::pow(_pt_l6,2));
      cost7 = pz7/std::sqrt(std::pow(pz7,2)+std::pow(_pt_l7,2));
      lcut = ((fabs(cost6)<=0.75 and _pt_l6>=1.) or (fabs(cost6)<=0.95 and fabs(cost6)>0.75 and fabs(pz6)>1.)) and
             ((fabs(cost7)<=0.75 and _pt_l7>=1.) or (fabs(cost7)<=0.95 and fabs(cost7)>0.75 and fabs(pz7)>1.));
      break;
    case 2:
      lcut = lmu1 and lmu2;
      break;
    case 3:
      lcut = lmu1 or lmu2;
      break;
  }

  // Cut on mass of final hadronic system (MX)
  if (_cuts.kinematics>1) {
    if (_mp3<_cuts.mxmin or _mp3>_cuts.mxmax) return 0.;
    if (_cuts.kinematics==3) {
      if (_mp5<_cuts.mxmin or _mp5>_cuts.mxmax) return 0.;
    }
  }

  // Cut on the proton's Q2 (first photon propagator T1)
  if ((_cuts.q2max!=-1. and _t1<-_cuts.q2max) or _t1>-_cuts.q2min) {
    lcut = false;
  }

  if (!lcut) { // Dismiss the cuts-failing events in the cross-section computation
    return 0.;
  }

  int intgp, intge;

  switch (_cuts.kinematics) { // FIXME inherited from CDF version
  default:
  case 0: // ep case
    intgp = intge = 1; // DESY
    //intgp = intge = 0; // CDF
    weight = sconst*_dj*this->PeriPP(intgp, intge);
    break;
  case 1: // elastic case
    intgp = intge = 2; // DESY
    //intgp = intge = 1; // CDF
    weight = sconst*_dj*this->PeriPP(intgp, intge);
    break;
  case 2: // single-dissociative case
    intgp = 3; // DESY
    intge = 2; // DESY
    /*intgp = 2; // CDF
      intge = 1; // CDF*/
    weight = sconst*_dj*this->PeriPP(intgp, intge)*std::pow(_dw31,2);
    break;
  case 3: // double-dissociative case
    intgp = intge = 3; // DESY
    //intgp = intge = 2; // CDF
    weight = sconst*_dj*this->PeriPP(intgp, intge)*std::pow(_dw31*_dw52,2);
    break;
  }

  return weight;
}

void
GamGamLL::FillKinematics(bool symmetrise_)
{
#ifdef DEBUG
  double gmux, gmuy, gmuw, gmunu;
#endif
  double ranphi, cp, sp;
  int rany, ransign, ranz, role;
  double plab_ip1[4], plab_ip2[4], plab_op1[4], plab_op2[4];
  double plab_ol1[4], plab_ol2[4], plab_ph1[4], plab_ph2[4];
  
  // Needed to parametrise a random rotation around z-axis
  rany = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
  ransign = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
  ranphi = ((double)rand()/RAND_MAX)*2.*pi;
  ranz = 1;
  if (symmetrise_) {
    ranz = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
    //_pp3 *= ranz;
    //_pp5 *= ranz;
  }
  cp = cos(ranphi);
  sp = sin(ranphi);
  
  // First incoming proton
  Particle ip1(1, _pdg1);
  plab_ip1[0] = 0.;
  plab_ip1[1] = 0.;
  plab_ip1[2] = _gamma*_p  +_betgam*_ep1;
  plab_ip1[3] = _gamma*_ep1+_betgam*_p;
  if (!ip1.P(0., 0., plab_ip1[2], plab_ip1[3])) {
    std::cerr << "Invalid incoming proton 1" << std::endl;
  }
  this->_ev->AddParticle(&ip1, true);
  
  // Second incoming proton
  Particle ip2(2, _pdg2);
  plab_ip2[0] = 0.;
  plab_ip2[1] = 0.;
  plab_ip2[2] = -_gamma*_p  +_betgam*_ep2;
  plab_ip2[3] =  _gamma*_ep2-_betgam*_p;
  if (!ip2.P(0., 0., plab_ip2[2], plab_ip2[3])) {
    std::cerr << "Invalid incoming proton 2" << std::endl;
  }
  this->_ev->AddParticle(&ip2, true);
  
  // First outgoing proton
  Particle op1(3, _pdg3);
  plab_op1[0] = _pp3*_st3*_cp3;
  plab_op1[1] = _pp3*_st3*_sp3;
  plab_op1[2] = _gamma*_pp3*_ct3+_betgam*_ep3;
  plab_op1[3] = _gamma*_ep3     +_betgam*_pp3*_ct3;
  if (!op1.P( plab_op1[0]*cp+rany*plab_op1[1]*sp,
             -plab_op1[0]*sp+rany*plab_op1[1]*cp,
              plab_op1[2],
              plab_op1[3])) {
    std::cerr << "Invalid outgoing proton 1" << std::endl;
  }
  if (_cuts.kinematics>1) {
    op1.status = -2;
    //std::cout << "M before : " << op1.M() << std::endl;
    op1.M(_mp3);
    //std::cout << "M as mp3 : " << _mp3 << std::endl;
    //std::cout << "M  after : " << op1.M() << std::endl;
  }
  else {
    op1.status = 1;
    op1.M(-1); //FIXME
  }
  this->_ev->AddParticle(&op1, true);
  
  // Second outgoing proton
  Particle op2(5, _pdg5);
  plab_op2[0] = _pp5*_st5*_cp5;
  plab_op2[1] = _pp5*_st5*_sp5;
  plab_op2[2] = _gamma*_pp5*_ct5+_betgam*_ep5;
  plab_op2[3] = _gamma*_ep5     +_betgam*_pp5*_ct5;
  if (!op2.P( plab_op2[0]*cp+rany*plab_op2[1]*sp,
             -plab_op2[0]*sp+rany*plab_op2[1]*cp,
              plab_op2[2],
              plab_op2[3])) {
    std::cerr << "Invalid outgoing proton 2" << std::endl;
  }
  if (_cuts.kinematics==3) {
    op2.status = -2;
    op2.M(_mp5);
  }
  else {
    op2.status = 1;
    /*std::cout << "----> " << _pp5*_pp5-_ep5*_ep5 << ", " << _gamma << ", " << _betgam << std::endl;
    std::cout << "second proton left undecayed : m = " << op2.M() << std::endl;
    op2.Dump();*/
    op2.M(-1); //FIXME
  }
  this->_ev->AddParticle(&op2, true);

  // First incoming photon
  // Equivalent in LPAIR : PLAB(x, 3)
  Particle ph1(41, 22);
  plab_ph1[0] = plab_ip1[0]-plab_op1[0];
  plab_ph1[1] = plab_ip1[1]-plab_op1[1];
  plab_ph1[2] = plab_ip1[2]-plab_op1[2];
  plab_ph1[3] = plab_ip1[3]-plab_op1[3];
  if (!ph1.P( plab_ph1[0]*cp+rany*plab_ph1[1]*sp,
             -plab_ph1[0]*sp+rany*plab_ph1[1]*cp,
	      plab_ph1[2],
	      plab_ph1[3])) {
    //std::cerr << "Invalid photon 1" << std::endl;
  }
  ph1.charge = 0;
  ph1.status = -1;
  this->_ev->AddParticle(&ph1);
  
  // Second incoming photon
  // Equivalent in LPAIR : PLAB(x, 4)
  Particle ph2(42, 22);
  plab_ph2[0] = plab_ip2[0]-plab_op2[0];
  plab_ph2[1] = plab_ip2[1]-plab_op2[1];
  plab_ph2[2] = plab_ip2[2]-plab_op2[2];
  plab_ph2[3] = plab_ip2[3]-plab_op2[3];
  if (!ph2.P( plab_ph2[0]*cp+rany*plab_ph2[1]*sp,
             -plab_ph2[0]*sp+rany*plab_ph2[1]*cp,
              plab_ph2[2],
              plab_ph2[3])) {
    //std::cerr << "Invalid photon 2" << std::endl;
  }
  ph2.charge = 0;
  ph2.status = -1;
  this->_ev->AddParticle(&ph2);

  // Central (two-photon) system
  Particle cs(4, 22);
  cs.status = -1;
  this->_ev->AddParticle(&cs);
  
  // First outgoing lepton
  role = (ransign<0) ? 6 : 7;
  Particle ol1(role, ransign*abs(_pdg6));
  plab_ol1[0] = _pl6*_st6*_cp6;
  plab_ol1[1] = _pl6*_st6*_sp6;
  plab_ol1[2] = _gamma*_pl6*_ct6+_betgam*_el6;
  plab_ol1[3] = _gamma*_el6     +_betgam*_pl6*_ct6;
  if (!ol1.P( plab_ol1[0]*cp+rany*plab_ol1[1]*sp,
             -plab_ol1[0]*sp+rany*plab_ol1[1]*cp,
              plab_ol1[2],
              plab_ol1[3])) {
    std::cerr << "Invalid outgoing lepton 1" << std::endl;
  }
  ol1.charge = ransign;
  ol1.status = 1;
  ol1.M(-1); //FIXME
  this->_ev->AddParticle(&ol1);
  
  // Second outgoing lepton
  role = (ransign<0) ? 7 : 6;
  Particle ol2(role, -ransign*abs(_pdg7));
  plab_ol2[0] = _pl7*_st7*_cp7;
  plab_ol2[1] = _pl7*_st7*_sp7;
  plab_ol2[2] = _gamma*_pl7*_ct7+_betgam*_el7;
  plab_ol2[3] = _gamma*_el7     +_betgam*_pl7*_ct7;
  if (!ol2.P( plab_ol2[0]*cp+rany*plab_ol2[1]*sp,
             -plab_ol2[0]*sp+rany*plab_ol2[1]*cp,
	      plab_ol2[2],
              plab_ol2[3])) {
    std::cerr << "Invalid outgoing lepton 2" << std::endl;
  }
  ol2.charge = -ransign;
  ol2.status = 1;
  ol2.M(-1); //FIXME
  this->_ev->AddParticle(&ol2);

  // Relations between particles

  this->_ev->GetOneByRole(3)->SetMother(this->_ev->GetOneByRole(1));
  this->_ev->GetOneByRole(5)->SetMother(this->_ev->GetOneByRole(2));
  this->_ev->GetOneByRole(41)->SetMother(this->_ev->GetOneByRole(1));
  this->_ev->GetOneByRole(42)->SetMother(this->_ev->GetOneByRole(2));
  this->_ev->GetOneByRole(4)->SetMother(this->_ev->GetOneByRole(41));
  this->_ev->GetOneByRole(4)->SetMother(this->_ev->GetOneByRole(42));
  this->_ev->GetOneByRole(6)->SetMother(this->_ev->GetOneByRole(4));
  this->_ev->GetOneByRole(7)->SetMother(this->_ev->GetOneByRole(4));

  //std::cout << "---> " << __PRETTY_FUNCTION__ << ", " << __FILE__ << std::endl;

#ifdef DEBUG
  gmux = -_t2/(_ep1*_eg2-_pp1*_p3_g2[2])/2.;
  gmuy = (_ep1*plab_ph2[3]-_pp1*plab_ph2[2])/(_ep2*plab_ph2[3]+_pp2*plab_ph2[2]);
  gmuw = std::pow(_ep1+plab_ph2[3], 2)-std::pow(_pp1+plab_ph2[2], 2);
  if (gmuw>=0.) {
    gmuw = std::sqrt(gmuw);
  }
  else {
    std::cerr << "[GamGamLL::FillKinematics] [FATAL] W**2 = " << gmuw << " < 0" << std::endl;
    gmuw = 0.;
  }
  gmunu = gmuy*2.*GetMassFromPDGId(2212)/_ep1/_ep2;
  std::cout << "[GamGamLL::FillKinematics] [DEBUG]" << std::endl
            << "   gmux = " << gmux << std::endl
            << "   gmuy = " << gmuy << std::endl
            << "   gmuw = " << gmuw << std::endl
            << "  gmunu = " << gmunu << std::endl;
#endif
  //this->_ev->Dump();
}

void
GamGamLL::SetKinematics(Kinematics cuts_)
{
  _cuts = cuts_;
  _cotth1 = 1./tan(cuts_.thetamax*pi/180.);
  _cotth2 = 1./tan(cuts_.thetamin*pi/180.);  
#ifdef DEBUG
  std::cout << "[GamGamLL::SetKinematics] [DEBUG]" << std::endl
            << "  cot(theta1) = " << _cotth1 << std::endl
            << "  cot(theta2) = " << _cotth2 << std::endl;
#endif
}

double
GamGamLL::PeriPP(int nup_, int ndown_)
{
  double rho = .585; // 1.05
  double cc1 = .86926; // .6303
  double cc2 = 2.23422; // 2.2049
  double dd1 = .12549; // .0468
  double cp = .96; // 1.23
  double bp = .63; // .61

  double dummy, psfw1, psfw2;
  double en, x, xt;
  double rhot, qqq, qdq;
  double t11, t12, t21, t22;
  double peripp;

  // Implicitely needed : full event kinematics, such as :
  // * _t1, _t2
  // * _mp1, _mp2 (--> _w1, _w2)
  // * _mp3, _mp5 (--> _w3, _w5)

#ifdef DEBUG
  std::cout << "[GamGamLL::PeriPP] [DEBUG]" << std::endl
	    << "   Nup  = " << nup_ << std::endl
	    << "  Ndown = " << ndown_ << std::endl;
#endif

  switch(nup_) {
  case 1: // DESY
      //case 0: // CDF
    _u1 = 1.;
    _u2 = 1.;
    break;
  case 2: // DESY
    //case 1: // CDF
    xt = std::pow(1.-_t1/.71, 2);
    _tau = _t1/(4.*_w1);
    _u1 = std::pow(2.79/xt, 2);
    _u2 = (1./std::pow(xt, 2)-_u1*_tau)/(1.-_tau);
    break;
  case 4: // DESY
    // does not exist in CDF version
    PSF(_t1, _w3, &dummy, &psfw1, &psfw2);
#ifdef DEBUG
    std::cout << "[GamGamLL::PeriPP] [DEBUG] Result of PSF : " << PSF(_t1, _w3, &dummy, &psfw1, &psfw2) << std::endl;
#endif
    _u1 = -psfw1*(2.*_mp1)/_t1;
    _u2 = psfw2/(2.*_mp1);
    break;
  default:
    x = _t1/(_t1-_w3);
    en = _w31-_t1;
    _tau = _t1/(4.*_w1);
    rhot = rho-_t1;
    _u1 = (-cc1*std::pow(rho/rhot, 2)*_w31-cc2*_w1*std::pow(1.-x, 4)/(x*(x*cp-2*bp)+1.))/_t1;
    _u2 = (-_tau*_u1-dd1*_w31*_t1*(rho/rhot)*std::pow(_w31/en, 2)/(rhot*_w1))/(1.-std::pow(en, 2)/(4.*_w1*_t1));
    break;
  }

  switch(ndown_) {
  case 1: // DESY
    //case 0: // CDF
    _v1 = 1.;
    _v2 = 1.;
    break;
  case 2: // DESY
    //case 1: // CDF
    xt = std::pow(1.-_t2/.71, 2);
    _tau = _t2/(4.*_w2);
    _v1 = std::pow(2.79/xt, 2);
    _v2 = (1./std::pow(xt, 2)-_v1*_tau)/(1.-_tau);
    break;
  default:
    x = _t2/(_t2-_w5);
    en = _w52-_t2;
    _tau = _t2/(4.*_w2);
    rhot = rho-_t2;
    _v1 = (-cc1*std::pow(rho/rhot, 2)*_w52-cc2*_w2*std::pow(1.-x, 4)/(x*(x*cp-2.*bp)+1))/_t2;
    _v2 = (-_tau*_v1-dd1*_w52*_t2*(rho/rhot)*std::pow(_w52/en, 2)/(rhot*_w2))/(1.-std::pow(en, 2)/(4.*_w2*_t2));
    break;
  }
#ifdef DEBUG
  std::cout << "[GamGamLL::PeriPP] [DEBUG]" << std::endl
            << "  u1 = " << _u1 << std::endl
            << "  u2 = " << _u2 << std::endl
            << "  v1 = " << _v1 << std::endl
            << "  v2 = " << _v2 << std::endl;
#endif

  qqq = std::pow(_q1dq, 2);
  qdq = 4.*_w6-_w4;
  t11 = 64.*(_bb*(qqq-_g4-qdq*(_t1+_t2+2.*_w6))-2.*(_t1+2.*_w6)*(_t2+2.*_w6)*qqq)*_t1*_t2;
  t12 = 128.*(-_bb*(_dd2+_g6)-2.*(_t1+2.*_w6)*(_sa2*qqq+std::pow(_a6, 2)))*_t1;
  t21 = 128.*(-_bb*(_dd4+_g5)-2.*(_t2+2.*_w6)*(_sa1*qqq+std::pow(_a5, 2)))*_t2;
  t22 = 512.*(_bb*(std::pow(_delta, 2)-_gram)-std::pow(_epsi-_delta*(qdq+_q1dq2), 2)-_sa1*std::pow(_a6, 2)-_sa2*std::pow(_a5, 2)-_sa1*_sa2*qqq);

  peripp = (((_u1*_v1*t11+_u2*_v1*t21+_u1*_v2*t12+_u2*_v2*t22)/(_t1*_t2*_bb))/(_t1*_t2*_bb))/4.;

#ifdef DEBUG
  std::cout << "[GamGamLL::PeriPP] [DEBUG]" << std::endl
            << "  t11 = " << t11 << std::endl
            << "  t12 = " << t12 << std::endl
            << "  t21 = " << t21 << std::endl
            << "  t22 = " << t22 << std::endl
            << "  tau = " << _tau << std::endl
            << "  --> PeriPP = " << peripp << std::endl;
#endif
  return peripp;
}
