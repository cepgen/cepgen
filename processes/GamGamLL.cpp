#include "GamGamLL.h"

GamGamLL::GamGamLL(int nOpt_) : GenericProcess("pp -> p(*) (gamma gamma -> l+ l-) p(*)"),
  _nOpt(nOpt_),
  _ep1(-1), _w1(-1), _ep2(-1), _w2(-1.),
  _w3(-1.), _w4(-1.), _w5(-1.), _w6(-1.), _w7(-1.),
  _p12(0.), _p13(0.), _p14(0.), _p15(0.), _p23(0.), _p24(0.), _p25(0.), _p34(0.), _p35(0.), _p45(0.),
  _p1k2(0.), _p2k1(0.),
  setp1(false), setp2(false), setp3(false), setp5(false), setll(false),
  _cotth1(-99999.), _cotth2(99999.)
{}

void
GamGamLL::AddEventContent()
{
  IncomingState is; OutgoingState os;
  is.insert(ParticleWithRole(Particle::IncomingBeam1,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::IncomingBeam2,    Particle::Proton));
  is.insert(ParticleWithRole(Particle::Parton1,          Particle::Photon));
  is.insert(ParticleWithRole(Particle::Parton2,          Particle::Photon));
  os.insert(ParticleWithRole(Particle::OutgoingBeam1,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::OutgoingBeam2,    Particle::Proton));
  os.insert(ParticleWithRole(Particle::CentralParticle1, Particle::Muon));
  os.insert(ParticleWithRole(Particle::CentralParticle2, Particle::Muon));
  GenericProcess::SetEventContent(is, os);
}

int
GamGamLL::GetNdim(ProcessMode process_mode_) const
{
  switch (process_mode_) {
    case ElasticElastic: default:                 return 7;
    case ElasticInelastic: case InelasticElastic: return 8;
    case InelasticInelastic:                      return 9;
  }
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
  
  DebugInsideLoop(Form("Optimised mode? %i", _nOpt));
  
  fJacobian = 0.;

  _w4 = std::pow(_mc4, 2);

  sig = _mc4+_mp5;
  sig1 = std::pow(sig, 2);
  sig2 = std::pow(sig, 2);

  DebugInsideLoop(Form("mc4 = %f\n\t"
                       "sig1 = %f\n\t"
                       "sig2 = %f", _mc4, sig1, sig2));

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

  DebugInsideLoop(Form("w1 = %f\n\t"
                       "w2 = %f\n\t"
                       "w3 = %f\n\t"
                       "w4 = %f\n\t"
                       "w5 = %f",
                       _w1, _w2, _w3, _w4, _w5));

  //Info(Form("w31 = %f\n\tw52 = %f\n\tw12 = %f", _w31, _w52, _w12));
  
  ss = fS+_w12;
  
  rl1 = std::pow(ss, 2)-4.*_w1*fS; // lambda(s, m1**2, m2**2)
  if (rl1<=0.) { Warning(Form("rl1 = %f <= 0", rl1)); return false; }
  _sl1 = std::sqrt(rl1);

  _s2 = ds2 = 0.;
  if (_nOpt==0) {
    smax = fS+_w3-2.*_mp3*fSqS;
    Map(x(2), sig1, smax, &_s2, &ds2, "s2");
    sig1 = _s2; //FIXME!!!!!!!!!!!!!!!!!!!!
  }
  
  DebugInsideLoop(Form("s2 = %f", _s2));

  //std::cout << "s=" << _s << ", w3=" << _w3 << ", sig1=" << sig1 << std::endl;
  sp = fS+_w3-sig1;

  _d3 = sig1-_w2;

  rl2 = std::pow(sp, 2)-4.*fS*_w3; // lambda(s, m3**2, sigma)
  if (rl2<=0.) { Warning(Form("rl2 = %f <= 0", rl2)); return false; }
  sl2 = std::sqrt(rl2);

  //std::cout << "ss=" << ss << ", sp=" << sp << ", sl1=" << _sl1 << ", sl2=" << sl2 << std::endl;
  _t1max = _w1+_w3-(ss*sp+_sl1*sl2)/(2.*fS); // definition from eq. (A.4) in [1]
  _t1min = (_w31*_d3+(_d3-_w31)*(_d3*_w1-_w31*_w2)/fS)/_t1max; // definition from eq. (A.5) in [1]

  // FIXME dropped in CDF version
  if (_t1max>-fCuts.q2min) { Warning(Form("t1max = %f > -q2min = %f", _t1max, -fCuts.q2min)); return false; }
  if (_t1min<-fCuts.q2max and fCuts.q2max>=0.) { Warning(Form("t1min = %f < -q2max = %f", _t1min, -fCuts.q2max)); return false; }
  if (_t1max<-fCuts.q2max and fCuts.q2max>=0.) _t1max = -fCuts.q2max;
  if (_t1min>-fCuts.q2min)                     _t1min = -fCuts.q2min;
  /////

  // t1, the first photon propagator, is defined here
  Map(x(0), _t1min, _t1max, &_t1, &dt1, "t1");
  // changes wrt mapt1 : dx->-dx
  dt1 = -dt1;
  
  DebugInsideLoop(Form("Definition of t1 = %f according to\n\t"
                       "(t1min, t1max) = (%f, %f)", _t1, _t1min, _t1max));

  _dd4 = _w4-_t1;
  d8 = _t1-_w2;

  t13 = _t1-_w1-_w3;

  _sa1 =-std::pow(_t1-_w31, 2)/4.+_w1*_t1;
  if (_sa1>=0.) {
    Warning(Form("_sa1 = %f >= 0", _sa1)); return false;
  }

  sl3 = std::sqrt(-_sa1);

  // one computes splus and (s2x=s2max)
  if (_w1!=0.) {
    sb =(fS*(_t1-_w31)+_w12*t13)/(2.*_w1)+_w3;
    sd = _sl1*sl3/_w1;
    se =(fS*(_t1*(fS+t13-_w2)-_w2*_w31)+_w3*(_w12*d8+_w2*_w3))/_w1;
    if (fabs((sb-sd)/sd)>=1.) {
      splus = sb-sd;
      s2max = se/splus;
    }
    else {
      s2max = sb+sd;
      splus = se/s2max;
    }
  }
  else { // 3
    s2max = (fS*(_t1*(fS+d8-_w3)-_w2*_w3)+_w2*_w3*(_w2+_w3-_t1))/(ss*t13);
    splus = sig2;
  }
  // 4
  s2x = s2max;
  
  DebugInsideLoop(Form("s2x = s2max = %f", s2x));
  
  if (_nOpt<0) { // 5
    if (splus>sig2) {
      sig2 = splus;
      DebugInsideLoop(Form("sig2 truncated to splus = %f", splus));
    }
    if (_nOpt<-1) {
      Map(x(2), sig2, s2max, &_s2, &ds2, "s2");
    }
    else { // _nOpt==-1
      Mapla(_t1, _w2, x(2), sig2, s2max, &_s2, &ds2);
    }
    s2x = _s2;
  }
  else if (_nOpt==0) { // 6
    s2x = _s2;
  }

  DebugInsideLoop(Form("s2x = %f", s2x));
  
  // 7
  r1 = s2x-d8;
  r2 = s2x-d6;
  
  rl4 = (std::pow(r1, 2)-4.*_w2*s2x)*(std::pow(r2, 2)-4.*_w5*s2x);
  if (rl4<=0.) { Warning(Form("rl4 = %f <= 0", rl4)); return false; }
  sl4 = std::sqrt(rl4);
  
  // t2max, t2min definitions from eq. (A.12) and (A.13) in [1]
  _t2max = _w2+_w5-(r1*r2+sl4)/(2.*s2x);
  _t2min = (_w52*_dd4+(_dd4-_w52)*(_dd4*_w2-_w52*_t1)/s2x)/_t2max;

  // t2, the second photon propagator, is defined here
  Map(x(1), _t2min, _t2max, &_t2, &dt2, "t2");
  // changes wrt mapt2 : dx->-dx

  dt2 = -dt2;

  _tau = _t1-_t2;
  r3 = _dd4-_t2;
  r4 = _w52-_t2;

  DebugInsideLoop(Form("r1 = %f\n\tr2 = %f\n\tr3 = %f\n\tr4 = %f", r1, r2, r3, r4));

  b = r3*r4-2.*(_t1+_w2)*_t2;
  c = _t2*d6*d8+(d6-d8)*(d6*_w2-d8*_w5);
  t25 = _t2-_w2-_w5;

  _sa2 = -std::pow(r4, 2)/4.+_w2*_t2;
  if (_sa2>=0.) { Warning(Form("_sa2 = %f >= 0", _sa2)); return false; }
  
  sl6 = 2.*std::sqrt(-_sa2);
  
  _g4 = -std::pow(r3, 2)/4.+_t1*_t2;
  if (_g4>=0.)  { Warning(Form("_g4 = %f >= 0", _g4)); return false; }
  
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
  if (_nOpt>1)       Map(x(2), s2min, s2max, &_s2, &ds2, "s2");
  else if (_nOpt==1) Mapla(_t1, _w2, x(2), s2min, s2max, &_s2, &ds2);

  ap = -std::pow(_s2+d8, 2)/4.+_s2*_t1;
  //Info(Form("s2 = %f, s2max = %f, splus = %f", _s2, s2max, splus));
  if (_w1!=0.) _dd1 = -_w1*(_s2-s2max)*(_s2-splus)/4.; // 10
  else         _dd1 = ss*t13*(_s2-s2max)/4.;
  // 11
  _dd2 = -_t2*(_s2-s2p)*(_s2-s2min)/4.;
  //Info(Form("t2 =%f\n\ts2 =%f\n\ts2p=%f\n\ts2min=%f\n\tdd2=%f",_t2, _s2, s2p, s2min, _dd2));
  // FIXME there should be a more beautiful way to check for nan!
  // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
  // FIXME dropped in CDF version
  if (!(_dd2>=0.) and !(_dd2<0.)) { // NaN
    Warning(Form("dd2 == NaN\n\tdd2 = %f\n\t"
                 "s2 = %f\ts2p = %f\n\t"
                 "s2min = %f\n\t"
                 "t2min = %f\tt2max = %f",
                 _dd2, _s2, s2p, s2min, _t2min, _t2max));
    return false;
  }
  
  yy4 = cos(pi*x(3));
  dd = _dd1*_dd2;
  _p12 = (fS-_w1-_w2)/2.;
  st = _s2-_t1-_w2;
  delb = (2.*_w2*r3+r4*st)*(4.*_p12*_t1-(_t1-_w31)*st)/(16.*ap);

  if (dd<=0.) { Warning(Form("dd = %f <= 0\n\tdd1 = %f\tdd2 = %f", dd, _dd1, _dd2)); return false; }

  _delta = delb-yy4*st*std::sqrt(dd)/(2.*ap);
  _s1 = _t2+_w1+(2.*_p12*r3-4.*_delta)/st;

  if (ap>=0.) { Warning(Form("ap = %f >= 0", ap)); return false; }

  fJacobian = ds2*dt1*dt2*std::pow(pi, 2)/(8.*_sl1*std::sqrt(-ap));

  DebugInsideLoop(Form("dj=%f", fJacobian));

  _gram = (1.-std::pow(yy4, 2))*dd/ap;

  _p13 = -t13/2.;
  _p14 = (_tau+_s1-_w3)/2.;
  _p15 = (fS+_t2-_s1-_w2)/2.;
  _p23 = (fS+_t1-_s2-_w1)/2.;
  _p24 = (_s2-_tau-_w5)/2.;
  _p25 = -t25/2.;
  _p34 = (_s1-_w3-_w4)/2.;
  _p35 = (fS+_w4-_s1-_s2)/2.;
  _p45 = (_s2-_w4-_w5)/2.;

  _p1k2 = (_s1-_t2-_w1)/2.;
  _p2k1 = st/2.;


  if (_w2!=0.) {
    sbb = (fS*(_t2-_w52)-_w12*t25)/(2.*_w2)+_w5;
    sdd = _sl1*sl6/(2.*_w2);
    see = (fS*(_t2*(fS+t25-_w1)-_w1*_w52)+_w5*(_w1*_w5-_w12*(_t2-_w1)))/_w2;
    if (sbb/sdd>=0.) {
      s1p = sbb+sdd;
      s1m = see/s1p;
      // FIXME there should be a more beautiful way to check for nan!
      // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
      // FIXME dropped in CDF version
      if (!(_dd2>=0.) && !(_dd2<0.)) { // NaN
        Warning(Form("dd2 == NaN\n\tdd2 = %f\n\t"
                     "s1 = %f\ts1p = %f\ts1m = %f\n\t"
                     "w2 = %f",
                     _dd2, _s1, s1p, s1m, _w2));
        return false;
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
    s1p = (fS*(_t2*(fS-_w5+_t2-_w1)-_w1*_w5)+_w1*_w5*(_w1+_w5-_t2))/(t25*(fS-_w12));
    _dd3 = -t25*(fS-_w12)*(s1p-_s1)/4.;
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
  if (!Pickin() or fJacobian==0.) {
    Warning(Form("Pickin failed! dj = %f", fJacobian));
    return false;
  }
  
  const double re = 1./(2.*fSqS);
  _ep1 = re*(fS+_w12);
  _ep2 = re*(fS-_w12);

  DebugInsideLoop(Form(" re = %f\n\t_w12 = %f", re, _w12));
  DebugInsideLoop(Form("Incoming particles' energy = %f, %f", _ep1, _ep2));

  _p = re*_sl1;

  _de3 = re*(_s2-_w3+_w12);
  _de5 = re*(_s1-_w5-_w12);

  // Final state energies
  const double ep3 = _ep1-_de3, ep5 = _ep2-_de5;
  _ec4 = _de3+_de5;

  if (_ec4<_mc4) {
    Warning(Form("_ec4 = %f < _mc4 = %f\n\t==> de3 = %f, de5 = %f", _ec4, _mc4, _de3, _de5));
    return false;
  }
  // What if the protons' momenta are not along the z-axis?
  _pc4 = std::sqrt((std::pow(_ec4, 2)-std::pow(_mc4, 2)));

  if (_pc4==0.) {
    Warning("_pzc4==0");
    return false;
  }
  const double pp3 = std::sqrt(std::pow(ep3, 2)-_w3), pt3 = std::sqrt(_dd1/fS)/_p,
               pp5 = std::sqrt(std::pow(ep5, 2)-_w5), pt5 = std::sqrt(_dd3/fS)/_p;

  DebugInsideLoop(Form("Central system's energy: E4 = %f\n\t"
                       "                 momentum: p4 = %f\n\t"
                       "                 invariant mass: m4 = %f\n\t"
                       "Outgoing particles' energy: E3 = %f\n\t"
                       "                            E5 = %f",
                       _ec4, _pc4, _mc4, ep3, ep5));

  const double st3 = pt3/pp3,
               st5 = pt5/pp5;

  DebugInsideLoop(Form("st3 = %f\n\tst5 = %f", st3, st5));

  // FIXME there should be a more beautiful way to check for nan!
  // (http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c)
  // FIXME dropped in CDF version
  if (!(_dd3>=0.) && !(_dd3<0.)) { Error("dd3 == NaN"); } // NaN
  if (!(_dd1>=0.) && !(_dd1<0.)) { Error("dd1 == NaN"); } // NaN
  /////

  if (st3>1.) { Warning(Form("st3 = %f > 1", st3)); return false; }
  if (st5>1.) { Warning(Form("st5 = %f > 1", st5)); return false; }

  double ct3 = std::sqrt(1.-std::pow(st3, 2)),
         ct5 = std::sqrt(1.-std::pow(st5, 2));

  if (_ep1*ep3<_p13) ct3 *= -1.;
  if (_ep2*ep5>_p25) ct5 *= -1.;

  DebugInsideLoop(Form("ct3 = %f\n\tct5 = %f", ct3, ct5));
  
  if (_dd5<0.) { Warning(Form("dd5 = %f < 0", _dd5)); return false; }

  // Centre of mass system kinematics (theta4 and phi4)
  _pt4 = std::sqrt(_dd5/fS)/_p;
  _st4 = _pt4/_pc4;

  if (_st4>1.) { Warning(Form("st4 = %f > 1", _st4)); return false; }
  
  _ct4 = std::sqrt(1.-std::pow(_st4, 2));
  if (_ep1*_ec4<_p14) _ct4 *= -1.;

  _al4 = 1.-_ct4;
  _be4 = 1.+_ct4;

  if (_ct4<0.) _be4 = std::pow(_st4, 2)/_al4;
  else         _al4 = std::pow(_st4, 2)/_be4;

  DebugInsideLoop(Form("ct4 = %f\n\tal4 = %f, be4 = %f", _ct4, _al4, _be4));

  const double rr  = std::sqrt(-_gram/fS)/(_p*_pt4);
  const double sp3 = rr/pt3, sp5 = -rr/pt5;
  //std::cout << "rr = " << rr << ", sp3 = " << fabs(sp3) << ", sp5 = " << fabs(sp5) << std::endl;

  if (fabs(sp3)>1.) { Warning(Form("sp3 = %f > 1", sp3)); return false; }
  if (fabs(sp5)>1.) { Warning(Form("sp5 = %f > 1", sp5)); return false; }

  const double cp3 = -std::sqrt(1.-std::pow(sp3, 2)), cp5 = -std::sqrt(1.-std::pow(sp5, 2));
  
  fP3lab = Particle::Momentum(pp3*st3*cp3, pp3*st3*sp3, pp3*ct3, ep3);
  fP5lab = Particle::Momentum(pp5*st5*cp5, pp5*st5*sp5, pp5*ct5, ep5);

  const double a1 = fP3lab.Px()-fP5lab.Px();

  DebugInsideLoop(Form("Kinematic quantities\n\t"
                       "cos(theta3) = %1.4f\tsin(theta3) = %1.4f\tcos( phi3 ) = %1.4f\tsin( phi3 ) = %1.4f\n\t"
                       "cos(theta4) = %1.4f\tsin(theta4) = %1.4f\n\t"
                       "cos(theta5) = %1.4f\tsin(theta5) = %1.4f\tcos( phi5 ) = %1.4f\tsin( phi5 ) = %1.4f\n\t"
                       "a1 = %f",
                       ct3, st3, cp3, sp3,
                       _ct4, _st4,
                       ct5, st5, cp5, sp5, a1));

  if (fabs(_pt4+fP3lab.Px()+fP5lab.Px())<fabs(fabs(a1)-_pt4)) {
    DebugInsideLoop(Form("|pp4+pp3*cos(phi3)+pp5*cos(phi5)| < | |a1|-pp4 |\n\t"
                         "pp4 = %f\tpp5 = %f\n\t"
                         "cos(phi3) = %f\tcos(phi5) = %f"
                         "a1 = %f",
                         _pt4, pt5, cp3, cp5, a1));
    return true;
  }
  /*if (a1<0.) _cp5 = -_cp5;
  else       _cp3 = -_cp3;*/
  if (a1<0.) fP5lab.SetP(0, -fP5lab.Px());
  else       fP3lab.SetP(0, -fP3lab.Px());
  return true;
}

double
GamGamLL::ComputeMX(double x_, double outmass_, double lepmass_, double *dw_)
{
  double mx2, dmx2;

  const double wx2min = std::pow(std::max(Particle::GetMassFromPDGId(Particle::Proton)+Particle::GetMassFromPDGId(Particle::PiPlus), fCuts.mxmin), 2),
               wx2max = std::pow(std::min(fSqS-outmass_-2.*lepmass_, fCuts.mxmax), 2);
  
  Map(x_, wx2min, wx2max, &mx2, &dmx2, "mx2");

  DebugInsideLoop(Form("mX^2 in range (%f, %f), x = %f\n\t"
                       "mX^2 = %f, d(mX^2) = %f\n\t"
                       "mX = %f, d(mX) = %f", wx2min, wx2max, x_, mx2, dmx2, sqrt(mx2), sqrt(dmx2)));

  *dw_ = sqrt(dmx2);
  return sqrt(mx2);
}

void
GamGamLL::BeforeComputeWeight()
{
  if (!GenericProcess::fIsPointSet) return;
  
  Particle *p1 = GetParticle(Particle::IncomingBeam1),
           *p2 = GetParticle(Particle::IncomingBeam2);

  _ep1 = p1->E();
  _mp1 = p1->M();
  _w1 = p1->M2();
  _pp1 = p1->GetMomentum().P();
  setp1 = true;

  _ep2 = p2->E();
  _mp2 = p2->M();
  _w2 = p2->M2();
  _pp2 = p2->GetMomentum().P();

  const double thetamin = EtaToTheta(fCuts.etamax),
               thetamax = EtaToTheta(fCuts.etamin);
  _cotth1 = 1./tan(thetamax*pi/180.);
  _cotth2 = 1./tan(thetamin*pi/180.);
  DebugInsideLoop(Form("cot(theta1) = %f\n\tcot(theta2) = %f", _cotth1, _cotth2));

  Particle* p;
  p = GetParticle(Particle::CentralParticle1); _w6 = p->M2();
  p = GetParticle(Particle::CentralParticle2); _w7 = p->M2();
  
  double m;
  switch (fCuts.kinematics) {
  case GenericProcess::ElasticElastic:
    _dw31 = _dw52 = 0.; break;
  case GenericProcess::ElasticInelastic: case GenericProcess::InelasticElastic:
    m = ComputeMX(x(7), GetParticle(Particle::IncomingBeam1)->M(), GetParticle(Particle::CentralParticle1)->M(), &_dw31);
    GetParticle(Particle::OutgoingBeam1)->SetM(m);
    GetParticle(Particle::OutgoingBeam2)->SetM(Particle::GetMassFromPDGId(GetParticle(Particle::OutgoingBeam2)->GetPDGId())); //FIXME
    break;
  case GenericProcess::InelasticInelastic:
    m = ComputeMX(x(7), GetParticle(Particle::IncomingBeam2)->M(), GetParticle(Particle::CentralParticle1)->M(), &_dw31);
    GetParticle(Particle::OutgoingBeam1)->SetM(m);
    m = ComputeMX(x(8), GetParticle(Particle::OutgoingBeam1)->M(), GetParticle(Particle::CentralParticle1)->M(), &_dw52);
    GetParticle(Particle::OutgoingBeam2)->SetM(m);
    break;
  }
  _mp3 = GetParticle(Particle::OutgoingBeam1)->M();
  _mp5 = GetParticle(Particle::OutgoingBeam2)->M();
  _w3 = std::pow(_mp3, 2);
  _w5 = std::pow(_mp5, 2);
  
  /*if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) {
    IsKinematicsDefined();
    if (fIsOutStateSet)  std::cout << "  --> Outgoing state is fully set" << std::endl;
    if (fIsKinematicSet) std::cout << "  --> Kinematics is fully set" << std::endl;
  }*/
}

double
GamGamLL::ComputeWeight()
{
  // WARNING ====> PP5-->_pt5
  //                P5-->_pp5
  double weight;
  double dw4;
  double wmin, wmax;
  double stg, cpg, spg, ctg;

  double qcx, qcz;
  double pc6x, pc6z;

  double b1, b2, b3;
  double c1, c2, c3;
  double h1, h2, hq;
  double r12, r13, r22, r23;

  bool lcut, lmu1, lmu2;

  weight = 0.;

  if (!fIsOutStateSet) { Warning("Output state not set!"); return 0.; }

  if (fCuts.wmax<0) fCuts.wmax = fS;

  // The minimal energy for the central system is its outgoing leptons' mass energy (or wmin_ if specified)
  wmin = std::pow(GetParticle(Particle::CentralParticle1)->M()+GetParticle(Particle::CentralParticle2)->M(),2);
  if (fabs(wmin)<fabs(fCuts.wmin)) wmin = fCuts.wmin;

  // The maximal energy for the central system is its CM energy with the outgoing particles' mass energy substracted (or _wmax if specified)
  wmax = std::pow(fSqS-_mp3-_mp5,2);
  DebugInsideLoop(Form("sqrt(s)=%f\n\tm(X1)=%f\tm(X2)=%f", fSqS, _mp3, _mp5));
  if (fabs(wmax)>fabs(fCuts.wmax)) wmax = fCuts.wmax;
  
  DebugInsideLoop(Form("wmin = %f\n\twmax = %f\n\twmax/wmin = %f", wmin, wmax, wmax/wmin));
  
  Map(x(4), wmin, wmax, &_w4, &dw4, "w4");
  _mc4 = std::sqrt(_w4);
  
  DebugInsideLoop(Form("Computed value for w4 = %f -> mc4 = %f", _w4, _mc4));
  
  if (!Orient()) return 0.;

  if (_t1>0.)  { Warning(Form("t1 = %f > 0", _t1)); return 0.; }
  if (_t2>0.)  { Warning(Form("t2 = %f > 0", _t2)); return 0.; }
  if (fJacobian==0.) { Warning(Form("dj = %f", fJacobian)); return 0.; }
  
  const double ecm6 = (_w4+_w6-_w7)/(2.*_mc4);
  const double pp6cm = std::sqrt(std::pow(ecm6, 2)-_w6);
  
  fJacobian *= dw4*pp6cm/(_mc4*sconstb*fS);
  
  // Let the most obscure part of this code begin...

  const double e1mp1 = _w1/(_ep1+_p),
               e3mp3 = _w3/(fP3lab.E()+fP3lab.P());

  const double al3 = std::pow(sin(fP3lab.Theta()), 2)/(1.+(fP3lab.Theta()));

  // 2-photon system kinematics ?!
  const double eg = (_w4+_t1-_t2)/(2.*_mc4);
  double pg = std::sqrt(std::pow(eg, 2)-_t1);
  const double pgx = -fP3lab.Px()*_ct4-_st4*(_de3-e1mp1+e3mp3+fP3lab.P()*al3),
               pgy = -fP3lab.Py(),
               pgz = _mc4*_de3/(_ec4+_pc4)-_ec4*_de3*_al4/_mc4-fP3lab.Px()*_ec4*_st4/_mc4+_ec4*_ct4/_mc4*(fP3lab.P()*al3+e3mp3-e1mp1);
  
  DebugInsideLoop(Form("pg3 = (%f, %f, %f)\n\t"
                       "pg3^2 = %f",
                       pgx, pgy, pgz,
                       std::sqrt(std::pow(pgx, 2)+std::pow(pgy, 2)+std::pow(pgz, 2))
                      ));
  
  const double pgp = std::sqrt(std::pow(pgx, 2)+std::pow(pgy, 2)), // outgoing proton (3)'s transverse momentum
               pgg = std::sqrt(std::pow(pgp, 2)+std::pow(pgz, 2)); // outgoing proton (3)'s momentum
  if (pgg>pgp*.9 && pgg>pg) { pg = pgg; } //FIXME ???
  
  // Phi angle for the 2-photon system ?!
  cpg = pgx/pgp;
  spg = pgy/pgp;
  
  // Theta angle for the 2-photon system ?!
  stg = pgp/pg;
  ctg = std::sqrt(1.-std::pow(stg, 2));
  if (pgz<0.) ctg *= -1.;
  
  double xx6 = x(5);
  
  const double amap = (_w4-_t1-_t2)/2.;
  const double bmap = std::sqrt((std::pow(_w4-_t1-_t2, 2)-4.*_t1*_t2)*(1.-4.*_w6/_w4))/2.;
  const double ymap = (amap+bmap)/(amap-bmap);
  const double beta = std::pow(ymap, (double)(2.*xx6-1.));
  xx6 = (amap/bmap*(beta-1.)/(beta+1.)+1.)/2.;
  if (xx6>1.) xx6 = 1.;
  if (xx6<0.) xx6 = 0.;
  
  DebugInsideLoop(Form("amap = %f\n\tbmap = %f\n\tymap = %f\n\tbeta = %f", amap, bmap, ymap, beta));
  
  // 3D rotation of the first outgoing lepton wrt the CM system
  const double theta6cm = acos(1.-2.*xx6);
  
  // match the Jacobian
  fJacobian *= (((amap+bmap*cos(theta6cm))*(amap-bmap*cos(theta6cm))/amap/bmap*log(ymap))/2.);
  
  DebugInsideLoop(Form("ctcm6 = %f\n\tstcm6 = %f", cos(theta6cm), sin(theta6cm)));
  
  const double phi6cm = 2.*pi*x(6);

  // First outgoing lepton's 3-momentum in the centre of mass system
  Particle::Momentum p6cm = Particle::Momentum::FromPThetaPhi(pp6cm, theta6cm, phi6cm);
  
  DebugInsideLoop(Form("p3cm6 = (%f, %f, %f)", p6cm.Px(), p6cm.Py(), p6cm.Pz()));

  h1 = stg*p6cm.Pz()+ctg*p6cm.Px();

  pc6z = ctg*p6cm.Pz()-stg*p6cm.Px();
  pc6x = cpg*h1-spg*p6cm.Py();
  
  qcx = 2.*pc6x;
  qcz = 2.*pc6z;
  // qcy == QCY is never defined
  
  double p6x, p6y, p6z;
  const double el6 = (_ec4*ecm6+_pc4*pc6z)/_mc4;
  h2 = (_ec4*pc6z+_pc4*ecm6)/_mc4;

  DebugInsideLoop(Form("h1 = %f\n\th2 = %f", h1, h2));

  // First outgoing lepton's 3-momentum
  p6x = _ct4*pc6x+_st4*h2;
  p6y = cpg*p6cm.Py()+spg*h1;
  p6z = _ct4*h2-_st4*pc6x;
  
  // first outgoing lepton's kinematics
  fP6cm = Particle::Momentum(p6x, p6y, p6z, el6);
  DebugInsideLoop(Form("E6(cm) = %f\n\tP6(cm) = (%f, %f, %f)", el6, p6x, p6y, p6z));
  
  hq = _ec4*qcz/_mc4;
  
  const Particle::Momentum qve(
    _ct4*qcx+_st4*hq,
    2.*p6y,
    _ct4*hq-_st4*qcx,
    _pc4*qcz/_mc4 // energy
  );
  
  // Available energy for the second lepton is the 2-photon system's energy with the first lepton's energy removed
  const double el7 = _ec4-el6;

  DebugInsideLoop(Form("Outgoing kinematics\n\t"
                       " first outgoing lepton: p = %f, E = %f\n\t"
                       "second outgoing lepton: p = %f, E = %f",
                       fP6cm.P(), fP6cm.E(), fP7cm.P(), fP6cm.E()));

  double p7x, p7y, p7z;
  // Second outgoing lepton's 3-momentum
  p7x = _pt4-p6x;
  p7y = -p6y;
  p7z = _pc4*_ct4-p6z;
  
  // second outgoing lepton's kinematics
  fP7cm = Particle::Momentum(p7x, p7y, p7z, el7);

  //fP6cm = Particle::Momentum(pl6*st6*cp6, pl6*st6*sp6, pl6*ct6, el6);
  //fP7cm = Particle::Momentum(pl7*st7*cp7, pl7*st7*sp7, pl7*ct7, el7);

  _q1dq = eg*(2.*ecm6-_mc4)-2.*pg*p6cm.Pz();
  _q1dq2 = (_w4-_t1-_t2)/2.;

  const double phi3 = fP3lab.Phi(), cp3 = cos(phi3), sp3 = sin(phi3),
               phi5 = fP5lab.Phi(), cp5 = cos(phi5), sp5 = sin(phi5);
  //std::cout << ">>> " << fP3lab.Pt() << "/" << fP5lab.Pt() << std::endl;

  _bb = _t1*_t2+(_w4*std::pow(sin(theta6cm), 2)+4.*_w6*std::pow(cos(theta6cm), 2))*std::pow(pg, 2);
  
  c1 = (qve.Px()*sp3-qve.Py()*cp3)*fP3lab.Pt();
  c2 = (qve.Pz()*_ep1-qve.E()*_p)*fP3lab.Pt();
  c3 = (_w31*std::pow(_ep1, 2)+2.*_w1*_de3*_ep1-_w1*std::pow(_de3, 2)+std::pow(fP3lab.Pt(), 2)*std::pow(_ep1, 2))/(fP3lab.E()*_p+fP3lab.Pz()*_ep1);
  
  b1 = (qve.Px()*sp5-qve.Py()*cp5)*fP5lab.Pt();
  b2 = (qve.Pz()*_ep2+qve.E()*_p)*fP5lab.Pt();
  b3 = (_w52*std::pow(_ep2, 2)+2.*_w2*_de5*_ep2-_w2*std::pow(_de5, 2)+std::pow(fP5lab.Pt()*_ep2, 2))/(_ep2*fP5lab.Pz()-fP5lab.E()*_p); //OK
  
  r12 =  c2*sp3+qve.Py()*c3;
  r13 = -c2*cp3-qve.Px()*c3;
  
  DebugInsideLoop(Form("qve = (%f, %f, %f, %f)", qve.E(), qve.Px(), qve.Py(), qve.Pz()));
  
  r22 =  b2*sp5+qve.Py()*b3;
  r23 = -b2*cp5-qve.Px()*b3;
  
  _epsi = _p12*c1*b1+r12*r22+r13*r23;

  _g5 = _w1*std::pow(c1, 2)+std::pow(r12, 2)+std::pow(r13, 2);
  _g6 = _w2*std::pow(b1, 2)+std::pow(r22, 2)+std::pow(r23, 2);

  _a5 = -(qve.Px()*cp3+qve.Py()*sp3)*fP3lab.Pt()*_p1k2-(_ep1*qve.E()-_p*qve.Pz())*(cp3*cp5+sp3*sp5)*fP3lab.Pt()*fP5lab.Pt()+(_de5*qve.Pz()+qve.E()*(_p+fP5lab.Pz()))*c3;
  _a6 = -(qve.Px()*cp5+qve.Py()*sp5)*fP5lab.Pt()*_p2k1-(_ep2*qve.E()+_p*qve.Pz())*(cp3*cp5+sp3*sp5)*fP3lab.Pt()*fP5lab.Pt()+(_de3*qve.Pz()-qve.E()*(_p-fP3lab.Pz()))*b3;
 
  DebugInsideLoop(Form("a5 = %f\n\ta6 = %f", _a5, _a6));
 
  //std::cout << _a5 << ">>>" << _a6 << std::endl;

  ////////////////////////////////////////////////////////////////
  // END of GAMGAMLL subroutine in the FORTRAN version
  ////////////////////////////////////////////////////////////////

  const Particle::Momentum cm = GetParticle(Particle::IncomingBeam1)->GetMomentum()
                               +GetParticle(Particle::IncomingBeam2)->GetMomentum();

  ////////////////////////////////////////////////////////////////
  // INFO from f.f
  ////////////////////////////////////////////////////////////////

  const double gamma = cm.E()/fSqS, betgam = cm.Pz()/fSqS;

  if (fCuts.mode==0) {
    Debug(Form("No cuts applied on the outgoing leptons kinematics!"));
  }
  // Kinematics computation for both leptons

  fP6cm.BetaGammaBoost(gamma, betgam);
  fP7cm.BetaGammaBoost(gamma, betgam);

  lcut = false; // Event discarded by default
  const double cott6 = fP6cm.Pz()/fP6cm.Pt(), cott7 = fP7cm.Pz()/fP7cm.Pt();

  // Cuts on outgoing leptons' kinematics

  lmu1 = cott6>=_cotth1
     and cott6<=_cotth2
     and (fP6cm.Pt()>=fCuts.ptmin or fCuts.ptmin<=0.)
     and (fP6cm.Pt()<=fCuts.ptmax or fCuts.ptmax<=0.)
     and (fP6cm.E()>=fCuts.emin   or fCuts.emin <=0.)
     and (fP6cm.E()<=fCuts.emax   or fCuts.emax <=0.);
  lmu2 = cott7>=_cotth1
     and cott7<=_cotth2
     and (fP7cm.Pt()>=fCuts.ptmin or fCuts.ptmin<=0.)
     and (fP7cm.Pt()<=fCuts.ptmax or fCuts.ptmax<=0.)
     and (fP7cm.E()>=fCuts.emin   or fCuts.emin <=0.)
     and (fP7cm.E()<=fCuts.emax   or fCuts.emax <=0.);

  switch (fCuts.mode) {
    case 0: default:
      lcut = true; break;
    case 1: // Vermaseren's hypothetical detector cuts
      {
        const double cost6 = fP6cm.Pz()/fP6cm.P(), cost7 = fP7cm.Pz()/fP7cm.P();
        lcut = ((fabs(cost6)<=0.75 and _pt_l6>=1.) or (fabs(cost6)<=0.95 and fabs(cost6)>0.75 and fabs(fP6cm.Pz())>1.)) and
               ((fabs(cost7)<=0.75 and _pt_l7>=1.) or (fabs(cost7)<=0.95 and fabs(cost7)>0.75 and fabs(fP7cm.Pz())>1.));
      }
      break;
    case 2: lcut = lmu1 and lmu2; break;
    case 3: lcut = lmu1  or lmu2; break;
  }

  // Cut on mass of final hadronic system (MX)
  if (fCuts.kinematics>1) {
    if (_mp3<fCuts.mxmin or _mp3>fCuts.mxmax) return 0.;
    if (fCuts.kinematics==4) {
      if (_mp5<fCuts.mxmin or _mp5>fCuts.mxmax) return 0.;
    }
  }

  // Cut on the proton's Q2 (first photon propagator T1)
  if ((fCuts.q2max!=-1. and _t1<-fCuts.q2max) or _t1>-fCuts.q2min) lcut = false;

  if (!lcut) { // Dismiss the cuts-failing events in the cross-section computation
    return 0.;
  }

  switch (fCuts.kinematics) { // FIXME inherited from CDF version
  default: case 0: weight = sconst*fJacobian*PeriPP(1, 2); break; // ep case
  case 1:          weight = sconst*fJacobian*PeriPP(2, 2); break; // elastic case
  case 2:  case 3: weight = sconst*fJacobian*PeriPP(3, 2)*std::pow(_dw31,2); break; // single-dissociative case
  case 4:          weight = sconst*fJacobian*PeriPP(3, 3)*std::pow(_dw31*_dw52,2); break; // double-dissociative case
  }

  return weight;
}

void
GamGamLL::FillKinematics(bool)
{
  const Particle::Momentum cm = GetParticle(Particle::IncomingBeam1)->GetMomentum()
                               +GetParticle(Particle::IncomingBeam2)->GetMomentum();
  const double gamma = cm.E()/fSqS, betgam = cm.Pz()/fSqS;

  double ranphi;
  int rany, ransign/*, ranz*/;
  
  // debugging variables
  double gmux, gmuy, gmuw, gmunu;
  
  // Needed to parametrise a random rotation around z-axis
  rany = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
  ransign = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
  ranphi = ((double)rand()/RAND_MAX)*2.*pi;
  /*ranz = 1;
  if (symmetrise_) {
    ranz = ((double)rand()>=.5*RAND_MAX) ? 1 : -1;
    //_pp3 *= ranz;
    //_pp5 *= ranz;
    }*/
  
  // First incoming proton
  Particle* ip1 = GetParticle(Particle::IncomingBeam1);
  Particle::Momentum plab_ip1(0., 0., _p, _ep1);
  plab_ip1.BetaGammaBoost(gamma, betgam);
  ip1->SetMomentum(plab_ip1);
  // Error("Invalid incoming proton 1");
  
  // Second incoming proton
  Particle* ip2 = GetParticle(Particle::IncomingBeam2);
  Particle::Momentum plab_ip2(0., 0., -_p, _ep2);
  plab_ip2.BetaGammaBoost(gamma, betgam);
  ip2->SetMomentum(plab_ip2);
  // Error("Invalid incoming proton 2");
  
  // First outgoing proton
  Particle* op1 = GetParticle(Particle::OutgoingBeam1);
  fP3lab.BetaGammaBoost(gamma, betgam);
  fP3lab.RotatePhi(ranphi, rany);
  op1->SetMomentum(fP3lab);
  // Error("Invalid outgoing proton 1");
  if (fCuts.kinematics>1) { // fragmenting remnants
    op1->status = Particle::Undecayed;
    op1->SetM(_mp3);
  }
  else { // stable proton
    op1->status = Particle::FinalState;
    op1->SetM(); //FIXME
  }
  
  // Second outgoing proton
  Particle* op2 = GetParticle(Particle::OutgoingBeam2);
  fP5lab.BetaGammaBoost(gamma, betgam);
  fP5lab.RotatePhi(ranphi, rany);
  op2->SetMomentum(fP5lab);
  // Error("Invalid outgoing proton 2");
  if (fCuts.kinematics==4) { // fragmenting remnants
    op2->status = Particle::Undecayed;
    op2->SetM(_mp5);
  }
  else { // stable proton
    op2->status = Particle::FinalState;
    op2->SetM(); //FIXME
  }

  // First incoming photon
  // Equivalent in LPAIR : PLAB(x, 3)
  Particle* ph1 = GetParticle(Particle::Parton1);
  Particle::Momentum plab_ph1 = plab_ip1-fP3lab;
  plab_ph1.RotatePhi(ranphi, rany);
  ph1->SetMomentum(plab_ph1);
  ////Error("Invalid photon 1");
  ph1->charge = 0;
  ph1->status = Particle::Incoming; // "incoming beam"
  
  // Second incoming photon
  // Equivalent in LPAIR : PLAB(x, 4)
  Particle* ph2 = GetParticle(Particle::Parton2);
  Particle::Momentum plab_ph2 = plab_ip2-fP5lab;
  plab_ph2.RotatePhi(ranphi, rany);
  ph2->SetMomentum(plab_ph2);
  ////Error("Invalid photon 2");
  ph2->charge = 0;
  ph2->status = Particle::Incoming; // "incoming beam"

  // Central (two-photon) system
  Particle* cs = GetParticle(Particle::CentralSystem);
  cs->status = Particle::Incoming;

  Particle::Role role_ol1, role_ol2;
  if (ransign<0) { role_ol1 = Particle::CentralParticle1; role_ol2 = Particle::CentralParticle2; }
  else           { role_ol1 = Particle::CentralParticle2; role_ol2 = Particle::CentralParticle1; }
  
  // First outgoing lepton
  Particle* ol1 = GetParticle(role_ol1);
  ol1->SetPDGId(ol1->GetPDGId(), ransign);
  fP6cm.RotatePhi(ranphi, rany);
  ol1->SetMomentum(fP6cm);
  // Error("Invalid outgoing lepton 1");
  ol1->status = Particle::FinalState;
  ol1->SetM(); //FIXME
  
  // Second outgoing lepton
  Particle* ol2 = GetParticle(role_ol2);
  ol2->SetPDGId(ol2->GetPDGId(), -ransign);
  fP7cm.RotatePhi(ranphi, rany);
  ol2->SetMomentum(fP7cm);
  ol2->status = Particle::FinalState;
  ol2->SetM(); //FIXME

  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) {
    gmux = -_t2/(_ep1*_eg2-_pp1*_p3_g2[2])/2.;
    gmuy = (_ep1*plab_ph2.E()-_pp1*plab_ph2.Pz())/(_ep2*plab_ph2.E()+_pp2*plab_ph2.Pz());
    gmuw = std::pow(_ep1+plab_ph2.E(), 2)-std::pow(_pp1+plab_ph2.Pz(), 2);
    if (gmuw>=0.) {
      gmuw = std::sqrt(gmuw);
    }
    else {
      Warning(Form("W^2 = %f < 0", gmuw));
      gmuw = 0.;
    }
    gmunu = gmuy*2.*Particle::GetMassFromPDGId(Particle::Proton)/_ep1/_ep2;
    DebugInsideLoop(Form(" gmux = %f\n\t"
                         " gmux = %f\n\t"
                         " gmuw = %f\n\t"
                         "gmunu = %f", gmux, gmuy, gmuw, gmunu));
  }
  //fEvent->Dump();
}

double
GamGamLL::PeriPP(int nup_, int ndown_)
{
  DebugInsideLoop(Form(" Nup  = %d\n\tNdown = %d", nup_, ndown_));

  FormFactors fp1, fp2;

  switch(nup_) {
    case 1:  fp1 = TrivialFormFactors(); break; // electron (trivial) form factor
    case 2:  fp1 = ElasticFormFactors(-_t1, _w1); break; // proton elastic form factor
    case 4: { // does not exist in CDF version
      double dummy, psfw1, psfw2;
      PSF(_t1, _w3, &dummy, &psfw1, &psfw2);
      fp1.FM = -psfw1*(2.*_mp1)/_t1;
      fp1.FE =  psfw2/(2.*_mp1);
      break;
    }
    default: fp1 = SuriYennieFormFactors(-_t1, _w1, _w3); break;
  }

  switch(ndown_) {
    case 1:  fp2 = TrivialFormFactors(); break; // electron (trivial) form factor
    case 2:  fp2 = ElasticFormFactors(-_t2, _w2); break; // proton elastic form factor
    default: fp2 = SuriYennieFormFactors(-_t2, _w2, _w5); break;
  }
  
  DebugInsideLoop(Form("u1 = %f\n\tu2 = %f\n\tv1 = %f\n\tv2 = %f", fp1.FM, fp1.FE, fp2.FM, fp2.FE));

  const double qqq = std::pow(_q1dq, 2),
               qdq = 4.*_w6-_w4;
  const double t11 = 64. *( _bb*(qqq-_g4-qdq*(_t1+_t2+2.*_w6))-2.*(_t1+2.*_w6)*(_t2+2.*_w6)*qqq)*_t1*_t2,
               t12 = 128.*(-_bb*(_dd2+_g6)-2.*(_t1+2.*_w6)*(_sa2*qqq+std::pow(_a6, 2)))*_t1,
               t21 = 128.*(-_bb*(_dd4+_g5)-2.*(_t2+2.*_w6)*(_sa1*qqq+std::pow(_a5, 2)))*_t2,
               t22 = 512.*( _bb*(std::pow(_delta, 2)-_gram)-std::pow(_epsi-_delta*(qdq+_q1dq2), 2)-_sa1*std::pow(_a6, 2)-_sa2*std::pow(_a5, 2)-_sa1*_sa2*qqq);

  const double peripp = ( fp1.FM*fp2.FM*t11
                         +fp1.FE*fp2.FM*t21
                         +fp1.FM*fp2.FE*t12
                         +fp1.FE*fp2.FE*t22) / pow(2.*_t1*_t2*_bb, 2);

  DebugInsideLoop(Form("t11 = %5.2f\tt12 = %5.2f\n\t"
                       "t21 = %5.2f\tt22 = %5.2f\n\t"
                       "=> PeriPP = %f", t11, t12, t21, t22, peripp));
  
  return peripp;
}
