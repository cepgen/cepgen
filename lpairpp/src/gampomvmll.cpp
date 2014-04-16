#include "gampomvmll.h"

GamPomVMLL::GamPomVMLL():
  _name("gamma,pomeron->VM->l+,l-"),
  _epsilm(0.0808),
  _wmin(20.), _wmax(0.),
  _q2min(4.), _q2max(100.),
  _ymin(0.), _ymax(1.),
  _b0(4.), _wb0(95.), _amxb0(14.), _anexp(0.),
  _wsig0(95.),
  _q2(0.),
  _gengam_first(true),
  _gephot_first(true),
  _fraggl_begin(true)
{
  pe = 100.;
  dme = 0.000511; // electron mass
  pp = 100.;
  dmp = 0.9; // proton mass

  itypvm = 553;
  ifragp = 1;
  ifragv = 1;

  _s = 4.*pe*pp;
  _ecm = std::sqrt(_s);

  _wmax = std::sqrt(_s+std::pow(dme, 2)+std::pow(dmp, 2));

  this->GenGam();

}

GamPomVMLL::~GamPomVMLL()
{
}

void
GamPomVMLL::GenGam()
{
  int igen, igent, igenl, iacc, iacct, iaccl;
  double smax, egammin;

  if (_gengam_first) {
    _gengam_first = false;
    igen = igent = igenl = iacc = iacct = iaccl = 0;

    smax = std::pow(_wmax, 2);
    egammin = std::pow(_wmin, 2)/4./pp;

    _w2 = std::pow(_wsig0, 2);
  }
}

double
GamPomVMLL::ComputeWeight()
{
  double w, wght;
  double yhat;
  bool look = false;
  double genmxt, mxt;

  // From common block /CGDIF/:                                                                                                                     
  //   IFRAGP, IFRAGV, DEMINP, AMASSV, B0, AMXB0, WB0,                                                                                                
  //   ALPHA1, ALPH1M, WMIN, PE, PP                                                                                                                   
  // From common block /CMASS/:                                                                                                                     
  //   DMP, DMPI0, DMVM, DWVM, DMNST, DWNST
  int ifragp, ifragv;
  double deminp, amassv;
  double alph1m;
  double dmp, dmpio, dmnst, dwnst;
  //
  double t, tmin, tmax, tmean;
  double pcm1;
  double b, bmin;
  double dmmin, dmxp, dmxv;

  if (_genmxt_begin) {
    _genmxt_begin = false;
    if (ifragp!=1 and ifragp!=-1 and ifragp!=2 and ifragv==0) {
      bmin = _b0+4.*_alpha1*log(_wmin/_wb0);
    }
    else if ((ifragp==1 or ifragp==-1 or ifragp==2) and ifragv!=0) {
      bmin = _b0+4.*_alpha1*log(4.*std::pow(_amxb0, 2)/(_wb0*_ecm));
    }
    else {
      bmin = _b0+4.*_alpha1*log(_amxb0/_wb0);
    }
    bmin = std::max(bmin, 0.5);
  }

  w = std::sqrt(_w2);
  
  // Generate masses at p and VM vertex
  
  if (ifragp==0) {
    dmxp = dmp;
  }
  else if (ifragp==1 or ifragp==-1 or ifragp==2) {
    dmxp = PXMass(dmp+deminp, _ecm);
  }
  else {
    dmxp = RanBW(dmnst, dwnst, dmp+deminp, dmnst+2.*dwnst);
  }
  
  if (ifragv!=0) {
    dmxv = VXMass(amassv, _ecm);
  }
  else {
    dmmin = _dmvm-3.*_dmvm;
    if (itypvm==100113 or itypvm==30113) {
      dmmin = std::max(dmmin, 1.2);
      }
    else if (itypvm==10333) {
      dmmin = std::max(dmmin, 1.4);
    }
    dmxv = RanBW(_dmvm, _dmvm, dmmin, _dmvm+10.*_dmvm);
  }
  
  // Return if generated masses are bigger than CM energy
  
  if (dmxp+dmxv>w-0.1) {
    t = b = yhat = _pcm3 = genmxt = 0.;
    wght = 1.;
      if (look) {
	//CALL SHS (2, 0, LOG10 (SNGL (DMXP)))
	//CALL SHS (3, 0, LOG10 (SNGL (DMXV)))
      }
      return wght;
  }
    
  // Calculate slope parameter b
  // Generate t with e**(b*t) distribution
  
  b = _b0+4.*_alpha1*log(w/_wb0);
  if (ifragp==1 or ifragp==-1 or ifragp==2) {
    b -= 4.*alph1m*log(dmxp/_amxb0);
  }
  if (ifragv!=0) {
    b -= 4.*_alpha1*log(dmxv/_amxb0);
  }
  
  if (b<.5) b = .5;
  //CALL GENERT (T, 0.0D0, S, B, 1D0*ANEXP)
  
  // Calculate actual minimal and maximal t for the generated masses
  // Note that t here is positive!
  // Formula (E.5) from Review of Particle Properties 1992, p. III.50
  // 1: gamma, 2: p, 3: VM(+X), 4: p remnant
  // The formula for Pcm1 is altered to take the imaginary photon mass
  // into account.
  
  pcm1 = std::sqrt(std::pow(_w2+_q2-std::pow(dmp, 2), 2)+4.*_q2*std::pow(dmp, 2))/w/2.;
  _pcm3 = std::sqrt((_w2-std::pow(dmxv+dmxp, 2))*(_w2-std::pow(dmxv-dmxp, 2)))/w/2.;
  tmean = ((-_q2-std::pow(dmp, 2))*(std::pow(dmxv, 2)-std::pow(dmxp, 2))/_w2+_w2+_q2-std::pow(dmp, 2)-std::pow(dmxv, 2)-std::pow(dmxp, 2))/2.;
  tmin = tmean-2.*pcm1*_pcm3;
  tmax = tmean+2.*pcm1*_pcm3;
  
  if (t<=tmax and t>=tmin) {
    mxt = 1;
    yhat = (t-tmin)/(4*pcm1*_pcm3);
  }
  else {
    mxt = 0;
  }
  
  if (look) {
    //CALL SHS (1, 0, SNGL (T))
    //CALL SHS (2, 0, LOG10 (SNGL (DMXP)))
    //CALL SHS (3, 0, LOG10 (SNGL (DMXV)))
    if (mxt==1) {
      //CALL SHS (1, 1, SNGL (T))
      //CALL SHS (2, 1, LOG10 (SNGL (DMXP)))
      //CALL SHS (3, 1, LOG10 (SNGL (DMXV)))
    }
  }

  wght = bmin/b;
  genmxt = mxt*wght;
  
  return genmxt;

}

double
GamPomVMLL::PXMass(double mmin_, double mmax_)
{
  int iter;
  double lmin, delta;
  double m2;
  double m2min, fact;
  double mmin2, mmax2;
  double y;

  mmin2 = std::pow(mmin_, 2);
  mmax2 = std::pow(mmax_, 2);

  if (fabs(_epsilm)<.001) {
    lmin = 2.*log(mmin_);
    delta = 2.*log(mmax_/mmin_);
  }
  else {
    m2min = std::pow(mmin_, -2*_epsilm);
    fact = std::pow(mmax_, -2*_epsilm)-m2min;
  }

  iter = 0;
  do {
    iter++;
    
    if (fabs(_epsilm)<.001) { // Basic spectrum: 1/M^2
      m2 = exp(random()*delta+lmin);
    }
    else { // Basic spectrum: 1/M^2(1+epsilon)
      m2 = std::pow(fact*random()+m2min, -1./_epsilm);
    }
    if (m2<mmin2) {
      std::cerr << "[GamPomVMLL::PXMass] ERROR: M2 = " << m2 << " < MMIN**2 = " << mmin2 << std::endl;
      //CALL ERRLOG (100, 'S: PXMASS: M2 < MMIN**2!')
      m2 = mmin2;
    }
    else if (m2>mmax2) {
      std::cerr << "[GamPomVMLL::PXMass] ERROR: M2 = " << m2 << " > MMAX**2 = " << mmax2 << std::endl;
      //CALL ERRLOG (101, 'S: PXMASS: M2 > MMAX**2!')
      m2 = mmax2;
    }

    // Old version with enhancements in lower mass region
    if      (m2>=4.)   y = 1.;
    else if (m2>=3.1)  y = 1.64-0.16*m2;
    else if (m2>=2.65) y = m2*(0.47-0.42*std::pow(m2-2.65, 2));
    else if (m2>=2.25) y = m2*(0.47+0.46*std::pow(m2-2.65, 2));
    else if (m2>=2.02) y = m2*(0.76-2.69*std::pow(m2-2.02, 2));
    else if (m2>=1.72) y = m2*(0.76-1.98*std::pow(m2-2.02, 2));
    else               y = 1.05*(m2-1.165);

  } while (1.6*random()>y and iter<=100);

  if (y>1.6) {
    std::cout << "[GamPomVMLL::PXMass] WARNING: Y = " << y << " for M2 = " << m2 << std::endl;
    //CALL ERRLOG (102, 'W: PXMASS: Y > 1.6')
  }

  if (iter>100) {
    std::cout << "[GamPomVMLL::PXMass] WARNING: more than 100 iterations!" << std::endl;
    //CALL ERRLOG (103, 'W: PXMASS: More than 100 iterations!')
  }

  return std::sqrt(m2);
}

double
GamPomVMLL::VXMass(double mmin_, double mmax_)
{
  double m2;
  double lmin, delta;
  double m2min, fact;
  double mmin2, mmax2;

  mmin2 = std::pow(mmin_, 2);
  mmax2 = std::pow(mmax_, 2);

  if (fabs(_epsilm)<.001) {
    lmin = 2.*log(mmin_);
    delta = 2.*log(mmax_/mmin_);
  }
  else {
    m2min = std::pow(mmin_, -2*_epsilm);
    fact = std::pow(mmax_, -2*_epsilm)-m2min;
  }

  if (fabs(_epsilm)<.001) { // Basic spectrum: 1/M^2
    m2 = exp(random()*delta+lmin);
  }
  else { // Basic spectrum: 1/M^2(1+epsilon)
    m2 = std::pow(fact*random()+m2min, -1./_epsilm);
  }
  if (m2<mmin2) {
    std::cerr << "[GamPomVMLL::VXMass] ERROR: M2 = " << m2 << " < MMIN**2 = " << mmin2 << std::endl;
    //CALL ERRLOG (100, 'S: VXMASS: M2 < MMIN**2!')
    m2 = mmin2;
  }
  else if (m2>mmax2) {
    std::cerr << "[GamPomVMLL::VXMass] ERROR: M2 = " << m2 << " > MMAX**2 = " << mmax2 << std::endl;
    //CALL ERRLOG (101, 'S: VXMASS: M2 > MMAX**2!')
    m2 = mmax2;
  }

  return std::sqrt(m2);

}

void
GamPomVMLL::FragGl()
{
  int ivvm, ipom, idifv, ivm; //FIXME : arguments

  double glumas, gluwid, dmasvm, dmasgl, dmass;
  double pcmgam[4], pvmvm[4], pcmglu[4], pt[3];
  double dmu1, dmu2, dmu3, dmu4;
  double c1, c2, c3;
  double t, tmin, tmax, yhat;
  double ctheta, stheta, phi, pin, pout, pgamf, ptf, b;
  int i, iglue;

  int idahep[10][2], mohep[10][2];
  int istat[10], itype[10];
  int npart;

  if (_fraggl_begin) {
    _fraggl_begin = false;
    glumas = GetMassFromPDGId(ifragv);
    gluwid = glumas/10.;
  }

  dmass = _ppcms8[idifv][4];

  if (dmass<_dmvm+glumas) {
    std::cerr << "[GamPomVMLL::FragGl] ERROR: not enough energy!" << std::endl;
    //CALL ERRLOG (90, 'F: FRAGGL: Not enough energy!')
    exit(0);
  }

  // Choose the actual VM and glueball masses
  do {
    dmasvm = RanBW(_dmvm, _dmvm, _dmvm-2.*_dwvm, _dmvm+2.*_dmvm);
    dmasgl = RanBW(glumas, gluwid, glumas-2.*gluwid, glumas+2.*gluwid);
  } while (dmasvm+dmasgl>=dmass);

  // Choose momentum transfer t
  // assume that b = 4GeV^-2 at a mass of 10GeV
  b = 4.+4.*_alpha1*log(dmass/10.);

  // Calculate actual minimal and maximal t for the generated masses
  // Note that t here is positive!
  // Formulae from PYTHIA manual (CERN-TH.7112/93), p. 99.
  // 1: virtual VM, 2: virtual pomeron, 3: real VVM, 4: glueball
  // (tmin and tmax are defined differently than for PYTHIA, since signs
  // are reversed)

  dmu1 = -std::pow(_ppcms8[ivvm][4]/dmass, 2);
  dmu2 = -std::pow(_ppcms8[ipom][4]/dmass, 2);
  dmu3 = std::pow(dmasvm/dmass, 2);
  dmu4 = std::pow(dmasgl/dmass, 2);

  c1 = 1.-(dmu1+dmu2+dmu3+dmu4)+(dmu1-dmu2)*(dmu3-dmu4);
  c2 = std::sqrt((std::pow(1.-dmu1-dmu2, 2)-4.*dmu1*dmu2)*(std::pow(1.-dmu3-dmu4, 2)-4.*dmu3*dmu4));
  c3 = (dmu3-dmu1)*(dmu4-dmu2)+(dmu1+dmu4-dmu2-dmu3)*(dmu1*dmu4-dmu2*dmu3);

  tmax = std::pow(_ppcms8[idifv][4], 2)*(c1+c2)/2.;
  tmin = std::pow(_ppcms8[idifv][4], 4)*c3/tmax;

  //CALL GENERT (T, TMIN, TMAX, B, 1D0*ANEXP)

  pin = dmass*std::sqrt(std::pow(1.-dmu1-dmu2, 2)-4.*dmu1*dmu2)/2.;
  pout = dmass*std::sqrt(std::pow(1.-dmu3-dmu4, 2)-4.*dmu3*dmu4)/2.;

  yhat = (t-tmin)/(4.*pin*pout);
  ctheta = 1.-2.*yhat;
  stheta = 2.*sqrt(yhat-std::pow(yhat, 2));

  // Calculate the 5-vectors of the VM and glueball in the gamma-pomeron CMS

  //CALL LORENF8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PPCMS8 (1, IVVM), PCMGAM)

  pgamf = pout*ctheta/std::sqrt(std::pow(pcmgam[0], 2)+std::pow(pcmgam[1], 2)+std::pow(pcmgam[2], 2));

  phi = 2.*pi*random();
  pt[0] = -cos(phi)*pcmgam[2];
  pt[1] =  sin(phi)*pcmgam[2];
  pt[2] = -sin(phi)*pcmgam[1]+cos(phi)*pcmgam[0];
  ptf = pout*stheta/std::sqrt(std::pow(pcmgam[2], 2)+std::pow(pt[2], 2));

  _pcmvm[3] = std::pow(_dmvm, 2);
  for (i=0; i<3; i++) {
    _pcmvm[i] = pgamf*pcmgam[i]+ptf*pt[i];
    _pcmvm[3]+= std::pow(_pcmvm[i], 2);
  }
  _pcmvm[3] = std::sqrt(_pcmvm[3]);

  if (fabs(std::pow(pout, 2)-std::pow(_pcmvm[0], 2)-std::pow(_pcmvm[1], 2)-std::pow(_pcmvm[2], 2))>std::pow(pout, 2)/100.) {
    std::cerr << "[GamPomVMLL::FragGl] WARNING: POUT <> |PCMVM|" << std::endl;
    //CALL ERRLOG (91, 'W: FRAGGL: POUT <> |PCMVM|')
  }

  pcmglu[3] = std::pow(dmasgl, 2);
  for (i=0; i<3; i++) {
    pcmglu[i] = -_pcmvm[i];
    pcmglu[3]+= std::pow(pcmglu[i], 2);
  }
  pcmglu[3] = std::sqrt(pcmglu[3]);

  npart = 0; //FIXME
 
  idahep[idifv][0] = npart+1;
  idahep[idifv][1] = npart+2;
  istat[idifv] = 2;

  // Glueball quantities
  iglue = npart+1;
  itype[iglue] = ifragv;
  //CALL LORENB8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PCMGLU, PPCMS8 (1, IGLUE))
  _ppcms8[iglue][4] = dmasgl;
  istat[iglue] = 1;
  idahep[iglue][0] = 0;
  idahep[iglue][1] = 0;
  mohep[iglue][0] = idifv;
  mohep[iglue][1] = 0;

  // Vector meson quantities
  ivm = npart+2;
  itype[ivm] = itypvm;
  //CALL LORENB8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PCMVM, PPCMS8 (1, IVM))
  _ppcms8[ivm][4] = dmasvm;
  istat[ivm] = 1;
  idahep[ivm][0] = 0;
  idahep[ivm][1] = 0;
  mohep[ivm][0] = idifv;
  mohep[ivm][1] = 0;

  npart += 2;

  // Perform glueball decay 
  //CALL DECGLU (IGLUE)

  // Glueball can decay to K* K*bar => call DECK0
  //CALL DECK0 (IGLUE+1, NPART)

}

void
GamPomVMLL::GEPhot(int igammd_)
{
  int iacc, isum, iacct, iaccl;
  double dsum, qsum, dsumt, qsumt, dsuml, qsuml;
  double eellab, elpr, s, esmp2, eel, wmin2, w12, ysqr;
  double dymax, dymin;
  double gq2min, gq2max;


  if (_gephot_first) {
    iacc = isum = iacct = iaccl = 0;
    dsum = qsum = dsumt = qsumt = dsuml = qsuml = 0.;

    eellab = std::sqrt(std::pow(pe, 2)+std::pow(dme, 2));
    elpr = std::sqrt(std::pow(pp, 2)+std::pow(dmp, 2))*eellab+pp*pe;
    s = 2.*elpr+std::pow(dme, 2)+std::pow(dmp, 2);
    esmp2 = std::pow(2.*elpr+std::pow(dme, 2), 2);

    if (igammd_>3) {
      eel = elpr/dmp;
    }
    else eel = eellab;

    wmin2 = std::pow(_wmin, 2);

    _gephot_first = false;
  }
}

bool GamPomVMLL::SetIncomingParticles(Particle, Particle){}
bool GamPomVMLL::SetOutgoingParticles(int, int){}
void GamPomVMLL::FillKinematics(bool){}
void GamPomVMLL::SetKinematics(Kinematics){}
void GamPomVMLL::ComputeCMenergy(){}
double GamPomVMLL::ComputeMX(double x_, double outmass_, double* dw_){}
void GamPomVMLL::StoreEvent(std::ofstream*,double){}
void GamPomVMLL::PrepareHadronisation(Particle *part_){}
