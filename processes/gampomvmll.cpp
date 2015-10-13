#include "gampomvmll.h"

GamPomVMLL::GamPomVMLL():
  // Parameters from steer.text
  _cthelb(-.9962), _eelmin(8.),
  _lambda(0.), _eprop(2.5), _xi(1.), _chi(1.),
  _epsilw(0.225), _epsilm(0.0808), _alpha1(0.), _alph1m(_alpha1),
  _igammd(1), _egamma(3.), _wmin(20.), _wmax(0.),
  _q2min(4.), _q2max(100.),
  _ymin(0.), _ymax(1.),
  _b0(4.), _wb0(95.), _amxb0(14.), _anexp(0.),
  _wsig0(95.),

  // Parameters from GDIINI
  /*_cthelb(-1.), _eelmin(0.),
  _lambda(0.), _eprop(2.), _xi(1.), _chi(0.),
  _epsilw(0.0808), _epsilm(0.0808), _alpha1(0.25), _alph1m(-1.),
  _igammd(1), _egamma(3.), _wmin(0.), _wmax(0.),
  _q2min(0.), _q2max(0.),
  _ymin(0.), _ymax(1.),
  _b0(10.), _wb0(14.), _amxb0(0.), _anexp(0.),
  _wsig0(14.),*/

  _gengam_w2(0.),

  _q2(0.),

  _genmxt_begin(true),
  _gengam_first(true),
  _gephot_first(true),
  _fraggl_begin(true),
  _photint_swei(0.), _photint_swei2(0.), _photint_sweit(0.), _photint_sweit2(0.), _photint_sweil(0.), _photint_sweil2(0.),
  _vmflux_f(0.), _vmflux_df(0.), _vmflux_fl(0.), _vmflux_dfl(0.), _vmflux_ft(0.), _vmflux_dft(0.)
{
  _name = "gamma,pomeron->VM->l+,l-";

  itypvm = Particle::UPSILON_1S;
  ifragp = 0;
  deminp = 0.236;
  ifragv = (Particle::ParticleCode)0;
  amassv = 0.;
  idifv = 1;
  ivvm = 2;
  ipom = 3;
  ivm = 4;

  _br = GetBRFromProcessId(itypvm); //FIXME what about other final states ???
  
}

GamPomVMLL::~GamPomVMLL()
{
}

void
GamPomVMLL::GDIBeg()
{
  double r;
  double wminmin;

  _dme = Particle::GetMassFromPDGId(Particle::ELECTRON);
  _dmp = Particle::GetMassFromPDGId(Particle::PROTON);
  _dmpi = Particle::GetMassFromPDGId(Particle::PI_PLUS);
  _dmpi0 = Particle::GetMassFromPDGId(Particle::PI_0);
  _dmn = Particle::GetMassFromPDGId(Particle::NEUTRON);
  _dmvm = Particle::GetMassFromPDGId((Particle::ParticleCode)itypvm);
  _dwvm = Particle::GetWidthFromPDGId((Particle::ParticleCode)itypvm); //FIXME
  _dml = _ev->GetOneByRole(OL1)->M();

  // For elastic N* production at p vertex initialize DMNST, DWNST
  if (abs(ifragp)>2) {
    _dmnst = Particle::GetMassFromPDGId(static_cast<Particle::ParticleCode>(ifragp));
    _dmnst = Particle::GetWidthFromPDGId(static_cast<Particle::ParticleCode>(ifragp));
    if (_dmnst<=0.) {
      //WRITE (ERTEXT, '(''F: GDIBEG: Mass of'',I8,''not known!'')'), IFRAGP
      //CALL ERRLOG (10, ERTEXT)
      std::cerr << "[GamPomVMLL::GDIBeg] ERROR: Mass of " << ifragp << " not known!" << std::endl;
      exit(0);
    }
  }

  // Check that beam particle is proton or antiproton
  if (abs(_ev->GetOneByRole(1)->pdgId)!=Particle::PROTON and abs(_ev->GetOneByRole(2)->pdgId)!=Particle::PROTON) {
    //WRITE (ERTEXT, '(''F: GDIBEG: Beam proton must be proton or antiproton. IBEAMP ='',I8)'), IBEAMP
    //CALL ERRLOG (11, ERTEXT)
    std::cerr << "[GamPomVMLL::GDIBeg] ERROR: Beam proton must be proton or antiproton. IBEAMP = " << _ev->GetOneByRole(1)->pdgId << " / " << _ev->GetOneByRole(2)->pdgId << std::endl;
    exit(0);
  }

  // If necessary, initialize LAMBDA
  if (_lambda<=0.) {
    if (itypvm==22) _lambda = Particle::GetMassFromPDGId(Particle::RHO_770_0);
    else _lambda = _dmvm;
    //std::cout << "[GamPomVMLL::GDIBeg] INFO: LAMBDA set to " << _lambda << std::endl;
  }

  // If necessary, initialize DEMINP, AMASSV
  if (deminp<_dmn+_dmpi0-_dmp) {
    deminp = _dmn+_dmpi0-_dmp+0.1;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: DEMINP set to " << deminp << std::endl;
  }
  if (abs(ifragp)>2 and deminp<_dmnst-2.*_dwnst-_dmp) {
    deminp = _dmnst-2.*_dwnst-_dmp;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: DEMINP set to " << deminp << std::endl;
  }

  if (amassv<2.*_dmpi) {
    if      ((itypvm/10)%10<2)  amassv = 1.0; // rho, omega states: mimimum is pi+ pi- state, stay away from rho
    else if ((itypvm/10)%10==3) amassv = 1.5; // phi states: minimum is K+ K- state
    else if ((itypvm/10)%10==4) amassv = 4.0; // psi states: minimum is J/psi pi+ pi- state
    else if ((itypvm/10)%10==5) amassv = 10.; // Upsilon states: minimum is Upsilon (1S) pi+ pi- state
    else {
      //CALL ERRLOG (12, 'F: GDIBEG: Unknown quark content of vector meson')
      std::cerr << "[GamPomVMLL::GDIBeg] ERROR: Unknown quark content of vector meson" << std::endl;
      exit(0);
    }
    std::cout << "[GamPomVMLL::GDIBeg] INFO: AMASSV set to " << amassv << std::endl;    
  }

  if ((int)ifragv>100) {
    r = Particle::GetMassFromPDGId(ifragv)+_dmvm;
    if (amassv<r) {
      amassv = r+1.;
      std::cout << "[GamPomVMLL::GDIBeg] INFO: AMASSV set to " << amassv << std::endl;
    }
  }

  // If necessary, initialize WMIN
  if (ifragp==0) wminmin = _dmp;
  else wminmin = _dmp+deminp;

  if (ifragv==(int)0) wminmin += _dmvm;
  else wminmin += amassv;

  if (_wmin<wminmin) {
    _wmin = wminmin+1.;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: WMIN set to " << _wmin << std::endl;
  }

  // If necessary, initialize WMAX
  if (_wmax<=_wmin) {
    _wmax = std::sqrt(4.*_pz1*_pz2+_ev->GetOneByRole(1)->M2()+_ev->GetOneByRole(2)->M2());
    std::cout << "[GamPomVMLL::GDIBeg] INFO: WMAX set to " << _wmax << std::endl;
  }

  // If necessary, initialize Q2MIN
  if (_eelmin>0. and _cthelb>-1. and _q2min<2.*_pz1*_eelmin*(1.+_cthelb)) {
    _q2min = 2.*_pz1*_eelmin*(1.+_cthelb);
    std::cout << "[GamPomVMLL::GDIBeg] INFO: Q2MIN set to " << _q2min << std::endl;
  }

  // If necessary, initialize Q2MAX
  if (_q2max<0.) {
    _q2max = fabs(_q2max);
    std::cout << "[GamPomVMLL::GDIBeg] INFO: Q2MAX set to " << _q2max << std::endl;
  }
  if (_q2max<=_q2min) {
    _q2max = 4.*_pz1*_pz2+_ev->GetOneByRole(IBE)->M2()+_ev->GetOneByRole(IBP)->M2();
    std::cout << "[GamPomVMLL::GDIBeg] INFO: Q2MAX set to " << _q2max << std::endl;
  }

  // If necessary, initialize AMXB0
  if (_amxb0<=0.) {
    if (ifragp==1 or ifragp==-1 or ifragp==2) {
      if (ifragv==(int)0) _amxb0 = _dmp;
      else _amxb0 = std::sqrt(_dmp+_dmvm);
    }
    else _amxb0 = _dmvm;
    //std::cout << "[GamPomVMLL::GDIBeg] INFO: AMXB0 set to " << _amxb0 << std::endl;
  }

  // If necessary, initialize BR
  if (_br==0.) {
    _br = 1.;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: BR set to " << _br << std::endl;
    //CALL ERRLOG (13, 'W: GDIBEG: BR was 0.0. Set to 1.0')
  }
  else if (_br>1.) {
    _br /= 100.;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: BR was > 1.0. Scaled down by 100 to " << _br << std::endl;
    //CALL ERRLOG (14, 'W: GDIBEG: BR was > 0.0. Scaled down by 100.')
  }

  // If necessary, initialize ALPH1M
  if (_alph1m<0.) {
    _alph1m = _alpha1;
    std::cout << "[GamPomVMLL::GDIBeg] INFO: ALPH1M set to " << _alph1m << std::endl;
    //CALL ERRLOG (15, 'I: GDIBEG: ALPH1M was set to ALPHA1')
  }
}

void
GamPomVMLL::GDIEvt()
{
  this->GenEvtDi();
}

void
GamPomVMLL::GenEvtDi()
{
#ifdef DEBUG
  std::cout << "[GamPomVMLL::GenEvtDi] [DEBUG] Generating the event" << std::endl;
#endif

  // GenBPr, GenBEl are done implicitely once the incoming particles are set for the process
  this->GenGam();
  
  // The following routines work in the gamma-p CM system and fill PPCMS8
  // Therefore: boost all 5-vectors into gamma-p CMS
  // for ()
  
  // Generate diffractive state
  this->GenDif();
  
}

void
GamPomVMLL::GenGam()
{
  int iacc, iter;
  double r;
  double sw, sw2, sw2bar;
  double sigwt, dsigwt;
  double wght;
  double wt;

  const unsigned int n = 10000;

  if (_gengam_first) {
    //_s = 4.*_ev->GetOneByRole(1)->pz*_ev->GetOneByRole(2)->pz; //FIXME FIXME FIXME
    _pz1 = fabs(_ev->GetOneByRole(1)->Pz()); //FIXME absolute value ???
    _pz2 = fabs(_ev->GetOneByRole(2)->Pz());
    _e1 = _ev->GetOneByRole(1)->E();
    _e2 = _ev->GetOneByRole(2)->E();
    _s = _ev->GetOneByRole(1)->M2()+_ev->GetOneByRole(2)->M2()+2*_e1*_e2-2*_pz1*_pz2;
    _ecm = std::sqrt(_s);
    _wmax = std::sqrt(_s+_ev->GetOneByRole(1)->M2()+_ev->GetOneByRole(2)->M2());

    this->GDIBeg();

    _gengam_first = false;


    _igen = _igent = _igenl = iacc = _iacct = _iaccl = 0;

    _event_smax = std::pow(_wmax, 2);
    _event_egammin = std::pow(_wmin, 2)/4./_ev->GetOneByRole(2)->Pz();

    _gengam_w2 = std::pow(_wsig0, 2);
    
    sw = sw2 = sw2bar = 0.;

    for (unsigned int i=0; i<n; i++) {
      r = this->GenMXT(&wt);
      //std::cout << "i = " << i << ", r = " << r << std::endl;
      sw += r;
      sw2 += std::pow(r, 2);
      sw2bar += std::pow(wt-r, 2);
    }

    if (sw<=0.) {
      std::cerr << "[GamPomVMLL::GenGam] ERROR: SW = " << sw << std::endl
		<< "  Cross section calculation impossible!" << std::endl;
      //CALL ERRLOG (60, 'F: GENGAM: SW < 0')
      exit(0);
    }

    std::cout << "[GamPomVMLL::GenGam] INFO: t/mx-combinations generated: " << n << std::endl
	      << "  Weight of t/mx-combinations accepted: " << sw << " (sw2 = " << sw2 << ", sw2bar = " << sw2bar << ")" << std::endl;
    
    _event_propmx = std::max(1., _xi*_q2min/(std::pow(_lambda, 2)+_xi*_chi*_q2min))/std::pow(1.+_q2min/std::pow(_lambda, 2), _eprop);
    sigwt = std::pow(_gengam_w2/_event_smax, (2.*_epsilw))/_event_propmx*sw*n;
    sw2bar = std::max(sw2bar, 1.);
    dsigwt = sigwt*std::sqrt(sw2*sw2bar/n)/sw;

    std::cout << "  sigwt = " << sigwt << ", dsigwt = " << dsigwt << ", propmx = " << _event_propmx << std::endl;

    // Print gamma-n cross section for W values in WVAL (skipped!)

    _gengam_first = false;
  }
  /*  
  iter = 0;
  do {
    iter++;
    wght = OneEvent();
  } while (wght<drand() and wght>=0.);
  
  //CALL SHSW (10, 3, SNGL (PCM (5)), 1.0)
  
  if (_event_heli==0) _iaccl++;
  else _iacct++;
  iacc++;
  */
}

double
GamPomVMLL::OneEvent()
{
  double weight;
  double pcm[4];
  double drlt;
  double wt;

  // Generate photons
  Particle pgam(41, Particle::PHOTON), pesc(5, _ev->GetOneByRole(IBE)->pdgId);
  //std::cout << "-> " << _igammd << std::endl;

  if (_igammd<0) { // Fixed photon energy
    //CALL FIXPHOT (PGAM, PESC, Q2, PPART8 (1, IBE), DBLE (EGAMMA))
    FixPhot(&pgam, &pesc, &_q2, *(_ev->GetOneByRole(1)), _egamma);
    _event_heli = Heli(0.);
  }
  else if (_igammd==0) { // Simple 1/k spectrum
    //CALL GENPHOT (PGAM, PESC, Q2, PPART8 (1, IBE), EGAMMIN, 0.0D0)
    GenPhot(&pgam, &pesc, &_q2, *(_ev->GetOneByRole(1)), _event_egammin, 0.);
    _event_heli = Heli(0.);
  }
  else {
    // 1 -> WWA
    // 2 -> full transverse
    // 3 -> full transverse and longitudinal spectrum
    // 4 -> full transverse and longitudinal spectrum in p rest frame
    GEPhot(&_q2, &_event_heli);
    //std::cout << "photon helicity = " << _event_heli << std::endl;
    //CALL GEPHOT (PPART8 (1, IBE), PPART8 (1, IBP), PGAM, PESC, Q2, HELI)
  }
  //std::cout << "q2=" << _q2 << std::endl;
  
  if (_event_heli==0) _igenl++;
  else _igent++;
  _igen++;
  
  // Determine actual CM energy
  
  pcm[3] = pgam.P(3)+_ev->GetOneByRole(2)->E();
  //std::cout << "pcm[3]=" << pcm[3] << std::endl;
  _gengam_w2 = std::pow(pcm[3], 2);
  for (int i=0; i<3; i++) {
    pcm[i] = pgam.P(i)+_ev->GetOneByRole(2)->P(i);
    _gengam_w2-= std::pow(pcm[i], 2);
  }
  //std::cout << "--> w2 = " << _gengam_w2 << std::endl;
  //exit(0);
  
  if (_gengam_w2<0.) {
    std::cerr << "[GamPomVMLL::GenGam] ERROR: W2 = " << _gengam_w2 << " < 0" << std::endl;
    return -1.;
  }
  pcm[4] = std::sqrt(_gengam_w2);
  
  // Determine weight (relative cross section) of the virtual vector meson
  //CALL SHSW (10, 0, SNGL (PCM (5)), 1.0)
  weight = 1./std::pow(1.+_q2/std::pow(_lambda, 2), _eprop);
  drlt = _xi*_q2/(std::pow(_lambda, 2)+_xi*_chi*_q2);
  //std::cout << "   q2 = " << _q2 << ", xi = " << _xi << ", lambda = " << _lambda << ", chi = " << _chi << std::endl;
  //std::cout << "--> heli = " << _event_heli << ", weight = " << weight << ", drlt = " << drlt << std::endl;
  
  if (_event_heli==0) {
    weight *= drlt;
    _photint_sweil += weight;
    _photint_sweil2 += std::pow(weight, 2);
  }
  else {
    _photint_sweit += weight;
    _photint_sweit2 += std::pow(weight, 2);
  }
  _photint_swei += weight;
  _photint_swei2 += std::pow(weight, 2);
  
  weight *= std::pow(_gengam_w2/_event_smax, 2.*_epsilw)/_event_propmx;
  
  //CALL SHSW (10, 1, SNGL (PCM (5)), WGHT)
  
  // Generate masses and t
  
  double genmxt = GenMXT(&wt);
  //std::cout << "weight = " << weight << ", genmxt=" << genmxt << ", " << _q2 << ", " << _gengam_w2 << std::endl;
  weight *= genmxt;
  //WGHT = WGHT*GENMXT (W2, Q2, .TRUE., DMXP,DMXV,T,DB,YHAT,POUT,WT)
  
  //CALL SHSW (10, 2, SNGL (PCM (5)), WGHT)
  
  if (weight>1.001) {
    std::cerr << "[GamPomVMLL::GenGam] ERROR: WEIGHT = " << weight << " > 1.001" << std::endl;	
    //CALL ERRLOG (62, 'W: GENGAM: WGHT > 1.001')
  }
  else if (wt>1.001) {
    std::cerr << "[GamPomVMLL::GenGam] ERROR: WT = " << wt << " > 1.001" << std::endl;	
    //CALL ERRLOG (63, 'W: GENGAM: WT > 1.001')
  }
  return weight;
  //return -1;
}
    
double
GamPomVMLL::ComputeWeight()
{
  // Event weight computation is performed in the gengam subroutine.
  // In this code, a first loop is introduced to extract the cross-section from an interpolation to
  // values observed at fixed Q2,t
  // Afterwards (line 20) the events loop is introduced and the common block containing the particles
  // kinematics is filled
  //
  // L. Forthomme, 19 sep 2014
  this->GenGam();

  return this->OneEvent();

  // From gdiffv.F
  // Set up for event generation (note: some of the parameters are set/calculated/changed during this step)
  //  this->GDIBeg();
  
  
  //this->GenEvtDi();
  //this->GenGam();
  return 0.;
}

double
GamPomVMLL::GenMXT(double* wght)
{
  double w;
  double genmxt, mxt;

  // From common block /CGDIF/:                                                                                                                     
  //   IFRAGP, IFRAGV, DEMINP, AMASSV, B0, AMXB0, WB0,                                                                                                
  //   ALPHA1, ALPH1M, WMIN, PE, PP                                                                                                                   
  // From common block /CMASS/:                                                                                                                     
  //   DMP, DMPI0, DMVM, DWVM, DMNST, DWNST
  //
  double tmin, tmax, tmean;
  double pcm1;
  double dmmin;

  //std::cout << "--> " << ifragp << " " << ifragv << std::endl;
  if (_genmxt_begin) {
    _genmxt_begin = false;
    _genmxt_bmin = 0;
    if ((ifragp!=1 and ifragp!=-1 and ifragp!=2) and ifragv==(int)0) {
      _genmxt_bmin = _b0+4.*_alpha1*log(_wmin/_wb0);
      //std::cout << "a-> " << _b0 << " " << _alpha1 << " " << _wmin << " " << _wb0 << std::endl;
    }
    else if ((ifragp==1 or ifragp==-1 or ifragp==2) and ifragv!=(int)0) {
      _genmxt_bmin = _b0+4.*_alpha1*log(4.*std::pow(_amxb0, 2)/(_wb0*_ecm));
      //std::cout << "b" << std::endl;
    }
    else {
      _genmxt_bmin = _b0+4.*_alpha1*log(_amxb0/_wb0);
      //std::cout << "--> " << _b0 << " " << _alpha1 << " " << _amxb0 << " " << _wb0 << std::endl;
    }
    _genmxt_bmin = std::max(_genmxt_bmin, 0.5);
    /*std::cout << "bmin=" << _genmxt_bmin << std::endl;
    exit(0);*/
  }

  w = std::sqrt(_gengam_w2);

  // Generate masses at p and VM vertex
  
  switch (ifragp) {
    case 0:
      _genmxt_dmxp = _dmp;
      break;
    case 1:
    case -1:
    case 2:
      _genmxt_dmxp = PXMass(_dmp+deminp, _ecm);
      break;
    default:
      _genmxt_dmxp = RanBW(_dmnst, _dwnst, _dmp+deminp, _dmnst+2.*_dwnst);
      break;
  }
  
  if (ifragv!=(int)0) {
    _genmxt_dmxv = VXMass(amassv, _ecm);
  }
  else {
    dmmin = _dmvm-3.*_dwvm;
    if (itypvm==100113 or itypvm==30113) {
      dmmin = std::max(dmmin, 1.2);
      }
    else if (itypvm==10333) {
      dmmin = std::max(dmmin, 1.4);
    }
    _genmxt_dmxv = RanBW(_dmvm, _dmvm, dmmin, _dmvm+10.*_dmvm);
    //std::cout << "dmxv = " << _genmxt_dmxv << ", dmvm = " << _dmvm << ", dmmin = " << dmmin << std::endl;
    if (_genmxt_dmxv<0) {
      //exit(0);
    }
  }
  
  // Return if generated masses are bigger than CM energy
  
  if (_genmxt_dmxp+_genmxt_dmxv>w-0.1) {
    _gengam_t = _genmxt_b = _gengam_yhat = _pcm3 = genmxt = 0.;
    /*std::cout << "genmxt: gen mass bigger than sqrt(s)" << std::endl;
    exit(0);*/
    *wght = 1.;
    return genmxt;
  }
    
  // Calculate slope parameter b
  // Generate t with e**(b*t) distribution
  
  _genmxt_b = _b0+4.*_alpha1*log(w/_wb0);
  if (ifragp==1 or ifragp==-1 or ifragp==2) {
    _genmxt_b -= 4.*_alph1m*log(_genmxt_dmxp/_amxb0);
  }
  if (ifragv!=(int)0) {
    _genmxt_b -= 4.*_alpha1*log(_genmxt_dmxv/_amxb0);
  }
  
  if (_genmxt_b<.5) _genmxt_b = .5;
  //std::cout << "before genert" << std::endl;
  _gengam_t = GenerT(0., _s, _genmxt_b, 1.*_anexp);
  //CALL GENERT (T, 0.0D0, S, B, 1D0*ANEXP)

  // Calculate actual minimal and maximal t for the generated masses
  // Note that t here is positive!
  // Formula (E.5) from Review of Particle Properties 1992, p. III.50
  // 1: gamma, 2: p, 3: VM(+X), 4: p remnant
  // The formula for Pcm1 is altered to take the imaginary photon mass
  // into account.
  
   pcm1 = std::sqrt(std::pow(_gengam_w2+_q2-std::pow(_dmp, 2), 2)+4.*_q2*std::pow(_dmp, 2))/w/2.;
  _pcm3 = std::sqrt((_gengam_w2-std::pow(_genmxt_dmxv+_genmxt_dmxp, 2))*(_gengam_w2-std::pow(_genmxt_dmxv-_genmxt_dmxp, 2)))/w/2.;
  tmean = ((-_q2-std::pow(_dmp, 2))*(std::pow(_genmxt_dmxv, 2)-std::pow(_genmxt_dmxp, 2))/_gengam_w2+_gengam_w2+_q2-std::pow(_dmp, 2)-std::pow(_genmxt_dmxv, 2)-std::pow(_genmxt_dmxp, 2))/2.;
  tmin = tmean-2.*pcm1*_pcm3;
  tmax = tmean+2.*pcm1*_pcm3;
  
  if (_gengam_t<=tmax and _gengam_t>=tmin) {
  //std::cout << "hahahahahahahahahaaaaa!" << std::endl;
    mxt = 1.;
    _gengam_yhat = (_gengam_t-tmin)/(4*pcm1*_pcm3);
  }
  else {
    mxt = 0.;
  }
  
  *wght = _genmxt_bmin/_genmxt_b;

  std::cout << "pcm1=" << pcm1 << std::endl;
  std::cout << "pcm3=" << _pcm3 << std::endl;
  std::cout << "w=" << w << std::endl;
  std::cout << "w2=" << _gengam_w2 << std::endl;
  std::cout << "dmxv=" << _genmxt_dmxv << std::endl;
  std::cout << "dmxp=" << _genmxt_dmxp << std::endl;
  std::cout << "t=" << _gengam_t << ", [" << tmin << ", <" << tmean << ">, " << tmax << "]"<< std::endl;
  std::cout << "genmxt: mxt=" << mxt << std::endl;

  //std::cout << _genmxt_bmin << " -> " << _genmxt_b << std::endl;
  genmxt = mxt*(*wght);
  //std::cout << "genmxt:" << genmxt << std::endl;
  
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
      m2 = std::exp(drand()*delta+lmin);
    }
    else { // Basic spectrum: 1/M^2(1+epsilon)
      m2 = std::pow(fact*drand()+m2min, -1./_epsilm);
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

  } while (1.6*drand()>y and iter<=100);

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
    m2 = exp(drand()*delta+lmin);
  }
  else { // Basic spectrum: 1/M^2(1+epsilon)
    m2 = std::pow(fact*drand()+m2min, -1./_epsilm);
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
  double glumas, gluwid, dmasvm, dmasgl, dmass;
  double pcmgam[4], pvmvm[4], pcmglu[4], pt[3];
  double dmu1, dmu2, dmu3, dmu4;
  double c1, c2, c3;
  double t, tmin, tmax, yhat;
  double ctheta, stheta, phi, pin, pout, pgamf, ptf, b;
  int i;

  int idahep[10][2];
  int npart;

  if (_fraggl_begin) {
    _fraggl_begin = false;
    glumas = Particle::GetMassFromPDGId(ifragv);
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
  dmu3 =  std::pow(dmasvm/dmass, 2);
  dmu4 =  std::pow(dmasgl/dmass, 2);

  c1 = 1.-(dmu1+dmu2+dmu3+dmu4)+(dmu1-dmu2)*(dmu3-dmu4);
  c2 = std::sqrt((std::pow(1.-dmu1-dmu2, 2)-4.*dmu1*dmu2)*(std::pow(1.-dmu3-dmu4, 2)-4.*dmu3*dmu4));
  c3 = (dmu3-dmu1)*(dmu4-dmu2)+(dmu1+dmu4-dmu2-dmu3)*(dmu1*dmu4-dmu2*dmu3);

  tmax = std::pow(_ppcms8[idifv][4], 2)*(c1+c2)/2.;
  tmin = std::pow(_ppcms8[idifv][4], 4)*c3/tmax;

  //CALL GENERT (T, TMIN, TMAX, B, 1D0*ANEXP)
  t = GenerT(tmin, tmax, b, 1.*_anexp);

  pin = dmass*std::sqrt(std::pow(1.-dmu1-dmu2, 2)-4.*dmu1*dmu2)/2.;
  pout = dmass*std::sqrt(std::pow(1.-dmu3-dmu4, 2)-4.*dmu3*dmu4)/2.;

  yhat = (t-tmin)/(4.*pin*pout);
  ctheta = 1.-2.*yhat;
  stheta = 2.*sqrt(yhat-std::pow(yhat, 2));

  // Calculate the 5-vectors of the VM and glueball in the gamma-pomeron CMS

  
  //void Lorenb(double u_, double ps_[], double pi_[], double pf_[]);
  //CALL LORENF8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PPCMS8 (1, IVVM), PCMGAM)

  pgamf = pout*ctheta/std::sqrt(std::pow(pcmgam[0], 2)+std::pow(pcmgam[1], 2)+std::pow(pcmgam[2], 2));

  phi = 2.*pi*drand();
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

  npart = _ev->NumParticles(); //FIXME
 
  idahep[idifv][0] = npart+1;
  idahep[idifv][1] = npart+2;
  //istat[idifv] = 2;

  // Glueball quantities
  /*iglue = npart+1;
  itype[iglue] = (int)ifragv;
	//double* glu;
	//glu3 = 
  //CALL LORENB8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PCMGLU, PPCMS8 (1, IGLUE))
  _ppcms8[iglue][4] = dmasgl;
  istat[iglue] = 1;
  idahep[iglue][0] = 0;
  idahep[iglue][1] = 0;
  mohep[iglue][0] = idifv;
  mohep[iglue][1] = 0;*/
  Particle glueball(43, ifragv);
  glueball.SetMother(_ev->GetOneByRole(4)); // 4 = idifv
  glueball.M(dmasgl);
  glueball.P(pcmglu);
  glueball.LorentzBoost(_ppcms8[idifv][4], _ppcms8[idifv]);
  //glueball.P() -> Lorentz boost
  glueball.status = 1;
  _ev->AddParticle(glueball);

  // Vector meson quantities
  /*ivm = npart+2;
  itype[ivm] = itypvm;
  //CALL LORENB8 (PPCMS8 (5, IDIFV), PPCMS8 (1, IDIFV), PCMVM, PPCMS8 (1, IVM))
  _ppcms8[ivm][4] = dmasvm;
  istat[ivm] = 1;
  idahep[ivm][0] = 0;
  idahep[ivm][1] = 0;
  mohep[ivm][0] = idifv;
  mohep[ivm][1] = 0;*/
  Particle VM(4, itypvm);
  VM.SetMother(_ev->GetOneByRole(42)); // 42 = idifv // FIXME???
  VM.M(dmasvm);
  VM.P(_pcmvm);
  VM.LorentzBoost(_ppcms8[idifv][4], _ppcms8[idifv]);
  VM.status = 1;
  _ev->AddParticle(VM);

  npart += 2;

  // Perform glueball decay 
  //CALL DECGLU (IGLUE)

  // Glueball can decay to K* K*bar => call DECK0
  //CALL DECK0 (IGLUE+1, NPART)

}

void
GamPomVMLL::GEPhot(double* q2_, int* heli_)
{
  std::vector<Particle> epa_result;
  std::vector<Particle>::iterator p;
  PhysicsBoundaries pb;

#ifdef DEBUG
  std::cout << "[GamPomVMLL::GEPhot] [DEBUG] Function called" << std::endl;
#endif
  
  pb.wmin = _wmin;
  pb.wmax = _wmax; // FIXME not present!
  pb.zmin = _ymin;
  pb.zmax = _ymax;
  pb.q2min = _q2min;
  pb.q2max = _q2max;

  //_ev->GetOneByRole(IP1)->Dump();
  //_ev->Dump();
  //std::cout << "Before EPA : " << _ev->NumParticles() << " particles" << std::endl;
  epa_result = EPA(_ev->GetOneByRole(IBE), _ev->GetOneByRole(IBP), _igammd, pb, q2_);
  for (p=epa_result.begin(); p!=epa_result.end(); p++) {
    if (p->role==2) { // electron/proton from EPA
      p->role = ISCE;
      *(_ev->GetOneByRole(ISCE)) = *p;
      continue;
    }
    if (p->role==3) { // photon from EPA
      p->role = IGAM;
      p->E(-1);
      *heli_ = p->helicity;
      _ev->AddParticle(*p);
      continue;
    }
  }
  //std::cout << "After EPA : " << _ev->NumParticles() << " particles" << std::endl;
  //_ev->Dump();

}

void
GamPomVMLL::GenDif()
{
  double ctheta, stheta;
  double r;
  double pcmvmx[5], pcmpx[5], pcmpom[5];
  double pt[3], ptf;
  double pgamf, pout, phi;

  // Check scattering angle in CMS
  if (_gengam_yhat<0.) {
    std::cerr << "[GamPomVMLL::GenDif] ERROR: YHAT < 0! YHAT = " << _gengam_yhat << std::endl;
    _gengam_yhat = 0.;
    //CALL ERRLOG (70, 'S: GENDIF: YHAT < 0!')
  }
  else if (_gengam_yhat>1.) {
    std::cerr << "[GamPomVMLL::GenDif] ERROR: YHAT > 1! YHAT = " << _gengam_yhat << std::endl;
    _gengam_yhat = 1.;
    //CALL ERRLOG (71, 'S: GENDIF: YHAT > 1!')
  }

  ctheta = 1.-2.*_gengam_yhat;
  stheta = 2.*std::sqrt(_gengam_yhat-std::pow(_gengam_yhat, 2));

  // Calculate the 5-vectors of the diffractive states in the CMS
  pgamf = pout*ctheta/std::sqrt(std::pow(_ppcms8[ivvm][0], 2)+std::pow(_ppcms8[ivvm][1], 2)+std::pow(_ppcms8[ivvm][2], 2));
  phi = 2.*pi*drand();
  pt[0] = -cos(phi)*_ppcms8[ivvm][2];
  pt[1] =  sin(phi)*_ppcms8[ivvm][2];
  pt[2] = -sin(phi)*_ppcms8[ivvm][1]+cos(phi)*_ppcms8[ivvm][0];
  ptf = pout*stheta/std::sqrt(std::pow(_ppcms8[ivvm][2], 2)+std::pow(pt[2], 2));

  pcmvmx[4] = _genmxt_dmxv;
  pcmvmx[3] = std::pow(_genmxt_dmxv, 2);
  for (int i=0; i<3; i++) {
    std::cout << "-> " << i << ", " << _ppcms8[ivvm][i] << std::endl;
    pcmvmx[i] = pgamf*_ppcms8[ivvm][i]+ptf*pt[i];
    pcmvmx[3]+= std::pow(pcmvmx[i], 2);
  }
  pcmvmx[3] = std::sqrt(pcmvmx[3]);

  if (fabs(std::pow(pout, 2)-std::pow(pcmvmx[0], 2)-std::pow(pcmvmx[1], 2)-std::pow(pcmvmx[2], 2))>std::pow(pout, 2)/100.) {
    std::cout << "[GamPomVMLL::GenDif] WARNING: POUT <> |PCMVMX|" << std::endl;
    //CALL ERRLOG (72, 'W: GENDIF: POUT <> |PCMVMX|')
    std::cout << "  POUT   = " << pout << std::endl;
    std::cout << "  PCMVMX = (" << pcmvmx[0] << ", " << pcmvmx[1] << ", " << pcmvmx[2] << ")" << std::endl;
  }

  pcmpx[4] = _genmxt_dmxp;
  pcmpx[3] = std::pow(_genmxt_dmxp, 2);
  for (int i=0; i<3; i++) {
    std::cout << i << " -> " << pcmvmx[i] << std::endl;
    pcmpx[i] = -pcmvmx[i];
    pcmpx[3]+= std::pow(pcmpx[i], 2);
  }
  pcmpx[3] = std::sqrt(pcmpx[3]);

  // Calculate momentum carried by the pomeron
  // the pomeron is thought to be a quasireal particle emitted by the proton
  // and absorbed by the virtual vector meson

  for (int i=0; i<4; i++) {
    pcmpom[i] = pcmvmx[i]-_ppcms8[ivvm][i];
  }
  pcmpom[4] = -std::sqrt(std::pow(pcmpom[0], 2)+std::pow(pcmpom[1], 2)+std::pow(pcmpom[2], 2)-std::pow(pcmpom[3], 2));

  // Virtual pomeron
  Particle pom(42, Particle::POMERON);
  pom.status = 3;
  pom.SetMother(_ev->GetOneByRole(2));
  pom.P(pcmpom[0], pcmpom[1], pcmpom[2], pcmpom[3]);
  //#ifdef DEBUG
  std::cout << "[GamPomVMLL::GenDif] [DEBUG] Virtual pomeron : " << pcmpom[4] << " <> " << pom.M() << std::endl;
  //#endif
  _ev->AddParticle(pom); // Pomeron

  // Diffractive proton state
  Particle dps(5, _ev->GetOneByRole(2)->pdgId);
  dps.status = 1;
  dps.SetMother(_ev->GetOneByRole(2));
  if (ifragp==1 or ifragp==-1 or ifragp==2) { // proton-dissociative case
    if (_genmxt_dmxp<1.48) dps.pdgId = (Particle::ParticleCode)12212;
    else if (_genmxt_dmxp<1.6) dps.pdgId = (Particle::ParticleCode)2124;
    else if (_genmxt_dmxp<1.9) {
      r = drand();
      if (r<.5) dps.pdgId = (Particle::ParticleCode)12216;
      else if (r<.83) dps.pdgId = (Particle::ParticleCode)22124;
      else dps.pdgId = (Particle::ParticleCode)42212;
    }
    else dps.pdgId = (Particle::ParticleCode)2210;
  }
  else if (ifragp!=0) dps.pdgId = (Particle::ParticleCode)abs(ifragp);
  dps.P(pcmpx[0], pcmpx[1], pcmpx[2], pcmpx[3]);
  std::cout << pcmpx[2] << std::endl;
  std::cout << "------> " << std::pow(pcmpx[3], 2)-std::pow(pcmpx[0], 2)-std::pow(pcmpx[1], 2)-std::pow(pcmpx[2], 2) << std::endl;
  dps.M(-1);
  //#ifdef DEBUG
  std::cout << "[GamPomVMLL::GenDif] [DEBUG] Diffractive proton : " << pcmpx[4] << " <> " << dps.M() << std::endl;
  //#endif
  _ev->AddParticle(dps);

  // Diffractive meson state
  Particle dms(8, (Particle::ParticleCode)itypvm);
  dms.SetMother(_ev->GetOneByRole(5));
  //FIXME dual mothers!
  if (ifragv!=0) {
    if (itypvm==22) dms.pdgId = Particle::REGGEON;
    else dms.pdgId = static_cast<Particle::ParticleCode>(10*((itypvm/10)%100));
  }
  dms.status = 1;
  dms.P(pcmvmx[0], pcmvmx[1], pcmvmx[2], pcmvmx[3]);
  //#ifdef DEBUG
  std::cout << "[GamPomVMLL::GenDif] [DEBUG] Diffractive meson : " << pcmvmx[4] << " <> " << dms.M() << std::endl;
  //#endif
  _ev->AddParticle(dms);
}

void
GamPomVMLL::FixPhot(Particle* phot_, Particle* ele_, double *q2_, Particle pel_, double egamma_)
{
  double y;
  double pgam[4], pe[4];

  y = egamma_/ele_->E();

  pe[3] = 0.;
  for (int i=0; i<3; i++) {
    pgam[i] = y*pel_.P(i);
    pe[i] = pel_.P(i)-pgam[i];
    pe[3] += std::pow(pe[i], 2);
  }

  pe[3] = std::sqrt(pe[3]+std::pow(_dme, 2));
  pgam[3] = pel_.E()-pe[3];
  *q2_ = std::pow(_dme, 2)+std::pow(y, 2)/(1.-y);
  (*phot_).P(pgam);
  (*ele_).P(pe);
}

void
GamPomVMLL::GenPhot(Particle* phot_, Particle* ele_, double *q2_, Particle pel_, double emin_, double q2max_)
{
  /*
   * Generate 1/k spectrum by methode of W. Bartel
   * R is fraction of EGAM of EMAX
   */
   
  double emax, r;
  double riter;
  double pe[5];
  double pgam[5];
   
  emax = ele_->P();
  riter = 0.;

  do {
    r = exp(drand()*log(emin_/emax));
    if (r>=1.) {
      std::cout << "[GamPomVMLL::GenPhot] Warning: R = " << r << " > 1" << std::endl;
    }
    pe[3] = 0.;
    *q2_ = 0.;
    for (int i=0; i<3; i++) {
      pgam[i] = r*pel_.P(i);
      pe[i] = pel_.P(i)-pgam[i];
      pe[3] += std::pow(pe[i], 2);
      *q2_ += std::pow(pgam[i], 2);
    }
    pe[3] = std::sqrt(pe[3]+std::pow(_dme, 2));
    pe[4] = _dme;
    pgam[3] = pel_.P(3)-pe[3];
    *q2_ -= std::pow(pgam[3], 2);
    
    riter += 1.;
  } while (fabs(*q2_)>fabs(q2max_) and q2max_!=0.);
  
  pgam[4] = -std::sqrt(fabs(*q2_));
   
  /* Fill CKINE:
      IHEL = 1
      FTRANS = -99999.0
      EPSIL = 0.0
   */
}

void
GamPomVMLL::VMFlux()
{
  
  if (_igammd==-1) {
    _vmflux_f = _vmflux_ft = 1.;
    _vmflux_fl = 0.;
    _vmflux_df = _vmflux_dft = _vmflux_dfl = 0.;
  }
  else if (_igammd==0 or _isum==0) {
    _vmflux_f = _vmflux_ft = 0.3;
    _vmflux_df = _vmflux_dft = 0.1;
    _vmflux_fl = _vmflux_dft = 0.;
  }
  
  if (_iacct>0) {
    _vmflux_ft = _dsumt/_isum*_iacct/_igent;
    _vmflux_dft = _vmflux_ft*std::sqrt((_qsumt/_dsumt-_dsumt/_isum)/(_isum-1)+(double)(_igent-_iacct)/_igent/_iacct);
  }
  else {
    _vmflux_ft = _vmflux_dft = 0.;
  }
  
  if (_iaccl>0) {
    _vmflux_fl = _dsuml/_isum*_iaccl/_igenl;
    _vmflux_dfl = _vmflux_fl*std::sqrt((_qsuml/_dsuml-_dsuml/_isum)/(_isum-1)+(double)(_igenl-_iaccl)/_igenl/_iaccl);
  }
  else {
    _vmflux_fl = _vmflux_dfl = 0.;
  }
  
  _vmflux_f = _vmflux_ft+_vmflux_fl;
  _vmflux_df = std::sqrt(std::pow(_vmflux_dft, 2)+std::pow(_vmflux_dfl, 2));
  
}

/*bool
GamPomVMLL::SetIncomingParticles(Particle, Particle){}

bool
GamPomVMLL::SetOutgoingParticles(int, int){}
*/
void
GamPomVMLL::FillKinematics(bool)
{
  /*
    ibe ---o------- isce
            \
             \ igam
              \_
              /xx\
	           (xxxx)-----
       	      \xx/
              / /
             / / ipom
            / /
    ibp ---ooo---- idifp
   */
  /*Particle beam_e(IBE);
    Particle beam_p(IBP);*/
  Particle scat_beam_e(ISCE);
  
  Particle diff_beam_p(IDIFP);
  Particle photon(IGAM);
  Particle pomeron(IPOM);
  Particle virt_vm(IVVM);
  Particle diff_qqbar(IDIFV);
  Particle glueball(IGLUE);
  Particle vm(IVM);
}
/*
void
GamPomVMLL::SetKinematics(Kinematics){}

void
GamPomVMLL::ComputeCMenergy(){}

void
GamPomVMLL::StoreEvent(std::ofstream*,double){}

void
GamPomVMLL::PrepareHadronisation(Particle *part_){}
*/
