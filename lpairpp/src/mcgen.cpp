#include "mcgen.h"

MCGen::MCGen(InputParameters ip_)
{
  std::string topo;
#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif
  srand(time(0));

  _ip = ip_;

  if (_ip.p1mod<=2 && _ip.p2mod<=2) {
    topo = "ELASTIC proton/proton";
    _ndim = 7;
  }
  else if (_ip.p1mod<=2 || _ip.p2mod<=2) {
    topo = "SINGLE-DISSOCIATIVE proton";
    _ndim = 8;
  }
  else {
    topo = "DOUBLE-DISSOCIATIVE protons";
    _ndim = 8;
  }
#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] Considered topology : " << topo << " case" << std::endl;
#endif

#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] Cuts mode : " << _ip.mcut << std::endl;
  switch(_ip.mcut) {
    case 1:
    case 2:
      std::cout << "[MCGen::MCGen] [DEBUG] Single leptons' transverse momentum condition : ";
      if (_ip.minpt<=0.) {
        std::cout << "no pT cut" << std::endl;
        break;
      }
      if (_ip.maxpt>0.) {
        std::cout << "pT in range [" 
                  << _ip.minpt << " GeV/c, " 
                  << _ip.maxpt << " GeV/c]" 
                  << std::endl;
        break;
      }
      std::cout << "pT > " << _ip.minpt << " GeV/c";
      if (_ip.mcut==1) {
        std::cout << " for at least one lepton" << std::endl;
      }
      else {
        std::cout << " for both the leptons" << std::endl;
      }
      break;
    case 0:
    default:
      std::cout << "[MCGen::MCGen] [DEBUG] No cuts applied on the total cross section" << std::endl;
      break;
  }
#endif
}

MCGen::~MCGen()
{
#ifdef DEBUG
  std::cout << "[MCGen::~MCGen] [DEBUG] MCGen destructed !" << std::endl;
#endif
}

void MCGen::Test()
{
  double x[7];
  int num = 5;
  double ext;
  _ip.mcut = 2;
  _ip.minpt = 0.1;
  //_ip.Dump();
  /*for (unsigned int i=0; i<sizeof(x)/sizeof(double); i++) {
    x[i] = 0.1;
  }
  std::cout << f(x, (int)(sizeof(x)/sizeof(double)), &_ip) << std::endl;*/
  ext = 1./num;
  int i1,i2,i3,i4,i5,i6,i7;
  for (i1=0; i1<=num; i1++) {
    x[0] = i1*ext;
    for (i2=0; i2<=num; i2++) {
      x[1] = i2*ext;
      for (i3=0; i3<=num; i3++) {
        x[2] = i3*ext;
        for (i4=0; i4<=num; i4++) {
          x[3] = i4*ext;
          for (i5=0; i5<=num; i5++) {
            x[4] = i5*ext;
            for (i6=0; i6<=num; i6++) {
              x[5] = i6*ext;
              for (i7=0; i7<=num; i7++) {
                x[6] = i7*ext;
                std::cout << " -- " << x[0]
                          << "\t" << x[1]
                          << "\t" << x[2]
                          << "\t" << x[3]
                          << "\t" << x[4]
                          << "\t" << x[5]
                          << "\t" << x[6]
                          << std::endl;
                f(x, (int)(sizeof(x)/sizeof(double)), &_ip);
              }
            }
          }
        }
      }
    }
  }
}

void MCGen::ComputeXsection(double* xsec_, double *err_)
{
  std::cout << "[MCGen::ComputeXsection] Starting the computation of the process cross-section" << std::endl;
  veg = new Vegas(_ndim,f, &_ip);
  veg->Integrate(xsec_, err_);
  std::cout << "[MCGen::ComputeXsection] Total cross-section = " << *xsec_ << " +/- " << *err_ << " pb" << std::endl;
  delete veg;
}

void MCGen::LaunchGen()
{
  veg = new Vegas(_ndim,f, &_ip);
  //veg->LaunchGeneration();
  veg->LaunchMyGeneration();
  delete veg;
}

int i = 0;

double f(double* x_, size_t ndim_, void* params_) {
  //double etot, ptot;
  double ff;
  int inp1pdg, inp2pdg, outp1pdg, outp2pdg;
  InputParameters *p;
  GamGamKinematics cuts;

  i += 1;
  p = (InputParameters*)params_;

  //FIXME at some point introduce non head-on colliding beams ?
  double inp1[3] = {0., 0.,  p->in1p};
  double inp2[3] = {0., 0., -p->in2p}; //FIXME at this stage or at the InputParameters::in2p definition ?

  ff = 0.;

#ifdef DEBUG
  std::cout << "=====================================" << std::endl;
  std::cout << "function f called ; some parameters :\n"
            << "\n  pz(p1) = " << p->in1p
            << "\n  pz(p2) = " << p->in2p
            << "\n   f(p1) = " << p->p1mod
            << "\n   f(p2) = " << p->p2mod
            << std::endl;
  std::cout << "=====================================" << std::endl;
#else
  //std::cout << "=============== f called ===============" << std::endl;
#endif

  if (p->p1mod<=2 && p->p2mod<=2) {
    inp1pdg = 2212;
    inp2pdg = 2212;
    outp1pdg = 2212;
    outp2pdg = 2212;
    cuts.kinematics = 1;
  }
  else {
    //FIXME for the inelastic case
    inp1pdg = 2212;
    inp2pdg = 2212;
    outp1pdg = 2212;
    outp2pdg = 2212;
    cuts.kinematics = 2; //FIXME
  }
  cuts.q2min = p->minq2;
  cuts.q2max = p->maxq2;
  cuts.mode = p->mcut;
  cuts.ptmin = p->minpt;
  cuts.ptmax = p->maxpt;
  cuts.thetamin = p->mintheta;
  cuts.thetamax = p->maxtheta;
  cuts.emin = p->minenergy;
  cuts.emax = p->maxenergy;
  cuts.mxmin = p->minmx;
  cuts.mxmax = p->maxmx;
  
  GamGam gg(ndim_, 0, x_);
  gg.SetIncomingKinematics(1, inp1, inp1pdg);
  gg.SetIncomingKinematics(2, inp2, inp2pdg);
  gg.SetOutgoingParticles(3, outp1pdg); // First outgoing proton
  gg.SetOutgoingParticles(5, outp2pdg); // Second outgoing proton
  gg.SetOutgoingParticles(6, p->pair); // Outgoing leptons
  gg.SetCuts(cuts);
  //cuts.Dump();
  //std::cout << "After the sets..." << std::endl;
  if (!gg.IsKinematicsDefined()) {
    std::cout << "[f] [ERROR] Kinematics is not properly set" << std::endl;
    return 0.;
  }
  ff = gg.ComputeXsec();
#ifdef DEBUG
  if (i==1) {
    std::cout << "--> f at first step = " << ff << std::endl;
    std::cout << "=========================" << std::endl;
    cuts.Dump();
    std::cout << "=========================" << std::endl;
  }
#endif
  
  /*double tmin,tmax;
  gg->GetT1extrema(tmin, tmax);
  std::cout << "T1\t" << gg->GetT1() << "\t" << tmin << "\t" << tmax << std::endl;
  gg->GetT2extrema(tmin, tmax);
  std::cout << "T2\t" << gg->GetT2() << "\t" << tmin << "\t" << tmax << std::endl;
  std::cout << "F\t" << ff << std::endl;*/

  if (ff<0.) {
    /*if (p->store)
      std::cout << "/////////////////// " << ff << std::endl;*/
    return 0.;
  }
  
  //if (p->store && p->ngen<p->maxgen) { // MC events generation
  if (p->store) { // MC events generation
    //p->ngen += 1;
    //std::cout << ">> new value for ngen = " << p->ngen << std::endl;
    /*if (p->store && p->ngen==0) {
      std::cout << "number of events generated : " << p->ngen << std::endl;
    }
    if ((p->ngen)%10000==0) {
      //std::cout << "[f] ngen = " << p->ngen << std::endl;
    }*/
    //std::cout << "number of events generated : " << p->ngen << std::endl;
    gg.FillKinematics(p->symmetrise);
    gg.GetEvent()->Store(p->file);

    double t1min, t1max, t2min, t2max;
    gg.GetT1extrema(t1min, t1max);
    gg.GetT2extrema(t2min, t2max);

    if (p->file_debug->is_open()) {
      *(p->file_debug) << (gg.GetEvent()->GetByRole(5)->p-gg.GetEvent()->GetByRole(1)->p)
                       << "\t" << gg.GetEvent()->GetByRole(3)->p
                       << "\t" << gg.GetEvent()->GetByRole(5)->p
		       << "\t" << gg.GetT1()
		       << "\t" << t1min
		       << "\t" << t1max
		       << "\t" << gg.GetT2()
		       << "\t" << t2min
		       << "\t" << t2max
		       << "\t" << gg.GetS1()
		       << "\t" << gg.GetS2()
		       << "\t" << gg.GetD3()
                       << std::endl;
    }
    //std::cout << "=============================" << std::endl;
    //gg.GetEvent()->GetByRole(41)->Dump();
    //gg.GetEvent()->GetByRole(42)->Dump();
    //gg.GetEvent()->Dump();
    //*(p->file) << ga1.px << "\t" << ga1.py << "\t" << ga1.pz << "\t" << ga1.e << std::endl;
  }
  return ff;
}

