#include "mcgen.h"

MCGen::MCGen(InputParameters ip_)
{
  unsigned int ndim;
  std::string topo;

#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif

  srand(time(0));
  _ip = ip_;

  if (_ip.p1mod<=2 && _ip.p2mod<=2) {
    topo = "ELASTIC proton/proton";
    ndim = 7;
  }
  else if (_ip.p1mod<=2 || _ip.p2mod<=2) {
    topo = "SINGLE-DISSOCIATIVE proton";
    ndim = 8;
  }
  else {
    topo = "DOUBLE-DISSOCIATIVE protons";
    ndim = 9;
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
  veg = new Vegas(ndim,f, &_ip);
}

MCGen::~MCGen()
{
  delete veg;
#ifdef DEBUG
  std::cout << "[MCGen::~MCGen] [DEBUG] MCGen destructed !" << std::endl;
#endif
}

void MCGen::ComputeXsection(double* xsec_, double *err_)
{
  std::cout << "[MCGen::ComputeXsection] Starting the computation of the process cross-section" << std::endl;
  veg->MyIntegrate(xsec_, err_);
  std::cout << "[MCGen::ComputeXsection] Total cross-section = " << *xsec_ << " +/- " << *err_ << " pb" << std::endl;
}

void MCGen::LaunchGeneration()
{
  veg->LaunchMyGeneration();
}

int i = 0;

double f(double* x_, size_t ndim_, void* params_) {
  //double etot, ptot;
  double ff;
  //int inp1pdg, inp2pdg;
  int outp1pdg, outp2pdg;
  InputParameters *p;
  GamGamKinematics cuts;

  i += 1;
  p = (InputParameters*)params_;

  //FIXME at some point introduce non head-on colliding beams ?

  //double inp1[3] = {0., 0.,  p->in1p};
  //double inp2[3] = {0., 0., -p->in2p}; //FIXME at this stage or at the InputParameters::in2p definition ?

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
#endif

  Particle *in1, *in2;

  switch(ndim_) {
  case 7:
  default:
    in1 = new Particle(1, 2212);
    in1->P(0., 0.,  p->in1p);
    in2 = new Particle(2, 2212);
    in2->P(0., 0., -p->in2p);
    outp1pdg = outp2pdg = 2212;
    cuts.kinematics = 1;
    break;
  case 8:
    in1 = new Particle(1, 2212);
    in1->P(0., 0.,  p->in1p);
    in2 = new Particle(2, 2212);
    in2->P(0., 0., -p->in2p);
    outp1pdg = 2;
    outp2pdg = 2212;
    cuts.kinematics = 2;
    break;
  case 9:
    in1 = new Particle(1, 2212);
    in1->P(0., 0.,  p->in1p);
    in2 = new Particle(2, 2212);
    in2->P(0., 0., -p->in2p);
    outp1pdg = outp2pdg = 2;
    cuts.kinematics = 3;
    break;
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
  gg.SetCuts(cuts);
  gg.SetIncomingKinematics(*in1, *in2);
  gg.SetOutgoingParticles(3, outp1pdg); // First outgoing proton
  gg.SetOutgoingParticles(5, outp2pdg); // Second outgoing proton
  gg.SetOutgoingParticles(6, p->pair); // Outgoing leptons
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
  
  if (ff<0.) {
    return 0.;
  }
  
  if (p->store) { // MC events generation
    gg.FillKinematics(p->symmetrise);
    gg.GetEvent()->Store(p->file);

    double t1min, t1max, t2min, t2max;

    if (p->file_debug->is_open()) {
      gg.GetT1extrema(t1min, t1max);
      gg.GetT2extrema(t2min, t2max);
      *(p->file_debug)
	<< (gg.GetEvent()->GetByRole(5)->p-gg.GetEvent()->GetByRole(1)->p)
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
	<< "\t" << gg.GetU1()
	<< "\t" << gg.GetU2()
	<< "\t" << gg.GetV1()
	<< "\t" << gg.GetV2()
	<< std::endl;
    }
  }

  delete in1;
  delete in2;
  return ff;
}

