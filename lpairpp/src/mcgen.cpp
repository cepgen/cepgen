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

  // For debugging purposes
  
  /*double x[ndim];
  for (unsigned int i=0; i<ndim; i++) {
    x[i] = 0.4;
  }
  std::cout << ">>> " << f(x, ndim, (void*)&ip_) << std::endl;
  exit(0);*/
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
  veg->Integrate(xsec_, err_);
  std::cout << "[MCGen::ComputeXsection] Total cross-section = " << *xsec_ << " +/- " << *err_ << " pb" << std::endl;
}

void MCGen::LaunchGeneration()
{
  veg->Generate();
}

int i = 0;

double f(double* x_, size_t ndim_, void* params_) {
  double ff;
  int outp1pdg, outp2pdg;
  InputParameters *p;
  GamGamKinematics kin;

  i += 1;
  p = (InputParameters*)params_;

  //FIXME at some point introduce non head-on colliding beams ?

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

  //FIXME electrons ?

  in1 = new Particle(1, 2212);
  in1->charge = 1;
  in1->P(0., 0.,  p->in1p);

  in2 = new Particle(2, 2212);
  in2->charge = 1;
  in2->P(0., 0., -p->in2p);
  
  switch(ndim_) {
  case 7:
  default:
    outp1pdg = outp2pdg = 2212;
    kin.kinematics = 1;
    break;
  case 8:
    outp1pdg = 2;
    outp2pdg = 2212;
    kin.kinematics = 2;
    break;
  case 9:
    outp1pdg = outp2pdg = 2;
    kin.kinematics = 3;
    break;
  }

  kin.q2min = p->minq2;
  kin.q2max = p->maxq2;
  //kin.q2min = 0; //FIXME
  //kin.q2max = -1; //FIXME
  kin.mode = p->mcut;
  kin.ptmin = p->minpt;
  kin.ptmax = p->maxpt;
  kin.thetamin = p->mintheta;
  kin.thetamax = p->maxtheta;
  kin.emin = p->minenergy;
  kin.emax = p->maxenergy;
  kin.mxmin = p->minmx;
  kin.mxmax = p->maxmx;
  //kin.mxmin = 0; //FIXME
  //kin.mxmax = 100000; //FIXME

  GamGam gg(ndim_, 0, x_);
  gg.SetKinematics(kin);
  gg.SetIncomingKinematics(*in1, *in2);
  gg.SetOutgoingParticles(3, outp1pdg); // First outgoing proton
  gg.SetOutgoingParticles(5, outp2pdg); // Second outgoing proton
  gg.SetOutgoingParticles(6, p->pair); // Outgoing leptons
  if (!gg.IsKinematicsDefined()) {
    std::cout << "[f] [ERROR] Kinematics is not properly set" << std::endl;
    return 0.;
  }
  ff = gg.ComputeWeight();
#ifdef DEBUG
  if (i==1) {
    std::cout << "--> f at first step = " << ff << std::endl;
    std::cout << "=========================" << std::endl;
    kin.Dump();
    std::cout << "=========================" << std::endl;
  }
#endif
  
  if (ff<0.) {
    return 0.;
  }
  
  if (p->store) { // MC events generation
    gg.FillKinematics(false);
    if (kin.kinematics>1) {
      gg.PrepareHadronisation(gg.GetEvent()->GetOneByRole(3));
      /*if (kin.kinematics==3) {
	gg.PrepareHadronisation(gg.GetEvent()->GetOneByRole(5));
	}*/
      //if (algo_=="pythia6") {
      //std::cout << "hadronisation using pythia6" << std::endl;
      Pythia6Hadroniser py;
      py.Hadronise(gg.GetEvent());
      //Jetset7Hadroniser js;
      //js.Hadronise(gg.GetEvent());
      //}
    }
    
    gg.GetEvent()->Dump();
    /*std::vector<Particle*> pp = gg.GetEvent()->GetParticles();
    std::vector<Particle*>::iterator p;
    for (p=pp.begin(); p!=pp.end(); p++) {
      std::cout << "--> " << (*p)->pdgId << ", " << (*p)->role << std::endl;
      }*/
    //gg.GetEvent()->Store(p->file);

    /*double t1min, t1max, t2min, t2max;

    if (p->file_debug->is_open()) {
      gg.GetT1extrema(t1min, t1max);
      gg.GetT2extrema(t2min, t2max);
      *(p->file_debug)
	<< (gg.GetEvent()->GetByRole(5)->p-gg.GetEvent()->GetOneByRole(1)->p)
	<< "\t" << gg.GetEvent()->GetOneByRole(3)->p
	<< "\t" << gg.GetEvent()->GetOneByRole(5)->p
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
	}*/
  }

  delete in1;
  delete in2;

  return ff;
}

