#include "../include/mcgen.h"

MCGen::MCGen(InputParameters ip_)
{
#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif
  _ip = ip_;
  if (_ip.p1mod<=2 && _ip.p2mod<=2) {
    std::cout << "[MCGen::MCGen] [DEBUG] elastic case" << std::endl;
    _ndim = 6;
  }
  else if (_ip.p1mod<=2 || _ip.p2mod<=2) {
    std::cout << "[MCGen::MCGen] [DEBUG] single-dissociative case" << std::endl;
    _ndim = 7;
  }
  else {
    std::cout << "[MCGen::MCGen] [DEBUG] double-dissociative case" << std::endl;
    _ndim = 8;
  }
  std::cout << "[MCGen::MCGen] [DEBUG] single lepton's pT in range [" << _ip.minpt << ", " << _ip.maxpt << "]" << std::endl;
  veg = new Vegas(_ndim,f, &_ip);
}

MCGen::~MCGen()
{
#ifdef DEBUG
  std::cout << "[MCGen::~MCGen] [DEBUG] MCGen destructed !" << std::endl;
#endif
  delete veg;
}

void MCGen::LaunchGen(int count_)
{
  /*double x[7];
  for (unsigned int i=0; i<sizeof(x)/sizeof(double); i++) {
    x[i] = 0.3;
  }
  f(x, 8, &_ip);*/
  veg->LaunchIntegration();
  veg->LaunchGeneration(count_);
}

void MCGen::AnalyzePhaseSpace(std::string outputFile_) {
  const unsigned int numPoints = 1e2;
  const unsigned int numScans = 1e2;
  
  unsigned int v, i, j, k;
  std::stringstream outName, gpCmd, ss;
  std::ofstream of;
  double xpoint;
  double scanValues[numScans];

  double x[numPoints][_ndim]; 
  double xsec[numPoints][_ndim];
  Gnuplot plot;
  
  InputParameters ip = _ip;
  ip.debug = true;
  
  for (v=0; v<numScans; v++) {
    scanValues[v] = (double)v/(double)numScans;
    outName.str(""); outName << outputFile_ << "_x" << scanValues[v] << ".png";
    gpCmd.str(""); gpCmd << "plot ";
    for (i=0; i<abs(_ndim); i++) {
      for (j=0; j<numPoints; j++) {
        xpoint = 0.+1.*(double)j/(numPoints-1.);
        for (k=0; k<abs(_ndim); k++) {
          if (k!=i) {
            x[j][k] = scanValues[v];
          }
          else {
            x[j][k] = xpoint;
          }
        }
        xsec[j][i] = f(x[j], _ndim, &ip);
      }
      if (i!=0) gpCmd << ", ";
      gpCmd << "'" << outputFile_ << "' u 1:" << i+2 << " w l title 'x(" << i << ")'";
    }
    of.open(outputFile_.c_str());
    for (i=0; i<numPoints; i++) {
      of << (double)i/(numPoints-1.);
      for (j=0; j<abs(_ndim); j++) {
        of << "\t" << xsec[i][j];
      }
      of << std::endl;
    }
    ss.str(""); ss << "Variation of x(i), with all others x(j) = " << scanValues[v];
    plot.SetTitle(ss.str());
    plot.SetOutputFile(outName.str());
    plot.SetXAxisTitle("x(i)");
    plot.SetLogy();
    of.close();
    plot << gpCmd.str();
  }
}

double f(double* x_, size_t ndim_, void* params_) {
  //double etot, ptot;
  double q2min, q2max;
  double ff;
  int inp1pdg, inp2pdg, outp1pdg, outp2pdg;
  InputParameters *p;
  GamGam *gg;
  
  p = (InputParameters*)params_;
  
  //FIXME at some point introduce non head-on colliding beams ?
  double inp1[3] = {0., 0., p->in1p};
  double inp2[3] = {0., 0., p->in2p};
  
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
    q2min = 0.;
    q2max = 1.e5;
    inp1pdg = 2212;
    inp2pdg = 2212;
    outp1pdg = 2212;
    outp2pdg = 2212;
  }
  else {
    //FIXME for the inelastic case
    q2min = 0.;
    q2max = 1.e5;
    inp1pdg = 2212;
    inp2pdg = 2212;
    outp1pdg = 2212;
    outp2pdg = 2212;
  }
  Cuts cuts;
  cuts.mode = p->mcut;
  cuts.ptmin = p->minpt;
  cuts.ptmax = p->maxpt;
  
  gg = new GamGam(ndim_, q2min, q2max, 0, x_);
  gg->SetIncomingKinematics(1, inp1, inp1pdg);
  gg->SetIncomingKinematics(2, inp2, inp2pdg);
  gg->SetOutgoingParticles(3, outp1pdg); // First outgoing proton
  gg->SetOutgoingParticles(5, outp2pdg); // Second outgoing proton
  gg->SetOutgoingParticles(6, p->pair); // Outgoing leptons
  gg->SetCuts(cuts);
  if (!gg->IsKinematicsDefined()) {
    std::cout << "[f] [ERROR] Kinematics is not properly set" << std::endl;
    return 0.;
  }
  else {
    ff = gg->ComputeXsec();
  }
  if (ff!=0. && p->generation) { // MC events generation
    Particle ip1 = gg->GetParticle(1);
    Particle ip2 = gg->GetParticle(2);
    Particle op1 = gg->GetParticle(3);
    Particle ga1 = gg->GetParticle(41);
    Particle ga2 = gg->GetParticle(42);
    Particle op2 = gg->GetParticle(5);
    Particle ol1 = gg->GetParticle(6);
    Particle ol2 = gg->GetParticle(7);
    if (p->debug) {
    }
    
    /**(p->file) << ol1.e 
       << "\t" << ol1.px
       << "\t" << ol1.py
       << "\t" << ol1.pz
       << "\t" << ol1.pt
       << "\t" << ol1.m
       << "\t" << ol1.pdgId
       << std::endl;
    *(p->file) << ol2.e
       << "\t" << ol2.px
       << "\t" << ol2.py
       << "\t" << ol2.pz
       << "\t" << ol2.pt
       << "\t" << ol2.m
       << "\t" << ol2.pdgId
       << std::endl;*/
    *(p->file) << ga1.px << "\t" << ga1.py << "\t" << ga1.pz << "\t" << ga1.e << std::endl;
  }
  
  delete gg;
  return ff;
}

