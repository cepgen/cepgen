#include "../include/mcgen.h"

MCGen::MCGen(InputParameters ip_)
{
#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif
  _ip = ip_;
  std::cout << "[MCGen::MCGen] [DEBUG] Considered case : " << std::endl;
  if (_ip.p1mod<=2 && _ip.p2mod<=2) {
    std::cout << "\telastic case" << std::endl;
    _ndim = 7;
  }
  else if (_ip.p1mod<=2 || _ip.p2mod<=2) {
    std::cout << "\tsingle-dissociative case" << std::endl;
    _ndim = 8;
  }
  else {
    std::cout << "\tdouble-dissociative case" << std::endl;
    _ndim = 9;
  }
  veg = new Vegas(_ndim,f, &_ip);
  std::cout << "[MCGen::MCGen] [DEBUG] Cuts mode : " << _ip.mcut << std::endl;
  switch(_ip.mcut) {
    case 0:
      std::cout << "[MCGen::MCGen] [DEBUG] No cuts applied on the total cross section" << std::endl;
      break;
    case 2:
      std::cout << "[MCGen::MCGen] [DEBUG] single lepton's pT in range [" << _ip.minpt << ", " << _ip.maxpt << "]" << std::endl;
      break;
    default:
      break;
  }
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
  std::cout << f(x, (int)(sizeof(x)/sizeof(double)), &_ip) << std::endl;*/
  if (_ip.debug) {
    _ip.plot[0]->SetHistogram(200, -50., 50., "px");
    _ip.plot[1]->SetHistogram(200, -50., 50., "py");
    _ip.plot[2]->SetHistogram(200, -50., 50., "pz");
    _ip.plot[3]->SetHistogram(100, 0., 100., "pt");
    _ip.plot[4]->SetHistogram(100, 0., 100., "invm");
    _ip.plot[5]->SetHistogram(100, 0., 100., "ptpair");
    _ip.plot[6]->SetHistogram(1000, -pi, pi, "sinth6");
    _ip.plot[7]->SetHistogram(1000, -pi, pi, "cosphi6");
    _ip.plot[8]->SetHistogram(1000, -pi, pi, "sinphi6");
  }
  veg->LaunchIntegration();
  veg->LaunchGeneration(count_);
  if (_ip.debug) {
    _ip.plot[0]->DrawHistogram();
    _ip.plot[1]->DrawHistogram();
    _ip.plot[2]->DrawHistogram();
    _ip.plot[3]->SetLogy();
    _ip.plot[3]->DrawHistogram();
    _ip.plot[4]->SetLogy();
    _ip.plot[4]->DrawHistogram();
    _ip.plot[5]->DrawHistogram();
    _ip.plot[6]->DrawHistogram();
    _ip.plot[7]->DrawHistogram();
    _ip.plot[8]->DrawHistogram();
  }
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
    for (i=0; i<(unsigned int)abs(_ndim); i++) {
      for (j=0; j<numPoints; j++) {
        xpoint = 0.+1.*(double)j/(numPoints-1.);
        for (k=0; k<(unsigned int)abs(_ndim); k++) {
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
      for (j=0; j<(unsigned int)abs(_ndim); j++) {
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
  ff = gg->ComputeXsec();

  if (ff!=0. && p->generation) { // MC events generation
    //gg->FillKinematics();
    
    Particle ip1 = gg->GetParticle(1);
    Particle ip2 = gg->GetParticle(2);
    Particle op1 = gg->GetParticle(3);
    Particle ga1 = gg->GetParticle(41);
    Particle ga2 = gg->GetParticle(42);
    Particle op2 = gg->GetParticle(5);
    Particle ol1 = gg->GetParticle(6);
    Particle ol2 = gg->GetParticle(7);
    if (p->debug) {
      if (ol1.pdgId<0) {
        p->plot[0]->Fill(ol1.px);
        p->plot[1]->Fill(ol1.py);
        p->plot[2]->Fill(ol1.pz);
        p->plot[3]->Fill(ol1.pt);
      }
      else {
        p->plot[0]->Fill(ol2.px);
        p->plot[1]->Fill(ol2.py);
        p->plot[2]->Fill(ol2.pz);
        p->plot[3]->Fill(ol2.pt);
      }
      double mass = std::sqrt(std::pow(ol1.m,2)+std::pow(ol2.m,2)+2*(ol1.e*ol2.e-ol1.px*ol2.px-ol1.py*ol2.py-ol1.pz*ol2.pz));
      p->plot[4]->Fill(mass);
      double ptpair = std::sqrt(std::pow(ol1.px+ol2.px, 2)+std::pow(ol1.py+ol2.py, 2));
      p->plot[5]->Fill(ptpair);
    }
    p->ngen += 1;
    if ((p->ngen)%10000==0) {
      std::cout << "[f] ngen = " << p->ngen << std::endl;
    }
    *(p->file) << ol1.e
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
       << std::endl;
    //*(p->file) << ga1.px << "\t" << ga1.py << "\t" << ga1.pz << "\t" << ga1.e << std::endl;
  }

  delete gg;
  return ff;
}

