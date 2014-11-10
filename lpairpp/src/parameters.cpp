#include "parameters.h"

Parameters::Parameters() :
  in1pdg(PROTON), in2pdg(PROTON),
  p1mod(2), p2mod(2),
  pair(MUON),
  mcut(0),
  minpt(0.5), maxpt(-1.),
  minenergy(1.), maxenergy(-1.),
  mintheta(5.), maxtheta(175.),
  //mintheta(0.), maxtheta(180.),
  minq2(0.), maxq2(1.e5),
  minmx(1.07), maxmx(320.),
  //ncvg(14000), itvg(10),
  ncvg(100000), itvg(10),
  ntreat(1), npoints(100),
  generation(true), store(false), debug(false),
  maxgen(1e5), ngen(0),
  gpdf(5), spdf(4), qpdf(12),
  hadroniser_max_trials(5),
  symmetrise(true)
{
  this->last_event = new Event();
  this->file = (std::ofstream*)NULL;
  this->hadroniser = (Hadroniser*)NULL;
  this->output_format = "lhe";
}

Parameters::~Parameters()
{
#ifdef DEBUG
  std::cout << "[Parameters::~Parameters] [DEBUG] Destructor called" << std::endl;
#endif
  delete last_event;
}

void Parameters::SetEtaRange(double etamin_, double etamax_)
{
  this->mintheta = 2.*atan(exp(-etamax_))/pi*180.;
  this->maxtheta = 2.*atan(exp(-etamin_))/pi*180.;
#ifdef DEBUG
  std::cout << "[Parameters::SetEtaRange] [DEBUG]"
	    << "\n\teta(min) = " << std::setw(5) << etamin_ << " -> theta(min) = " << this->mintheta
	    << "\n\teta(max) = " << std::setw(5) << etamax_ << " -> theta(max) = " << this->maxtheta
	    << std::endl;
#endif
}

void Parameters::Dump()
{
  std::string cutsmode, particles;

  switch(mcut) {
    case 1:
      cutsmode = "Vermaseren"; break;
    case 2:
      cutsmode = "both leptons"; break;
    case 3:
      cutsmode = "single lepton"; break;
    case 0:
    default:
      cutsmode = "none"; break;
  }
  switch(pair) {
    case 11:
      particles = "electrons"; break;
    case 13:
    default: 
      particles = "muons"; break;
    case 15:
      particles = "taus"; break;
  }
  const int wb = 65;
  const int wt = 40;
  int wp = wb-wt-2;
  std::cout 
    << std::left
    << "[Parameters::Dump] BEGINNING dump " << std::setfill('=') << std::setw(wb-32) << "" << std::endl
    << std::endl
    << " _" << std::setfill('_') << std::setw(wb) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "_ " << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Process to generate"  << std::setw(wp) << process->GetName() << std::endl
    << "| " << std::setw(wt) << "Events generation ? " << std::setw(wp) << generation << " |" << std::endl
    << "| " << std::setw(wt) << "Number of events to generate" << std::setw(wp) << maxgen << " |" << std::endl
    << "| " << std::setw(wt) << "Events storage ? " << std::setw(wp) << store << " |" << std::endl
    << "| " << std::setw(wt) << "Debugging mode ? " << std::setw(wp) << debug << " |" << std::endl
    << "| " << std::setw(wt) << "Output file opened ? " << std::setw(wp) << (file!=(std::ofstream*)NULL && file->is_open()) << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(wb) << "_/¯ INCOMING- AND OUTGOING KINEMATICS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Incoming protons-like particles " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Mode" << std::setw(3) << p1mod << ", " << std::setw(3) << p2mod << std::setw(wp-8) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Incoming particles" << std::setw(5) << in1pdg << ", " << std::setw(5) << in2pdg << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Momenta [GeV/c]" << std::setw(5) << in1p << ", " << std::setw(5) << in2p << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing leptons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Pair" << std::setw(2) << pair << " -> " << std::setw(wp-6) << particles << " |" << std::endl
    << "| " << std::setw(wt) << "Cuts mode" << std::setw(2) << mcut << " -> " << std::setw(wp-6) << cutsmode << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' pT in range [GeV/c]" << "[" << std::setw(4) << minpt << ", " << std::setw(4) << maxpt << "]" << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' energy in range [GeV]" << "[" << std::setw(4) << minenergy << ", " << std::setw(4) << maxenergy << "]" << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw(3) << maxtheta << "]" << std::setw(wp-10) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing remnants " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl;
  if (this->hadroniser!=(Hadroniser*)NULL)
    std::cout << "| " << std::setw(wt) << "Hadronisation algorithm" << std::setw(12) << hadroniser->GetName() << std::setw(wp-12) << "" << " |" << std::endl;
  std::cout << "| " << std::setw(wt) << "Minimal mass [GeV/c**2]" << std::setw(wp) << minmx << " |" << std::endl
    << "| " << std::setw(wt) << "Maximal mass [GeV/c**2]" << std::setw(wp) << maxmx << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(wb) << "_/¯ VEGAS INTEGRATION PARAMETERS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Maximum number of iterations" << std::setw(wp) << itvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of function calls" << std::setw(wp) << ncvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of points to try per bin" << std::setw(wp) << npoints << " |" << std::endl
    << "| " << std::setw(wt) << "Integration smoothed (TREAT) ? " << std::setw(wp) << ntreat << " |" << std::endl
    << "|_" << std::right << std::setfill('_') << std::setw(wb) << "_|" << std::left << std::endl
    //<< " -" << std::right << std::setfill('-') << std::setw(wb) << "- " << std::left << std::endl
    << std::endl
    << "[Parameters::Dump] END of dump " << std::setfill('=') << std::setw(wb-29) << "" << std::endl;
}

bool Parameters::ReadConfigFile(std::string inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_.c_str(), std::fstream::in);
  if (!f.is_open()) {
    return false;
  }
#ifdef DEBUG
  std::cout << "[Parameters::ReadConfigFile] [DEBUG] File " << inFile_ << " succesfully opened !" << std::endl
            << "======================================================" << std::endl
            << "Configuration file content : " << std::endl
            << "======================================================" << std::endl
            << std::left;
#endif
  while (f >> key >> value) {
    //std::cout << std::setw(60) << "[" << key << "] = " << value << std::endl;
    //if (strncmp(key.c_str(), "#")==0) continue; // FIXME need to ensure there is no extra space before !
    if (key=="IEND") {
      int iend = (int)atoi(value.c_str());
      if (iend>1) {
        this->generation = true;
      }
    }
    else if (key=="NCVG") {
      this->ncvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of function calls" << this->ncvg << std::endl;
#endif
    }
    else if (key=="NCSG") {
      this->npoints = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of points to probe" << this->npoints << std::endl;
#endif
    }
    else if (key=="ITVG") {
      this->itvg = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of Vegas iterations" << this->itvg << std::endl;
#endif
    }
    else if (key=="INPP") {
      this->in1p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * First incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="PMOD") {
      this->p1mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * First incoming particles' mode" << this->p1mod << " --> ";
      switch (this->p1mod) {
        case 1:          std::cout << "electron"; break;
        case 2: default: std::cout << "elastic proton"; break;
        case 11:         std::cout << "dissociating proton [structure functions]"; break;
        case 12:         std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
        case 101:        std::cout << "dissociating proton [parton model, only valence quarks]"; break;
        case 102:        std::cout << "dissociating proton [parton model, only sea quarks]"; break;
        case 103:        std::cout << "dissociating proton [parton model, valence and sea quarks]"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="INPE") {
      this->in2p = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Second incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
#endif
    }
    else if (key=="EMOD") {
      this->p2mod = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Second incoming particles' mode" << this->p2mod << " --> ";
      switch (this->p2mod) {
        case 1:           std::cout << "electron"; break;
        case 2: default:  std::cout << "elastic proton [EPA]"; break;
        case 11:          std::cout << "dissociating proton [structure functions]"; break;
        case 12:          std::cout << "dissociating proton [structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
        case 101:         std::cout << "dissociating proton [parton model, only valence quarks]"; break;
        case 102:         std::cout << "dissociating proton [parton model, only sea quarks]"; break;
        case 103:         std::cout << "dissociating proton [parton model, valence and sea quarks]"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PAIR") {
      this->pair = (ParticleId)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Outgoing leptons' PDG id   " << this->pair << " --> ";
      switch (this->pair) {
        case ELECTRON: default: std::cout << "electrons"; break;
        case MUON:              std::cout << "muons"; break;
        case TAU:               std::cout << "taus"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="MCUT") {
      this->mcut = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Set of cuts to apply on the total generation  " << this->mcut << " --> ";
      switch (this->mcut) {
        case 0: default: std::cout << "no cuts"; break;
        case 3:          std::cout << "cuts on at least one outgoing lepton"; break;
        case 2:          std::cout << "cuts on both the outgoing leptons"; break;
        case 1:          std::cout << "Vermaseren's hypothetical detector cuts"; break;
      }
      std::cout << std::endl;
#endif
    }
    else if (key=="PTCT") {
      this->minpt = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Single outgoing lepton's minimal transverse momentum" << this->minpt << " GeV/c" << std::endl;
#endif
    }
    else if (key=="ECUT") {
      this->minenergy = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Single outgoing lepton's minimal energy" << this->minenergy << " GeV" << std::endl;
#endif
    }
    else if (key=="NTRT") {
      this->ntreat = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of TREAT calls" << this->ntreat << std::endl;
#endif
    }
    else if (key=="NGEN") {
      this->maxgen = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Number of events to generate" << this->maxgen << std::endl;
#endif
    }
    else if (key=="THMN") {
      this->mintheta = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal polar production angle for the leptons" << this->mintheta << std::endl;
#endif
    }
    else if (key=="THMX") {
      this->maxtheta = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal polar production angle for the leptons" << this->maxtheta << std::endl;
#endif
    }
    else if (key=="Q2MN") {
      this->minq2 = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal Q^2 for the incoming photons" << this->minq2 << " GeV^2" << std::endl;
#endif
    }
    else if (key=="Q2MX") {
      this->maxq2 = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal Q^2 for the incoming photons" << this->maxq2 << " GeV^2" << std::endl;
#endif
    }
    else if (key=="MXMN") {
      this->minmx = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Minimal invariant mass of proton remnants" << this->minmx << " GeV/c^2" << std::endl;
#endif
    }
    else if (key=="MXMX") {
      this->maxmx = (double)atof(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Maximal invariant mass of proton remnants" << this->maxmx << " GeV/c^2" << std::endl;
#endif
    }
    else if (key=="GPDF") {
      this->gpdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * GPDF" << this->gpdf << std::endl;
#endif
    }
    else if (key=="SPDF") {
      this->spdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * SPDF" << this->spdf << std::endl;
#endif
    }
    else if (key=="QPDF") {
      this->qpdf = (double)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * QPDF" << this->qpdf << std::endl;
#endif
    }
    else {
      std::cout << std::setw(60) << "[Parameters::ReadConfigFile] <WARNING> Unrecognized argument : [" << key << "] = " << value << std::endl;
    }
  }
  f.close();
  std::cout << std::right << "======================================================" << std::endl;
  return true;
}

bool Parameters::StoreConfigFile(std::string outFile_)
{
  std::ofstream f;
  f.open(outFile_.c_str(), std::fstream::out | std::fstream::trunc);
  if (!f.is_open()) {
    return false;
  }
  // ...
  if (this->itvg!=-1) f << "ITVG  " << this->itvg << std::endl;
  if (this->minenergy!=-1) f << "ECUT  " << this->minenergy << std::endl;
  if (this->minenergy!=-1) f << "PTCT  " << this->minpt << std::endl;
  if (this->ntreat!=-1) f << "NTRT  " << this->ntreat << std::endl;
  // ...
  f.close();
  return true;
}

