#include "inputparameters.h"

InputParameters::InputParameters() :
  p1mod(2), p2mod(2),
  pair(13),
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
  maxgen(1e5),
  gpdf(5), spdf(4), qpdf(12),
  symmetrise(true)
{
  /*for (int i=0; i<MAX_HISTOS; i++) {
    this->plot[i] = new Gnuplot();
  }*/
}

InputParameters::~InputParameters()
{
  /*for (int i=0; i<MAX_HISTOS; i++) {
    delete[] this->plot[i];
  }//FIXME ???*/
  //delete[] this->plot;
}

HEPRUP
InputParameters::GetEventsInfo()
{
  HEPRUP out(1);
  if (this->p1mod==1) out.idbmup[0] = 11;
  if (this->p2mod==1) out.idbmup[1] = 11;
  out.ebmup[0] = this->in1p;
  out.ebmup[1] = this->in2p;

  return out;
}

void InputParameters::SetEtaRange(double etamin_, double etamax_)
{
  this->mintheta = 2.*atan(exp(-etamax_))/pi*180.;
  this->maxtheta = 2.*atan(exp(-etamin_))/pi*180.;
#ifdef DEBUG
  std::cout << "[InputParameters::SetEtaRange] [DEBUG]"
	    << "\n\teta(min) = " << std::setw(5) << etamin_ << " -> theta(min) = " << this->mintheta
	    << "\n\teta(max) = " << std::setw(5) << etamax_ << " -> theta(max) = " << this->maxtheta
	    << std::endl;
#endif
}

void InputParameters::Dump()
{
  std::string cutsmode, particles;

  switch(mcut) {
    case 1:
      cutsmode = "Vermaseren"; break;
    case 2:
      cutsmode = "both"; break;
    case 3:
      cutsmode = "single"; break;
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
  std::cout 
    << std::left
    << "[InputParameters::Dump] BEGINNING dump ===============" << std::endl << std::endl
    << " _" << std::setfill('_') << std::setw(52) << "_/¯ INCOMING- AND OUTGOING KINEMATICS ¯\\_" << std::setfill(' ') << "_ " << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Incoming protons-like particles " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Mode" << std::setw(4) << p1mod << ", " << std::setw(4) << p2mod << " |" << std::endl
    << "| " << std::setw(40) << "Momenta [GeV/c]" << std::setw(4) << in1p << ", " << std::setw(4) << in2p << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Outgoing leptons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Pair" << std::setw(2) << pair << " " << std::setw(7) << particles << " |" << std::endl
    << "| " << std::setw(40) << "Cuts mode" << std::setw(1) << mcut << " (" << std::setw(6) << cutsmode << ")" << " |" << std::endl
    << "| " << std::setw(40) << "pT [GeV/c]" << "[" << std::setw(3) << minpt << ", " << std::setw(3) << maxpt << "]" << " |" << std::endl
    << "| " << std::setw(40) << "Energy [GeV]" << "[" << std::setw(3) << minenergy << ", " << std::setw(3) << maxenergy << "]" << " |" << std::endl
    << "| " << std::setw(40) << "Polar angle theta [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw(3) << maxtheta << "]" << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(50) << " Outgoing remnants " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Minimal mass [GeV/c**2]" << std::setw(10) << minmx << " |" << std::endl
    << "| " << std::setw(40) << "Maximal mass [GeV/c**2]" << std::setw(10) << maxmx << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(52) << "_/¯ VEGAS INTEGRATION PARAMETERS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Maximum number of iterations" << std::setw(10) << itvg << " |" << std::endl
    << "| " << std::setw(40) << "Number of function calls" << std::setw(10) << ncvg << " |" << std::endl
    << "| " << std::setw(40) << "Number of points to try per bin" << std::setw(10) << npoints << " |" << std::endl
    << "| " << std::setw(40) << "Is the integration smoothed (TREAT) ? " << std::setw(10) << ntreat << " |" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(52) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(52) << " |" << std::left << std::endl
    << "| " << std::setw(40) << "Events generation ? " << std::setw(10) << generation << " |" << std::endl
    << "| " << std::setw(40) << "Number of events to generate" << std::setw(10) << maxgen << " |" << std::endl
    << "| " << std::setw(40) << "Events storage ? " << std::setw(10) << store << " |" << std::endl
    << "| " << std::setw(40) << "Debugging mode ? " << std::setw(10) << debug << " |" << std::endl
    << "| " << std::setw(40) << "Is Output file opened ? " << std::setw(10) << file->is_open() << " |" << std::endl
    //<< "| " << std::setw(40) << "Is Debug file opened ? " << std::setw(10) << file_debug->is_open() << " |" << std::endl
    << "|_" << std::right << std::setfill('_') << std::setw(52) << "_|" << std::left << std::endl
    //<< " -" << std::right << std::setfill('-') << std::setw(52) << "- " << std::left << std::endl
    << std::endl
    << "[InputParameters::Dump] END of dump ==================" << std::endl;
}

bool InputParameters::ReadConfigFile(std::string inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_.c_str(), std::fstream::in);
  if (!f.is_open()) {
    return false;
  }
#ifdef DEBUG
  std::cout << "[InputParameters::ReadConfigFile] [DEBUG] File " << inFile_ << " succesfully opened !" << std::endl
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
      this->pair = (int)atoi(value.c_str());
#ifdef DEBUG
      std::cout << std::setw(60) << " * Outgoing leptons' PDG id   " << this->pair << " --> ";
      switch (this->pair) {
        case 11: default: std::cout << "electrons"; break;
        case 13:          std::cout << "muons"; break;
        case 15:          std::cout << "taus"; break;
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
      std::cout << std::setw(60) << "[InputParameters::ReadConfigFile] <WARNING> Unrecognized argument : [" << key << "] = " << value << std::endl;
    }
  }
  f.close();
  std::cout << std::right << "======================================================" << std::endl;
  return true;
}

bool InputParameters::StoreConfigFile(std::string outFile_)
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

