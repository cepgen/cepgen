#include "parameters.h"

Parameters::Parameters() :
  in1pdg(Particle::Proton), in2pdg(Particle::Proton),
  remnant_mode(Process::SuriYennie),
  pair(Particle::Muon),
  process_mode(Process::ElasticElastic),
  mcut(0),
  minpt(0.5), maxpt(-1.),
  minenergy(1.), maxenergy(-1.),
  //minrapidity(-5.), maxrapidity(5.),
  mineta(-5.), maxeta(5.),
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
  PrintDebug("Destructor called");
  
  delete last_event;
}

void Parameters::SetThetaRange(double thetamin_, double thetamax_)
{
  this->mineta = -log(tan(thetamax_/180.*pi/2.));
  this->maxeta = -log(tan(thetamin_/180.*pi/2.));

  PrintDebug(Form("eta(min) = %5.2f => theta(min) = %5.2f"
                  "eta(max) = %5.2f => theta(max) = %5.2f",
                  mineta, thetamin_, maxeta, thetamax_));
}

void Parameters::Dump()
{
  std::string cutsmode, particles, pmode;
  std::ostringstream os;

  switch (mcut) {
    case 1:           cutsmode = "Vermaseren"; break;
    case 2:           cutsmode = "both leptons"; break;
    case 3:           cutsmode = "single lepton"; break;
    case 0:  default: cutsmode = "none"; break;
  }
  switch (pair) {
    case 11:          particles = "electrons"; break;
    case 13: default: particles = "muons"; break;
    case 15:          particles = "taus"; break;
  }
  switch (process_mode) {
    case 1:  default: pmode = "elastic-elastic"; break;
    case 2:           pmode = "elastic-inelastic"; break;
    case 3:           pmode = "inelastic-elastic"; break;
    case 4:           pmode = "inelastic-inelastic"; break;
  }
  const int wb = 65;
  const int wt = 40;
  int wp = wb-wt-2;
  os 
    << std::left
    << " _" << std::setfill('_') << std::setw(wb) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "_ " << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Process to generate"  << std::setw(wp) << process->GetName() << std::endl
    << "| " << std::setw(wt) << "Events generation ? " << std::setw(wp) << generation << " |" << std::endl
    << "| " << std::setw(wt) << "Number of events to generate" << std::setw(wp) << maxgen << " |" << std::endl
    << "| " << std::setw(wt) << "Events storage ? " << std::setw(wp) << store << " |" << std::endl
    << "| " << std::setw(wt) << "Debugging mode ? " << std::setw(wp) << debug << " |" << std::endl
    << "| " << std::setw(wt) << "Output file opened ? " << std::setw(wp) << (file!=(std::ofstream*)NULL && file->is_open()) << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Vegas integration parameters" << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Maximum number of iterations" << std::setw(wp) << itvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of function calls" << std::setw(wp) << ncvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of points to try per bin" << std::setw(wp) << npoints << " |" << std::endl
    << "| " << std::setw(wt) << "Integration smoothed (TREAT) ? " << std::setw(wp) << ntreat << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(wb) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Incoming particles " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Subprocess' mode" << std::setw(20) << pmode << std::setw(wp-20) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Incoming particles" << std::setw(5) << in1pdg << ", " << std::setw(5) << in2pdg << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Momenta [GeV/c]" << std::setw(5) << in1p << ", " << std::setw(5) << in2p << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Incoming photons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Virtuality in range [GeV^2]" << "[" << std::setw(4) << minq2 << ", " << std::setw(6) << maxq2 << "]" << std::setw(wp-14) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing leptons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Pair" << std::setw(2) << pair << " -> " << std::setw(wp-6) << particles << " |" << std::endl
    << "| " << std::setw(wt) << "Cuts mode" << std::setw(2) << mcut << " -> " << std::setw(wp-6) << cutsmode << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' pT in range [GeV/c]" << "[" << std::setw(4) << minpt << ", " << std::setw(4) << maxpt << "]" << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' energy in range [GeV]" << "[" << std::setw(4) << minenergy << ", " << std::setw(4) << maxenergy << "]" << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Pseudorapidity in range" << "[" << std::setw(4) << mineta << ", " << std::setw(4) << maxeta << "]" << std::setw(wp-12) << "" << " |" << std::endl
    //<< "| " << std::setw(wt) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw(3) << maxtheta << "]" << std::setw(wp-10) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing remnants " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl;
  if (this->hadroniser!=(Hadroniser*)NULL)
    os << "| " << std::setw(wt) << "Hadronisation algorithm" << std::setw(12) << hadroniser->GetName() << std::setw(wp-12) << "" << " |" << std::endl;
  os << "| " << std::setw(wt) << "Minimal mass [GeV/c^2]" << std::setw(wp) << minmx << " |" << std::endl
             << "| " << std::setw(wt) << "Maximal mass [GeV/c^2]" << std::setw(wp) << maxmx << " |" << std::endl
             << "|_" << std::right << std::setfill('_') << std::setw(wb) << "_|" << std::left;
  PrintInfo(os.str());
}

bool Parameters::ReadConfigFile(std::string inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_.c_str(), std::fstream::in);
  if (!f.is_open()) {
    return false;
  }

  PrintDebug(Form("File %s succesfully opened!", inFile_));
  std::ostringstream os;
  os << "======================================================" << std::endl
     << "Configuration file content : " << std::endl
     << "======================================================" << std::endl
     << std::left;

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
      os << std::setw(60) << " * Number of function calls" << this->ncvg << std::endl;
    }
    else if (key=="NCSG") {
      this->npoints = (int)atoi(value.c_str());
      os << std::setw(60) << " * Number of points to probe" << this->npoints << std::endl;
    }
    else if (key=="ITVG") {
      this->itvg = (int)atoi(value.c_str());
      os << std::setw(60) << " * Number of Vegas iterations" << this->itvg << std::endl;
    }
    else if (key=="INPP") {
      this->in1p = (double)atof(value.c_str());
      os << std::setw(60) << " * First incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
    }
    else if (key=="MODE") {
      this->process_mode = static_cast<Process::ProcessMode>(atoi(value.c_str()));
      os << std::setw(60) << " * Subprocess' mode" << this->process_mode << " --> ";
      switch (this->process_mode) {
        case Process::ElasticElastic: default: os << "elastic-elastic"; break;
        case Process::ElasticInelastic:        os << "elastic-inelastic"; break;
        case Process::InelasticElastic:        os << "inelastic-elastic"; break;
        case Process::InelasticInelastic:      os << "inelastic-inelastic"; break;
      }
      os << std::endl;
    }
    else if (key=="PMOD" or key=="EMOD") {
      this->remnant_mode = static_cast<Process::StructureFunctions>(atoi(value.c_str()));
      os << std::setw(60) << " * Outgoing primary particles' mode" << this->remnant_mode << " --> ";
      switch (this->remnant_mode) {
        case Process::Electron:        os << "electron"; break;
        case Process::SuriYennie:      os << "dissociating proton [SY structure functions]"; break;
        case Process::SuriYennieLowQ2: os << "dissociating proton [SY structure functions, for MX < 2 GeV, Q^2 < 5 GeV^2]"; break;
        case Process::SzczurekUleshchenko: os << "dissociating proton [SU structure functions]"; break;
        case Process::FioreVal:        os << "dissociating proton [parton model, only valence quarks]"; break;
        case Process::FioreSea:        os << "dissociating proton [parton model, only sea quarks]"; break;
        case Process::Fiore:           os << "dissociating proton [parton model, valence and sea quarks]"; break;
      }
      os << std::endl;
    }
    else if (key=="INPE") {
      this->in2p = (double)atof(value.c_str());
      os << std::setw(60) << " * Second incoming particles' momentum" << this->in1p << " GeV/c" << std::endl;
    }
    else if (key=="PAIR") {
      this->pair = (Particle::ParticleCode)atoi(value.c_str());
      os << std::setw(60) << " * Outgoing leptons' PDG id   " << this->pair << " --> ";
      switch (this->pair) {
        case Particle::Electron: default: os << "electrons"; break;
        case Particle::Muon:              os << "muons"; break;
        case Particle::Tau:               os << "taus"; break;
      }
      os << std::endl;
    }
    else if (key=="MCUT") {
      this->mcut = (int)atoi(value.c_str());
      os << std::setw(60) << " * Set of cuts to apply on the total generation  " << this->mcut << " --> ";
      switch (this->mcut) {
        case 0: default: os << "no cuts"; break;
        case 3:          os << "cuts on at least one outgoing lepton"; break;
        case 2:          os << "cuts on both the outgoing leptons"; break;
        case 1:          os << "Vermaseren's hypothetical detector cuts"; break;
      }
      os << std::endl;
    }
    else if (key=="PTCT") {
      this->minpt = (double)atof(value.c_str());
      os << std::setw(60) << " * Single outgoing lepton's minimal transverse momentum" << this->minpt << " GeV/c" << std::endl;
    }
    else if (key=="ECUT") {
      this->minenergy = (double)atof(value.c_str());
      os << std::setw(60) << " * Single outgoing lepton's minimal energy" << this->minenergy << " GeV" << std::endl;
    }
    else if (key=="NTRT") {
      this->ntreat = (double)atoi(value.c_str());
      os << std::setw(60) << " * Number of TREAT calls" << this->ntreat << std::endl;
    }
    else if (key=="NGEN") {
      this->maxgen = (double)atoi(value.c_str());
      os << std::setw(60) << " * Number of events to generate" << this->maxgen << std::endl;
    }
    else if (key=="THMN") {
      //this->mintheta = (double)atof(value.c_str());
      //this->SetThetaRange((double)atof(value.c_str()), 0.); // FIXME FIXME
      os << std::setw(60) << " * Minimal polar production angle for the leptons" << EtaToTheta(mineta) << std::endl;
    }
    else if (key=="THMX") {
      //this->maxtheta = (double)atof(value.c_str());
      //this->SetThetaRange(0., (double)atof(value.c_str())); //FIXME FIXME
      os << std::setw(60) << " * Maximal polar production angle for the leptons" << EtaToTheta(maxeta) << std::endl;
    }
    else if (key=="Q2MN") {
      this->minq2 = (double)atof(value.c_str());
      os << std::setw(60) << " * Minimal Q^2 for the incoming photons" << this->minq2 << " GeV^2" << std::endl;
    }
    else if (key=="Q2MX") {
      this->maxq2 = (double)atof(value.c_str());
      os << std::setw(60) << " * Maximal Q^2 for the incoming photons" << this->maxq2 << " GeV^2" << std::endl;
    }
    else if (key=="MXMN") {
      this->minmx = (double)atof(value.c_str());
      os << std::setw(60) << " * Minimal invariant mass of proton remnants" << this->minmx << " GeV/c^2" << std::endl;
    }
    else if (key=="MXMX") {
      this->maxmx = (double)atof(value.c_str());
      os << std::setw(60) << " * Maximal invariant mass of proton remnants" << this->maxmx << " GeV/c^2" << std::endl;
    }
    else if (key=="GPDF") {
      this->gpdf = (double)atoi(value.c_str());
      os << std::setw(60) << " * GPDF" << this->gpdf << std::endl;
    }
    else if (key=="SPDF") {
      this->spdf = (double)atoi(value.c_str());
      os << std::setw(60) << " * SPDF" << this->spdf << std::endl;
    }
    else if (key=="QPDF") {
      this->qpdf = (double)atoi(value.c_str());
      os << std::setw(60) << " * QPDF" << this->qpdf << std::endl;
    }
    else {
      PrintInfo(Form("<WARNING> Unrecognized argument : [%s] = %s", key, value));
    }
  }
  f.close();
  
  PrintInfo(os.str());
  
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

