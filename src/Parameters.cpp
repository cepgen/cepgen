#include "Parameters.h"

Parameters::Parameters() :
  process(0), process_mode(GenericProcess::ElasticElastic),
  remnant_mode(GenericProcess::SuriYennie),
  in1p(3500.), in2p(3500.),
  in1pdg(Particle::Proton), in2pdg(Particle::Proton),
  pair(Particle::Muon),
  mcut(2),
  minpt(0.), maxpt(-1.),
  minptdiff(0.), maxptdiff(-1.),
  minenergy(0.), maxenergy(-1.),
  mineta(-5.), maxeta(5.),
  minqt(0.), maxqt(500.),
  minq2(0.), maxq2(1.e5),
  minmx(1.07), maxmx(320.),
  ncvg(100000), itvg(10), npoints(100), first_run(true),
  generation(true), store(false), maxgen(1e5),
  symmetrise(true), ngen(0),
  gpdf(5), spdf(4), qpdf(12),
  hadroniser(0),
  hadroniser_max_trials(5)
{
  this->last_event = new Event();
  this->file = (std::ofstream*)NULL;
}

Parameters::~Parameters()
{
  delete last_event;
}

void Parameters::SetThetaRange(double thetamin_, double thetamax_)
{
  this->mineta = -log(tan(thetamax_/180.*Constants::Pi/2.));
  this->maxeta = -log(tan(thetamin_/180.*Constants::Pi/2.));

  Debug(Form("eta(min) = %5.2f => theta(min) = %5.2f"
             "eta(max) = %5.2f => theta(max) = %5.2f",
             mineta, thetamin_, maxeta, thetamax_));
}

void Parameters::Dump()
{
  std::string cutsmode, particles;
  std::ostringstream os;

  switch (mcut) {
    case 1:           cutsmode = "Vermaseren"; break;
    case 2:           cutsmode = "both leptons"; break;
    case 3:           cutsmode = "single lepton"; break;
    case 0:  default: cutsmode = "none"; break;
  }
  switch (pair) {
    case Particle::Electron:      particles = "electrons"; break;
    case Particle::Muon: default: particles = "muons"; break;
    case Particle::Tau:           particles = "taus"; break;
  }
  const int wb = 77;
  const int wt = 40;
  int wp = wb-wt-2;
  os 
    << std::left
    << std::endl
    << " _" << std::setfill('_') << std::setw(wb) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "_ " << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl;
  if (process)
    os << "| " << std::setw(wt) << "Process to generate"  << std::setw(wp) << process->GetName() << std::endl;
  os
    << "| " << std::setw(wt) << "Events generation ? " << std::setw(wp) << generation << " |" << std::endl
    << "| " << std::setw(wt) << "Number of events to generate" << std::setw(wp) << maxgen << " |" << std::endl
    << "| " << std::setw(wt) << "Events storage ? " << std::setw(wp) << store << " |" << std::endl
    << "| " << std::setw(wt) << "Debugging mode ? " << std::setw(wp) << (Logger::GetInstance()->Level) << " |" << std::endl
    << "| " << std::setw(wt) << "Output file opened ? " << std::setw(wp) << (file!=(std::ofstream*)NULL && file->is_open()) << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Vegas integration parameters" << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Maximum number of iterations" << std::setw(wp) << itvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of function calls" << std::setw(wp) << ncvg << " |" << std::endl
    << "| " << std::setw(wt) << "Number of points to try per bin" << std::setw(wp) << npoints << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|_" << std::setfill('_') << std::setw(wb) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill(' ') << "_|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Incoming particles " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Subprocess' mode" << std::setw(20) << process_mode << std::setw(wp-20) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Incoming particles" << std::setw(7) << in1pdg << ", " << std::setw(7) << in2pdg << std::setw(wp-16) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Momenta (GeV/c)" << std::setw(5) << in1p << ", " << std::setw(5) << in2p << std::setw(wp-12) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Incoming photons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Virtuality in range (GeV^2)" << "[" << std::setw(4) << minq2 << ", " << std::setw(6) << maxq2 << "]" << std::setw(wp-14) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing leptons " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "| " << std::setw(wt) << "Pair" << std::setw(2) << (int)pair << " -> " << std::setw(wp-6) << pair << " |" << std::endl
    << "| " << std::setw(wt) << "Cuts mode" << std::setw(2) << mcut << " -> " << std::setw(wp-6) << cutsmode << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' pT in range (GeV/c)" << std::right << "[" << std::setw(5) << minpt << ", " << std::setw(5) << maxpt << "]" << std::left << std::setw(wp-14) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Lepton(s)' energy in range (GeV)" << std::right << "[" << std::setw(5) << minenergy << ", " << std::setw(5) << maxenergy << "]" << std::left << std::setw(wp-14) << "" << " |" << std::endl
    << "| " << std::setw(wt) << "Pseudorapidity in range" << std::right << "[" << std::setw(5) << mineta << ", " << std::setw(5) << maxeta << "]" << std::left << std::setw(wp-14) << "" << " |" << std::endl
    //<< "| " << std::setw(wt) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw(3) << maxtheta << "]" << std::setw(wp-10) << "" << " |" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl
    << "|-" << std::setfill('-') << std::setw(wb-2) << " Outgoing remnants " << std::setfill(' ') << "-|" << std::endl
    << "| " << std::right << std::setw(wb) << " |" << std::left << std::endl;
  if (hadroniser)
    os << "| " << std::setw(wt) << "Hadronisation algorithm" << std::setw(12) << hadroniser->GetName() << std::setw(wp-12) << "" << " |" << std::endl;
  os << "| " << std::setw(wt) << "Minimal mass (GeV/c^2)" << std::setw(wp) << minmx << " |" << std::endl
             << "| " << std::setw(wt) << "Maximal mass (GeV/c^2)" << std::setw(wp) << maxmx << " |";
  Information(os.str());
}

bool Parameters::ReadConfigFile(const char* inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_, std::fstream::in);
  if (!f.is_open()) {
    return false;
  }

  Debug(Form("File '%s' succesfully opened!", inFile_));
  std::ostringstream os;
  os << "Configuration file content : " << "\n\t";

  while (f >> key >> value) {
    //std::cout << std::setw(60) << "[" << key << "] = " << value << std::endl;
    //if (strncmp(key.c_str(), "#")==0) continue; // FIXME need to ensure there is no extra space before !
    if (key[0]=='#') continue;
    if (key=="IEND") {
      int iend = (int)atoi(value.c_str());
      if (iend>1) {
        this->generation = true;
      }
    }
    else if (key=="DEBG") {
      Logger::GetInstance()->Level = static_cast<Logger::LoggingLevel>(atoi(value.c_str()));
    }
    else if (key=="NCVG") {
      this->ncvg = (int)atoi(value.c_str());
      os << " * Number of function calls: " << this->ncvg << "\n\t";
    }
    else if (key=="NCSG") {
      this->npoints = (int)atoi(value.c_str());
      os << " * Number of points to probe: " << this->npoints << "\n\t";
    }
    else if (key=="ITVG") {
      this->itvg = (int)atoi(value.c_str());
      os << " * Number of Vegas iterations: " << this->itvg << "\n\t";
    }
    else if (key=="INPP") {
      this->in1p = (double)atof(value.c_str());
      os << " * First incoming particles' momentum: " << this->in1p << " GeV/c\n\t";
    }
    else if (key=="PROC") {
      if (value=="lpair") {
      	this->process = new GamGamLL;
      	os << " * Process: LPAIR\n\t";
      }
      else if (value=="pptoll") {
        this->process = new PPtoLL;
        os << " * Process: PPTOLL\n\t";
      }
    }
    else if (key=="HADR") {
      if (value=="pythia6") {
        this->hadroniser = new Pythia6Hadroniser;
        os << " * Hadroniser: Pythia6\n\t";
      }
      if (value=="jetset7") {
        this->hadroniser = new Jetset7Hadroniser;
        os << " * Hadroniser: Jetset7\n\t";
      }
#ifdef PYTHIA8
      if (value=="pythia8") {
        this->hadroniser = new Pythia8Hadroniser;
        os << " * Hadroniser: Pythia8\n\t";
      }
#endif
    }
    else if (key=="MODE") {
      this->process_mode = static_cast<GenericProcess::ProcessMode>(atoi(value.c_str()));
      os << " * Subprocess' mode: " << (unsigned int)this->process_mode << " --> " << this->process_mode << "\n\t";
    }
    else if (key=="PMOD" or key=="EMOD") {
      this->remnant_mode = static_cast<GenericProcess::StructureFunctions>(atoi(value.c_str()));
      os << " * Outgoing primary particles' mode: " << (unsigned int)this->remnant_mode
      	 << " --> " << this->remnant_mode << "\n\t";
    }
    else if (key=="INPE") {
      this->in2p = (double)atof(value.c_str());
      os << " * Second incoming particles' momentum: " << this->in1p << " GeV/c\n\t";
    }
    else if (key=="PAIR") {
      this->pair = (Particle::ParticleCode)atoi(value.c_str());
      os << " * Outgoing leptons' PDG id: " << (int)this->pair
         << " --> " << this->pair << "\n\t";
    }
    else if (key=="MCUT") {
      this->mcut = (int)atoi(value.c_str());
      os << " * Set of cuts to apply on the total generation  : " << this->mcut
         << " --> ";
      switch (this->mcut) {
        case 0: default: os << "no cuts"; break;
        case 3:          os << "cuts on at least one outgoing lepton"; break;
        case 2:          os << "cuts on both the outgoing leptons"; break;
        case 1:          os << "Vermaseren's hypothetical detector cuts"; break;
      }
      os << "\n\t";
    }
    else if (key=="PTCT") {
      this->minpt = (double)atof(value.c_str());
      os << " * Single outgoing lepton's minimal transverse momentum: " << this->minpt << " GeV/c\n\t";
    }
    else if (key=="ECUT") {
      this->minenergy = (double)atof(value.c_str());
      os << " * Single outgoing lepton's minimal energy: " << this->minenergy << " GeV\n\t";
    }
    else if (key=="NGEN") {
      this->maxgen = (double)atoi(value.c_str());
      os << " * Number of events to generate: " << this->maxgen << "\n\t";
    }
    else if (key=="THMN") {
      //this->mintheta = (double)atof(value.c_str());
      //this->SetThetaRange((double)atof(value.c_str()), 0.); // FIXME FIXME
      os << " * Minimal polar production angle for the leptons" << EtaToTheta(mineta) << "\n\t";
    }
    else if (key=="THMX") {
      //this->maxtheta = (double)atof(value.c_str());
      //this->SetThetaRange(0., (double)atof(value.c_str())); //FIXME FIXME
      os << " * Maximal polar production angle for the leptons" << EtaToTheta(maxeta) << "\n\t";
    }
    else if (key=="ETMN") {
      this->mineta = (double)atof(value.c_str());
      os << " * Minimal pseudo-rapidity for the leptons: " << this->mineta << "\n\t";
    }
    else if (key=="ETMX") {
      this->maxeta = (double)atof(value.c_str());
      os << " * Maximal pseudo-rapidity for the leptons: " << this->maxeta << "\n\t";
    }
    else if (key=="Q2MN") {
      this->minq2 = (double)atof(value.c_str());
      os << " * Minimal Q^2 for the incoming photons: " << this->minq2 << " GeV^2\n\t";
    }
    else if (key=="Q2MX") {
      this->maxq2 = (double)atof(value.c_str());
      os << " * Maximal Q^2 for the incoming photons: " << this->maxq2 << " GeV^2\n\t";
    }
    else if (key=="MXMN") {
      this->minmx = (double)atof(value.c_str());
      os << " * Minimal invariant mass of proton remnants: " << this->minmx << " GeV/c^2\n\t";
    }
    else if (key=="MXMX") {
      this->maxmx = (double)atof(value.c_str());
      os << " * Maximal invariant mass of proton remnants: " << this->maxmx << " GeV/c^2\n\t";
    }
    else if (key=="GPDF") {
      this->gpdf = (double)atoi(value.c_str());
      os << " * GPDF: " << this->gpdf << "\n\t";
    }
    else if (key=="SPDF") {
      this->spdf = (double)atoi(value.c_str());
      os << " * SPDF: " << this->spdf << "\n\t";
    }
    else if (key=="QPDF") {
      this->qpdf = (double)atoi(value.c_str());
      os << " * QPDF: " << this->qpdf << "\n\t";
    }
    else {
      Information(Form("<WARNING> Unrecognized argument : [%s] = %s", key.c_str(), value.c_str()));
    }
  }
  f.close();
  
  Information(os.str());
  
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
  // ...
  f.close();
  return true;
}

