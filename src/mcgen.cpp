#include "mcgen.h"

MCGen::MCGen() :
  fVegas(0), fCrossSection(-1.), fCrossSectionError(-1.)
{
  Debug("Generator initialized");
  
  try { this->PrintHeader(); }  catch (Exception& e) { e.Dump(); }
  
  srand(time(0)); // Random number initialization
  
  this->parameters = new Parameters;
}

MCGen::MCGen(Parameters *ip_) :
  parameters(ip_), fVegas(0)
{}

MCGen::~MCGen()
{
  Debug("Destructor called");
  
  if (fVegas) delete fVegas;
  delete parameters;
}

void
MCGen::PrintHeader()
{
  std::string tmp;
  std::ostringstream os;
  std::ifstream hf("README");
  if (!hf.good()) throw Exception(__PRETTY_FUNCTION__, "Failed to open README file", JustWarning);
  while (true) {
    if (!hf.good()) break;
    getline(hf, tmp);
    os << "\n " << tmp;
  }
  hf.close();
  Info(os.str().c_str());
}


void
MCGen::BuildVegas()
{
  if (Logger::GetInstance()->Level>=Logger::Debug) {
    std::string topo;
    switch (parameters->process_mode) {
      case Process::ElasticElastic:
        topo = "elastic proton/proton"; break;
      case Process::ElasticInelastic:
      case Process::InelasticElastic:
        topo = "single-dissociative proton"; break;
      case Process::InelasticInelastic:
        topo = "double-dissociative protons"; break;
    }
    Debug(Form("Considered topology: %s case", topo.c_str()));
  }
  
  fVegas = new Vegas(GetNdim(), f, parameters);
}

void
MCGen::ComputeXsection(double* xsec_, double *err_)
{
  if (!fVegas) BuildVegas();

  Info("Starting the computation of the process cross-section");
  
  fVegas->Integrate(xsec_, err_);
  
  fCrossSection = *xsec_;
  fCrossSectionError = *err_;
  fHasCrossSection = true;
  
  Info(Form("Total cross section: %f +/- %f pb", *xsec_, *err_));
}

Event*
MCGen::GenerateOneEvent()
{
  bool good = false;
  if (!fHasCrossSection) {
    double xsec, err;
    ComputeXsection(&xsec, &err);
  }
  while (!good) {
    good = fVegas->GenerateOneEvent();
  }

  last_event = this->parameters->last_event;
  return (Event*)this->last_event;
}

void
MCGen::LaunchGeneration()
{
  // LHE file preparation
  if (!this->parameters->file->is_open()) { Debug("Output file is not opened"); }
  else                                    { Debug("Output file successfully opened"); }
  
  *(this->parameters->file) << "<LesHouchesEvents version=\"1.0\">" << std::endl;
  *(this->parameters->file) << "<header>This file was created from the output of the CLPAIR generator</header>" << std::endl;
  *(this->parameters->file) << "<init>" << std::endl
	       << "2212 2212 "
	       << std::setprecision(2) << this->parameters->in1p << " "
	       << std::setprecision(2) << this->parameters->in2p << " "
	       << "0 0 10042 10042 2 1" << std::endl
	       << fCrossSection << " " << fCrossSectionError << " 0.26731120000E-03 0" << std::endl
	       << "</init>" << std::endl;

  fVegas->Generate();
  
  *(this->parameters->file) << "</LesHouchesEvents>" << std::endl;
}

double f(double* x_, size_t ndim_, void* params_)
{
  double ff;
  Parameters *p;
  Kinematics kin;
  Particle *in1, *in2;
  Timer tmr;
  bool hadronised;
  double num_hadr_trials;
  std::ostringstream os;

  p = static_cast<Parameters*>(params_);

  //if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) {
  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) {
    os.str(""); for (unsigned int i=0; i<ndim_; i++) { os << x_[i] << " "; }
    DebugInsideLoop(Form("Computing dim-%d point ( %s)", ndim_, os.str().c_str()));
  }

  tmr.reset();

  //FIXME at some point introduce non head-on colliding beams ?

  ff = 0.;

  DebugInsideLoop(Form("Function f called -- some parameters:\n\t"
                            "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                            "  remnant mode: %d",
                            p->in1p, p->in2p, p->remnant_mode));
  
  kin.kinematics = p->process_mode;
  kin.q2min = p->minq2;
  kin.q2max = p->maxq2;
  kin.mode = p->mcut;
  kin.ptmin = p->minpt;
  kin.ptmax = p->maxpt;
  kin.etamin = p->mineta;
  kin.etamax = p->maxeta;
  kin.emin = p->minenergy;
  kin.emax = p->maxenergy;
  kin.mxmin = p->minmx;
  kin.mxmax = p->maxmx;
  
  p->process->SetKinematics(kin);
  p->process->SetPoint(ndim_, x_);
  p->process->GetEvent()->clear(); // need to move this sw else ?

  // Start by adding primary particles
  //FIXME electrons ?
  in1 = new Particle(1, p->in1pdg);
  in1->charge = p->in1pdg/abs(p->in1pdg);
  in1->P(0., 0.,  p->in1p);

  in2 = new Particle(2, p->in2pdg);
  in2->charge = p->in2pdg/abs(p->in2pdg);
  in2->P(0., 0., -p->in2p);  
  
  p->process->SetIncomingParticles(*in1, *in2);

  // Then add outgoing leptons
  p->process->SetOutgoingParticles(6, p->pair);
  p->process->SetOutgoingParticles(7, p->pair);

  // Then add outgoing protons or remnants
  switch (kin.kinematics) {
    case Process::ElasticElastic:
      p->process->SetOutgoingParticles(3, Particle::Proton, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, Particle::Proton, 2); // Second outgoing proton
      break;
    case Process::ElasticInelastic:
      //p->process->SetOutgoingParticles(3, Particle::Proton, 1); // First outgoing proton
      //p->process->SetOutgoingParticles(5, Particle::uQuark, 2); // Second outgoing proton remnant
      break;
    case Process::InelasticElastic:
      p->process->SetOutgoingParticles(3, Particle::uQuark, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, Particle::Proton, 2); // Second outgoing proton remnant
      break;
    case Process::InelasticInelastic:
      p->process->SetOutgoingParticles(3, Particle::uQuark, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, Particle::uQuark, 2); // Second outgoing proton remnant
      break;
  }

  // Check that everything is there
  if (!p->process->IsKinematicsDefined()) return 0.;
  
  try { ff = p->process->ComputeWeight(); } catch (Exception& e) { e.Dump(); }
  
  if (ff<0.) return 0.;
  
  if (p->store) { // MC events generation
    p->process->FillKinematics(false);
    p->process->GetEvent()->time_generation = tmr.elapsed();

    if (kin.kinematics!=Process::ElasticElastic) {

      Debug(Form("Event before calling the hadroniser (%s)", p->hadroniser->GetName().c_str()));
      if (Logger::GetInstance()->Level>=Logger::Debug) p->process->GetEvent()->Dump();
      
      num_hadr_trials = 0;
      do {
        try { hadronised = p->hadroniser->Hadronise(p->process->GetEvent()); } catch (Exception& e) { e.Dump(); }

        if (num_hadr_trials>0) { Debug(Form("Hadronisation failed. Trying for the %dth time", num_hadr_trials+1)); }
        
        num_hadr_trials++;
      } while (!hadronised and num_hadr_trials<=p->hadroniser_max_trials);
      p->process->GetEvent()->num_hadronisation_trials = num_hadr_trials;

      Debug(Form("Event hadronisation succeeded after %d trial(s)", p->process->GetEvent()->num_hadronisation_trials));

      if (num_hadr_trials>p->hadroniser_max_trials) return 0.; //FIXME
      
      Debug(Form("Event after calling the hadroniser (%s)", p->hadroniser->GetName().c_str()));
      if (Logger::GetInstance()->Level>=Logger::Debug) p->process->GetEvent()->Dump();
    }
    p->process->GetEvent()->time_total = tmr.elapsed();
    
    Debug(Form("Generation time:       %5.6f sec\n\t"
               "Total time (gen+hadr): %5.6f sec",
               p->process->GetEvent()->time_generation,
               p->process->GetEvent()->time_total));

    *(p->last_event) = *(p->process->GetEvent());
    //p->process->GetEvent()->Store(p->file);

  }
  //p->process->ClearEvent();
  p->process->GetEvent()->clear(); // need to move this sw else ?

  delete in1;
  delete in2;

  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) {
    os.str(""); for (unsigned int i=0; i<ndim_; i++) { os << Form("%10.8f ", x_[i]); }
    Debug(Form("f value for  dim-%d point ( %s): %4.4f", ndim_, os.str().c_str(), ff));
  }
  
  return ff;
}

