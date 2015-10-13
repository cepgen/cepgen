#include "mcgen.h"

MCGen::MCGen() :
  veg(0), _xsec(-1.), _xsec_error(-1.)
{
  this->PrintHeader();

#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif

  srand(time(0)); // Random number initialization
  
  this->parameters = new Parameters;
}

MCGen::MCGen(Parameters *ip_) :
  parameters(ip_)
{}

MCGen::~MCGen()
{
  if (veg) delete veg;
  delete parameters;
#ifdef DEBUG
  std::cout << "[MCGen::~MCGen] [DEBUG] Destructor called" << std::endl;
#endif
}

void
MCGen::PrintHeader()
{
  const int bw = 64;
  const int lw = 43;
  int sp = (bw-lw)/2;

  std::cout << std::setw(bw+3) << std::setfill('-') << "" << std::endl;
  std::cout << std::setfill(' ') << std::left
	    << "| " << std::setw(sp) << "" << "             #                             " << std::setw(sp+1) << "" << "|" << std::endl  
	    << "| " << std::setw(sp) << "" << " ####        #       #####    ##   # ##### " << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(sp) << "" << "#    #       #       #    #  #  #  # #    #" << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(sp) << "" << "#      ##### #       #    # #    # # #    #" << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(sp) << "" << "#            #       #####  ###### # ##### " << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(sp) << "" << "#    #       #       #      #    # # #   # " << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(sp) << "" << " ####        ####### #      #    # # #    #" << std::setw(sp+1) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "Copyright (C) 2014  Laurent Forthomme" << "|" << std::endl
	    << "| " << std::setw(bw) << "                   <laurent.forthomme@uclouvain.be>" << "|" << std::endl
	    << "| " << std::setw(bw) << "              2005  Nicolas Schul" << "|" << std::endl
	    << "| " << std::setw(bw) << "              XXXX  Bryan (f.f in CDF version)" << "|" << std::endl
	    << "| " << std::setw(bw) << "         1991-1992  Olaf Duenger" << "|" << std::endl
	    << "| " << std::setw(bw) << "              199X  Dariusz Bocian" << "|" << std::endl
	    << "| " << std::setw(bw) << "              1996  MGVH (gmubeg.f in DESY version)" << "|" << std::endl
	    << "| " << std::setw(bw) << "              1994  ZEUS offline group" << "|" << std::endl
	    << "| " << std::setw(bw) << "              197X  Jos Vermaseren" << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "This program is free software: you can redistribute it and/or" << "|" << std::endl
	    << "| " << std::setw(bw) << "modify it under the terms of the GNU General Public License as" << "|" << std::endl
	    << "| " << std::setw(bw) << "published by the Free Software Foundation, either version 3 of" << "|" << std::endl
	    << "| " << std::setw(bw) << "the License, or any later version." << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "This program is distributed in the hope that it will be useful," << "|" << std::endl
	    << "| " << std::setw(bw) << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << "|" << std::endl
	    << "| " << std::setw(bw) << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << "|" << std::endl
	    << "| " << std::setw(bw) << "GNU General Public License for more details." << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl
	    << "| " << std::setw(bw) << "You should have received a copy of the GNU General Public" << "|" << std::endl
	    << "| " << std::setw(bw) << "License along with this program.  If not, see" << "|" << std::endl
	    << "| " << std::setw(bw) << "<http://www.gnu.org/licenses/>." << "|" << std::endl
	    << "| " << std::setw(bw) << "" << "|" << std::endl;
  std::cout << std::setw(bw+3) << std::setfill('-') << "" << std::endl;
}


void
MCGen::BuildVegas()
{
#ifdef DEBUG
  std::string topo;
  if (parameters->process_mode==1) {
    topo = "ELASTIC proton/proton";
  }
  else if (parameters->process_mode==2 or parameters->process_mode==3) {
    topo = "SINGLE-DISSOCIATIVE proton";
  }
  else if (parameters->process_mode==4) {
    topo = "DOUBLE-DISSOCIATIVE protons";
  }
  std::cout << "[MCGen::MCGen] [DEBUG] Considered topology : " << topo << " case" << std::endl;
#endif

  veg = new Vegas(parameters->process->GetNdim(parameters->process_mode), f, parameters);
}

void
MCGen::ComputeXsection(double* xsec_, double *err_)
{
  if (!veg) this->BuildVegas();

  std::cout << "[MCGen::ComputeXsection] Starting the computation of the process cross-section" << std::endl;
  veg->Integrate(xsec_, err_);
  this->_xsec = *xsec_;
  this->_xsec_error = *err_;
  std::cout << "[MCGen::ComputeXsection] Total cross-section = " << *xsec_ << " +/- " << *err_ << " pb" << std::endl;
  this->_xsec_comp = true;
}

Event*
MCGen::GenerateOneEvent()
{
  bool good = false;
  if (!this->_xsec_comp) {
    double xsec, err;
    this->ComputeXsection(&xsec, &err);
  }
  while (!good) {
    good = veg->GenerateOneEvent();
  }

  this->last_event = this->parameters->last_event;
  return (Event*)this->last_event;
}

void
MCGen::LaunchGeneration()
{
  // LHE file preparation
  if (!this->parameters->file->is_open()) {
    std::cerr << "[MCGen::LaunchGeneration] [ERROR] output file is not opened !" << std::endl;
  }
  //#ifdef DEBUG
  else {
    std::cout << "[MCGen::LaunchGeneration] [DEBUG] output file is correctly opened !" << std::endl;
  }
  //#endif
  *(this->parameters->file) << "<LesHouchesEvents version=\"1.0\">" << std::endl;
  *(this->parameters->file) << "<header>This file was created from the output of the CLPAIR generator</header>" << std::endl;
  *(this->parameters->file) << "<init>" << std::endl
	       << "2212 2212 "
	       << std::setprecision(2) << this->parameters->in1p << " "
	       << std::setprecision(2) << this->parameters->in2p << " "
	       << "0 0 10042 10042 2 1" << std::endl
	       << this->_xsec << " " << this->_xsec_error << " 0.26731120000E-03 0" << std::endl
	       << "</init>" << std::endl;

  veg->Generate();
  
  *(this->parameters->file) << "</LesHouchesEvents>" << std::endl;
}

double f(double* x_, size_t ndim_, void* params_)
{
  /*((Parameters*)params_)->process->SetPoint(ndim_, x_);
  return ((Parameters*)params_)->process->ComputeWeight(); //FIXME*/
  
  double ff;
  Parameters *p;
  Kinematics kin;
  Particle *in1, *in2;
  Timer tmr;
  bool hadronised;
  double num_hadr_trials;

  p = (Parameters*)params_;

  tmr.reset();

  //FIXME at some point introduce non head-on colliding beams ?

  ff = 0.;

#ifdef DEBUG
  std::cout << "=====================================" << std::endl
	    << "function f called ; some parameters :" << std::endl
            << "  pz(p1) = " << p->in1p << std::endl
            << "  pz(p2) = " << p->in2p << std::endl
            << "  remnant mode = " << p->remnant_mode << std::endl
	    << "=====================================" << std::endl;
#endif

  //FIXME electrons ?
  in1 = new Particle(1, p->in1pdg);
  in1->charge = p->in1pdg/abs(p->in1pdg);
  in1->P(0., 0.,  p->in1p);

  in2 = new Particle(2, p->in2pdg);
  in2->charge = p->in2pdg/abs(p->in2pdg);
  in2->P(0., 0., -p->in2p);

  p->process->SetIncomingParticles(*in1, *in2);

  kin.kinematics = p->process_mode;

  switch (kin.kinematics) {
    case 1:
      p->process->SetOutgoingParticles(3, PROTON, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, PROTON, 2); // Second outgoing proton
    case 2:
      p->process->SetOutgoingParticles(3, PROTON, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, QUARK_U, 2); // Second outgoing proton remnant
    case 3:
      p->process->SetOutgoingParticles(3, QUARK_U, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, PROTON, 2); // Second outgoing proton remnant
    case 4:
      p->process->SetOutgoingParticles(3, QUARK_U, 1); // First outgoing proton
      p->process->SetOutgoingParticles(5, QUARK_U, 2); // Second outgoing proton remnant
  }
  
  kin.q2min = p->minq2;
  kin.q2max = p->maxq2;
  //kin.q2min = 0; //FIXME
  //kin.q2max = -1; //FIXME
  kin.mode = p->mcut;
  kin.ptmin = p->minpt;
  kin.ptmax = p->maxpt;
  /*if (p->mintheta>=0. and p->maxtheta>=0.) { // FIXME need conversion from rapidity
    kin.thetamin = p->mintheta;
    kin.thetamax = p->maxtheta;
  }
  else {*/
    kin.etamin = p->mineta;
    kin.etamax = p->maxeta;
  //}
  kin.emin = p->minenergy;
  kin.emax = p->maxenergy;
  kin.mxmin = p->minmx;
  kin.mxmax = p->maxmx;
  
  p->process->SetPoint(ndim_, x_);
  p->process->SetKinematics(kin);
  p->process->SetOutgoingParticles(6, p->pair); // Outgoing leptons
  if (!p->process->IsKinematicsDefined()) {
    std::cout << "[f] [ERROR] Kinematics is not properly set" << std::endl;
    p->process->GetEvent()->Dump();
    return 0.;
  }
  ff = p->process->ComputeWeight();
  
  if (ff<0.) return 0.;
  
  if (p->store) { // MC events generation
    p->process->FillKinematics(false);
    p->process->GetEvent()->time_generation = tmr.elapsed();

    if (kin.kinematics>1) {
#ifdef DEBUG
      std::cout << "[f] [DEBUG] Event before calling the hadroniser (" << p->hadroniser->GetName() << ")" << std::endl;
      p->process->GetEvent()->Dump();
#endif
      num_hadr_trials = 0;
      do {
        hadronised = p->hadroniser->Hadronise(p->process->GetEvent());
#ifdef DEBUG
        if (num_hadr_trials>0) {
          std::cout << "[f] [DEBUG] Hadronisation failed. Trying for the " << num_hadr_trials+1 << "th time" << std::endl;
        }
#endif
        num_hadr_trials++;
      } while (!hadronised and num_hadr_trials<=p->hadroniser_max_trials);
      p->process->GetEvent()->num_hadronisation_trials = num_hadr_trials;
#ifdef DEBUG
      std::cout << "[f] [DEBUG] Event hadronisation succeded after " << p->process->GetEvent()->num_hadronisation_trials << " trial(s)" << std::endl;
#endif

      if (num_hadr_trials>p->hadroniser_max_trials) return 0.; //FIXME
#ifdef DEBUG
      std::cout << "[f] [DEBUG] Event after calling the hadroniser (" << p->hadroniser->GetName() << ")" << std::endl;
      p->process->GetEvent()->Dump();
#endif
    }
    p->process->GetEvent()->time_total = tmr.elapsed();
    
#ifdef DEBUG
  std::cout << "[f] [DEBUG]" << std::endl
	    << "       Generation time : " << std::setprecision(8) << p->process->GetEvent()->time_generation << " sec" << std::endl;
	    << "  Total (+ hadr.) time : " << std::setprecision(8) << p->process->GetEvent()->time_total << " sec" << std::endl;
#endif

    //*(p->file) << p->process->GetEvent()->GetLHERecord();
    *(p->last_event) = *(p->process->GetEvent());
    //p->process->GetEvent()->Store(p->file);

  }
  //p->process->ClearEvent();
  p->process->GetEvent()->clear(); // need to move this sw else ?

  delete in1;
  delete in2;

  return ff;
}

