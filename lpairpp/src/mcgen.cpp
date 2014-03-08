#include "mcgen.h"

MCGen::MCGen() :
  _xsec(-1.), _xsec_error(-1.)
{
  this->PrintHeader();

#ifdef DEBUG
  std::cout << "[MCGen::MCGen] [DEBUG] MCGen initialized !" << std::endl;
#endif

  srand(time(0));
  
  this->parameters = new Parameters;
}

MCGen::MCGen(Parameters *ip_) 
{
  this->parameters = ip_;
  this->BuildVegas();
}

MCGen::~MCGen()
{
  delete veg;
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
            << "| " << "Version "<< std::setw(bw-8) << SVN_REV << "|" << std::endl
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
  unsigned int ndim;
  std::string topo;

  if (this->parameters->p1mod<=2 && this->parameters->p2mod<=2) {
    topo = "ELASTIC proton/proton";
    ndim = 7;
  }
  else if (this->parameters->p1mod<=2 || this->parameters->p2mod<=2) {
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
  std::cout << "[MCGen::MCGen] [DEBUG] Cuts mode : " << this->parameters->mcut << std::endl;
  switch(this->parameters->mcut) {
    case 1:
    case 2:
      std::cout << "[MCGen::MCGen] [DEBUG] Single leptons' transverse momentum condition : ";
      if (this->parameters->minpt<=0.) {
        std::cout << "no pT cut" << std::endl;
        break;
      }
      if (this->parameters->maxpt>0.) {
        std::cout << "pT in range [" 
                  << this->parameters->minpt << " GeV/c, " 
                  << this->parameters->maxpt << " GeV/c]" 
                  << std::endl;
        break;
      }
      std::cout << "pT > " << this->parameters->minpt << " GeV/c";
      if (this->parameters->mcut==1) {
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
  veg = new Vegas(ndim,f, this->parameters);
}

void
MCGen::ComputeXsection(double* xsec_, double *err_)
{
  this->BuildVegas();
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
  float time;
  std::clock_t c_start, c_stop;

  c_start = std::clock();

  if (!this->_xsec_comp) {
    double xsec, err;
    this->ComputeXsection(&xsec, &err);
  }

  while (!good) {
    good = veg->GenerateOneEvent();
  }
  c_stop = std::clock();
#ifdef DEBUG
  //if (c_stop!=c_start)
  std::cout << "[MCGen::GenerateOneEvent] [DEBUG]" << std::endl
	    << "  Generation time (CPU time) : " << std::setprecision(8) << (c_stop-c_start) << " cycles" << std::endl;
#endif
  time = (c_stop-c_start)/CLOCKS_PER_SEC*1000.;
  this->last_event = this->parameters->last_event;
  this->last_event->time_cpu = time;
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

double f(double* x_, size_t ndim_, void* params_) {
  double ff;
  int outp1pdg, outp2pdg;
  Parameters *p;
  Kinematics kin;
  Particle *in1, *in2;

  p = (Parameters*)params_;

  //FIXME at some point introduce non head-on colliding beams ?

  ff = 0.;

#ifdef DEBUG
  std::cout << "=====================================" << std::endl
	    << "function f called ; some parameters :" << std::endl
            << "  pz(p1) = " << p->in1p << std::endl
            << "  pz(p2) = " << p->in2p << std::endl
            << "   f(p1) = " << p->p1mod << std::endl
            << "   f(p2) = " << p->p2mod << std::endl
	    << "=====================================" << std::endl;
#endif


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

  p->process->GetEvent()->clear();
  p->process->SetPoint(ndim_, x_);
  p->process->SetKinematics(kin);
  p->process->SetIncomingParticles(*in1, *in2);
  p->process->SetOutgoingParticles(3, outp1pdg); // First outgoing proton
  p->process->SetOutgoingParticles(5, outp2pdg); // Second outgoing proton
  p->process->SetOutgoingParticles(6, p->pair); // Outgoing leptons
  if (!p->process->IsKinematicsDefined()) {
    std::cout << "[f] [ERROR] Kinematics is not properly set" << std::endl;
    return 0.;
  }
  ff = p->process->ComputeWeight();
  
  if (ff<0.) {
    return 0.;
  }
  if (p->store) { // MC events generation
    p->process->FillKinematics(false);
    if (kin.kinematics>1) {
#ifdef DEBUG
      std::cout << "[f] [DEBUG] Event before calling the hadroniser (" << p->hadroniser->GetName() << ")" << std::endl;
      p->process->GetEvent()->Dump();
#endif
      p->hadroniser->Hadronise(p->process->GetEvent());
#ifdef DEBUG
      std::cout << "[f] [DEBUG] Event after calling the hadroniser (" << p->hadroniser->GetName() << ")" << std::endl;
      p->process->GetEvent()->Dump();
#endif
    }
    
    //*(p->file) << p->process->GetEvent()->GetLHERecord();
    *(p->last_event) = *(p->process->GetEvent());
    //p->process->GetEvent()->Store(p->file);

  }

  delete in1;
  delete in2;

  return ff;
}

