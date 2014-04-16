#include <iostream>

#include "include/mcgen.h"

/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 */
int main(int argc, char* argv[]) {
  double xsec, err;
  MCGen mg;
  Event ev;
  //GamGamLL proc;
  GamPomVMLL proc;
  //Herwig6Hadroniser had;
  //Pythia6Hadroniser had;
  Jetset7Hadroniser had;

  if (argc==1) {
    std::cout << "[Main] [DEBUG] No config file provided. Setting the default parameters." << std::endl;
    mg.parameters->in1p = 4000.;
    mg.parameters->in2p = 4000.;
    mg.parameters->pair = 13;
    mg.parameters->p1mod = 2;
    mg.parameters->p2mod = 2;
    mg.parameters->mcut = 2;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 15.;
    //mg.parameters->SetEtaRange(-2.5, 2.5);
    //mg.parameters->ncvg = 5e3; //FIXME
    mg.parameters->generation = true;
    mg.parameters->maxgen = 2;
    //mg.parameters->maxgen = 1e5;
    mg.parameters->hadroniser = &had;
    mg.parameters->process = &proc;
  }
  else {
#ifdef DEBUG
    std::cout << "[Main] [DEBUG] Reading config file stored in " << argv[1] << std::endl;
#endif
    if (!mg.parameters->ReadConfigFile(std::string(argv[1]))) {
      std::cout << "=== Error reading the configuration !" << std::endl;
      std::cout << "  Please check your input file (" << std::string(argv[1]) << ")" << std::endl;
      return -1;
    }
  }

  mg.parameters->Dump();
  
  mg.ComputeXsection(&xsec, &err);
  if (mg.parameters->generation) {
    for (int i=0; i<mg.parameters->maxgen; i++) {
      ev = *mg.GenerateOneEvent();
      std::cout << ev.GetLHERecord();
    }
    //mg.LaunchGeneration();
  }
  
  //mg.parameters->StoreConfigFile("lastrun.card");

  return 0;
}

