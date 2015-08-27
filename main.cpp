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
  GamGamLL proc;
  //GamPomVMLL proc;
  //PPtoLL proc;
  //Herwig6Hadroniser had;
  Pythia6Hadroniser had;
  //Jetset7Hadroniser had;//FIXME FIXME FIXME buggy !
  std::ofstream output;
  //std::ofstream output2;
  
  if (argc==1) {
    std::cout << "[Main] [DEBUG] No config file provided. Setting the default parameters." << std::endl;
    mg.parameters->in1p = 4000.;
    mg.parameters->in2p = 4000.;
    mg.parameters->pair = MUON;
    mg.parameters->remnant_mode = 11;
    mg.parameters->mcut = 2;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 5.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    //mg.parameters->ncvg = 5e3; //FIXME
    mg.parameters->generation = true;
    mg.parameters->maxgen = 2;
    //mg.parameters->maxgen = 1e5;
    mg.parameters->hadroniser = &had;
    mg.parameters->process = &proc;
    mg.parameters->process_mode = 2;
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

  // We might want to cross-check visually the validity of our run
  mg.parameters->Dump();
  
  // Let there be cross-section...
  mg.ComputeXsection(&xsec, &err);

  // The events generation starts here !

  Particles particles;
  Particles::const_iterator part;
  int status;

  if (mg.parameters->generation) {

    for (int i=0; i<mg.parameters->maxgen; i++) {
      if (i%10000==0)
        std::cout << "Generating event #" << i+1 << std::endl;
      ev = *mg.GenerateOneEvent();
      ev.Dump();
    }

  }

  //mg.parameters->StoreConfigFile("lastrun.card");

  return 0;
}

