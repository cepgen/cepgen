#include <iostream>

#include "include/mcgen.h"

using namespace std;

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
  //GamPomVMLL proc;
  //PPtoLL proc;
  ofstream output;
  //ofstream output2;
  
  Logger::GetInstance()->Level = Logger::Debug;
  //Logger::GetInstance()->OutputStream = ofstream("log.txt");
  
  if (argc==1) {
    Info("No config file provided. Setting the default parameters.");
    
    mg.parameters->hadroniser = new Pythia6Hadroniser;
    mg.parameters->process = new GamGamLL;
    mg.parameters->process_mode = Process::InelasticElastic;
    mg.parameters->remnant_mode = Process::SuriYennie;
    //mg.parameters->itvg = 2;
    
    mg.parameters->in1p = 4000.;
    mg.parameters->in2p = 4000.;
    mg.parameters->pair = Particle::Muon;
    mg.parameters->mcut = 2;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 5.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    //mg.parameters->ncvg = 5e3; //FIXME
    mg.parameters->generation = true;
    mg.parameters->maxgen = 2;
    //mg.parameters->maxgen = 1e5;
  }
  else {
    Debug(Form("Reading config file stored in %s", argv[1]));
    if (!mg.parameters->ReadConfigFile(string(argv[1]))) {
      Info(Form("Error reading the configuration!\n\t"
                "Please check your input file (%s)", argv[1]));
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
        cout << "Generating event #" << i+1 << endl;
      ev = *mg.GenerateOneEvent();
      ev.Dump();
    }

  }

  //mg.parameters->StoreConfigFile("lastrun.card");

  return 0;
}

