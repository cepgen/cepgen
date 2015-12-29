#include <iostream>

#include "include/MCGen.h"

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
  
  //Logger::GetInstance()->Level = Logger::Debug;
  //Logger::GetInstance()->Level = Logger::DebugInsideLoop;
  //Logger::GetInstance()->OutputStream = ofstream("log.txt");
  
  if (argc==1) {
    Info("No config file provided. Setting the default parameters.");
    
    //mg.parameters->hadroniser = new Pythia6Hadroniser;
    mg.parameters->process = new GamGamLL;
    mg.parameters->process_mode = GenericProcess::InelasticElastic;
    mg.parameters->remnant_mode = GenericProcess::SuriYennie;
    //mg.parameters->itvg = 2;
    
    mg.parameters->in1p = 4000.;
    mg.parameters->in2p = 4000.;
    mg.parameters->pair = Particle::Muon;
    mg.parameters->mcut = 2;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 5.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    mg.parameters->ncvg = 5e4; //FIXME
    mg.parameters->generation = true;
    mg.parameters->maxgen = 2;
    //mg.parameters->maxgen = 1e5;
  }
  else {
    Debug(Form("Reading config file stored in %s", argv[1]));
    if (!mg.parameters->ReadConfigFile(argv[1])) {
      Info(Form("Error reading the configuration!\n\t"
                "Please check your input file (%s)", argv[1]));
      return -1;
    }
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->Dump();

  //for (unsigned int i=0; i<mg.GetNdim()/2; i++) { x[i*2] = 0.2; }
  //for (unsigned int i=0; i<mg.GetNdim(); i++) { x[i] = (i+1)*0.1; }
  //double x[9] = { 0.28950314, 0.13319218, 0.62257359, 0.01415460, 0.43977024, 0.60647596, 0.24439744, 0.40325279, 0.5 };
  //double x[9] = { 0.333, 0.336, 0.662, 0.069, 0.444, 0.403, 0.319, 0.321, 0.5};
  /*Logger::GetInstance()->Level = Logger::DebugInsideLoop;
  double x[9] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
  mg.ComputePoint(x);
  exit(0);*/
  
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

