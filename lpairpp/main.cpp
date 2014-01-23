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
  InputParameters ip;
  double xsec, err;
  std::ofstream of, fd;
  //HEPRUP hr;
  EventsList *ev;

  ev = new EventsList(&of, 1000);

  if (argc==1) {
    std::cout << "[Main] [DEBUG] No config file provided. Setting the default parameters." << std::endl;
    ip.in1p = 3500.;
    ip.in2p = 3500.;
    ip.pair = 13;
    ip.p1mod = 11;
    ip.p2mod = 2;
    ip.mcut = 2;
    ip.minenergy = 0.; //FIXME
    ip.minpt = 5.;
    ip.maxgen = 1e1;
    ip.itvg = 1;
    ip.ncvg = 5e3; //FIXME
    //ip.maxgen = 1e5;
    //ip.SetEtaRange(-2.5, 2.5);
  }
  else {
#ifdef DEBUG
    std::cout << "[Main] [DEBUG] Reading config file stored in " << argv[1] << std::endl;
#endif
    if (!ip.ReadConfigFile(std::string(argv[1]))) {
      std::cout << "=== Error reading the configuration !" << std::endl;
      std::cout << "  Please check your input file (" << std::string(argv[1]) << ")" << std::endl;
      return -1;
    }
  }

  //hrp.SetParameters(ip);

  ip.generation = true;
  of.open("test");
  fd.open("test_q2");
  ip.file = &of;
  ip.file_debug = &fd;
  ip.Dump();

  MCGen mg(ip);

  mg.ComputeXsection(&xsec, &err);
  if (ip.generation) {
    //ev->Info();
    //ip.eventslist = ev;
    ip.eventslist = new EventsList(&of, 1000);
    mg.LaunchGeneration();
    ip.eventslist->Info();
    //std::cout << "--> " << ip.eventslist->NumEvents() << std::endl;
  }
  
  ip.StoreConfigFile("lastrun.card");

  delete ev;
  of.close();
  fd.close();
  return 0;
}

