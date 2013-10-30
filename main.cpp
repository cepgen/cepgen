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
  if (argc==1) {
    std::cout << "[Main] [DEBUG] No config file provided. Setting the default parameters." << std::endl;
    ip.in1p = 3500.;
    ip.in2p = -3500.;
    /*ip.in1p = 4000.;
    ip.in2p = -4000.;*/
    ip.pair = 13;
    ip.p1mod = 2;
    ip.p2mod = 2;
    ip.mcut = 2;
    ip.minenergy = 1.;
    ip.minpt = 5.;
    ip.maxgen = 1e5;
    //ip.debug = true;
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
  std::srand(std::time(0));
  /*std::cout << (double)rand()/(double)RAND_MAX << std::endl;
  std::cout << (double)rand()/(double)RAND_MAX << std::endl;
  std::cout << (double)rand()/(double)RAND_MAX << std::endl;*/
  ip.generation = true;
  std::ofstream of("test");
  std::ofstream fd("test_q2");
  ip.file = &of;
  ip.file_debug = &fd;
  /*ip.minenergy = 0.;
  ip.mintheta = 5.;
  ip.maxtheta = 175.;*/
  ip.Dump();
  MCGen mg(ip);
  /*std::ofstream out;
  out.open("haha");
  for (int i=0; i<50000; i++) {
    mg.ComputeXsection(&xsec, &err);
    //if (i%100==0) {
      std::cout << "test --> " << i << "\t" << xsec << std::endl;
    //}
    out << xsec << "\t" << err << std::endl;
  }
  out.close();*/
  std::cout << "before computing the cross-section" << std::endl;
  mg.ComputeXsection(&xsec, &err);
  std::cout << " --> xsect = " << xsec << " +/- " << err << std::endl;

  if (ip.generation) {
    mg.LaunchGen();
  }
  
  ip.StoreConfigFile("lastrun.card");
  //mg.Test();
  //mg.AnalyzePhaseSpace("testing/psprobe");*/
  return 0;
}

