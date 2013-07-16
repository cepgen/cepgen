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
  if (argc==1) {
    std::cout << "[Main] [DEBUG] No config file provided. Setting the default parameters." << std::endl;
    ip.in1p = 3500.;
    ip.in2p = 3500.;
    /*ip.in1p = 4000.;
    ip.in2p = 4000.;*/
    ip.pair = 13;
    ip.p1mod = 2;
    ip.p2mod = 2;
    ip.mcut = 2;
    ip.minpt = 0.;
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
  MCGen mg(ip);
  double xsec, err;
  mg.ComputeXsection(&xsec, &err);

  if (ip.generation) {
    //mg.LaunchGen(5e4);
    mg.Test();
  }
  //mg.Test();
  //mg.AnalyzePhaseSpace("testing/psprobe");*/
  return 0;
}

