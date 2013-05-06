#include <iostream>

#include "include/mcgen.h"
/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * Main caller for the present Monte Carlo generator. Loads the configuration
 * files' variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 */
int main(int argc, char* argv[]) {
  InputParameters ip;
  if (argc==1) {
    ip.in1p = 3500.;
    ip.in2p = 3500.;
    ip.pair = 13;
    ip.p1mod = 2;
    ip.p2mod = 2;
    ip.mcut = 2;
    ip.minpt = 5.;
    ip.debug = true;
  }
  else {
    std::cout << "Main : reading config file stored in " << argv[1] << std::endl;
    if (!ip.ReadConfigFile(std::string(argv[1]))) {
      return -1;
    }
  }
  MCGen mg(ip);
  mg.LaunchGen(1e4);
  //mg.AnalyzePhaseSpace("testing");
  return 0;
}

