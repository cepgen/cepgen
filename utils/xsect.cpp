#include <iostream>
#include <fstream>

#include "mcgen.h"
/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 */
int main(int argc, char* argv[]) {
  InputParameters ip;
  double xsec, err, minpt;
  double max;
  std::ofstream tmp;
  int it;

  max = 10.;

  if (argc>1) {
    it = atoi(argv[1]);
  }
  else {
    it = 10;
  }

  ip.in1p = 3500.;
  ip.in2p = 3500.;
  ip.pair = 13;
  ip.p1mod = 2;
  ip.p2mod = 2;
  ip.mcut = 2;

  tmp.open("tmp/xsec.dat");
  for (int i=0; i<it; i++) {
    minpt = (double)i/(double)it*max;
    ip.minpt = minpt;
    MCGen *mg = new MCGen(ip);
    mg->ComputeXsection(&xsec, &err);
    std::cout << minpt << "\t" << xsec << "\t" << err << std::endl;
    tmp << minpt << "\t" << xsec << "\t" << err << std::endl;
    delete mg;
  }
  return 0;
}

