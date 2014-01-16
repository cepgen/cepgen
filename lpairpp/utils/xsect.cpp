#include <iostream>
#include <fstream>

#include "../include/mcgen.h"
/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 */
int main(int argc, char* argv[]) {
  InputParameters ip;
  double xsec, err, minpt;
  double max;
  std::ofstream tmp;
  int it;

  // max = 10.;
  max = 10.;

  if (argc>1) {
    it = atoi(argv[1]);
  }
  else {
    it = 100;
  }

  ip.in1p = 3500.;
  ip.in2p = 3500.;
  ip.pair = 13;
  ip.p1mod = 2;
  //ip.p1mod = 11;
  ip.p2mod = 2;
  ip.mcut = 2;
  ip.minenergy = 0.;
  //ip.itmx = 5;
  ip.generation = false;
  std::cout << "test" << std::endl;
  ip.Dump();

  tmp.open("tmp/xsec_lpairpp_elastic.dat");
  //tmp.open("tmp/xsec_lpairpp_singleinelastic.dat");

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

