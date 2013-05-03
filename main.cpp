#include <iostream>

#include "include/mcgen.h"
/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 */
int main(int argc, char* argv[]) {
  InputParameters ip;
  ip.in1p = 3500.;
  ip.in2p = 3500.;
  /*ip.in1p = 270.;
  ip.in2p = 270.;*/
  ip.pair = 13;
  ip.p1mod = 2;
  ip.p2mod = 2;
  ip.mcut = 2;
  ip.minpt = 5.;
  ip.plot->SetOutputFile("haha");
  MCGen mg(ip);
  mg.LaunchGen(1e8);
  //mg.AnalyzePhaseSpace("testing");
}

