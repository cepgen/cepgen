#include <iostream>
#include <fstream>

#include "../include/MCGen.h"
/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 */
int main(int argc, char* argv[]) {
  Parameters ip;
  double xsec, err, minpt, sqs;
  double min, max;
  std::ofstream tmp;
  int it;
  MCGen *mg;

  // max = 10.;
  //max = 10.;
  min = 0.;
  max = 50.;
  it = 100;
  //max = 14000.;
  /*min = 10.;
  max = 60.;
  it = 50;*/
  
  if (argc>1) it = atoi(argv[1]);

  ip.in1p = 3500.;
  ip.in2p = 3500.;
  ip.process = new GamGamLL;
  //ip.hadroniser = new Pythia6Hadroniser;
  //ip.process_mode = GenericProcess::ElasticElastic;
  //ip.process_mode = GenericProcess::InelasticElastic;
  ip.process_mode = GenericProcess::InelasticInelastic;
  ip.pair = Particle::Muon;
  ip.remnant_mode = GenericProcess::SuriYennie;
  //ip.SetThetaRange(5., 175.);
  //ip.SetThetaRange(0., 180.);
  ip.mineta =-2.5;
  ip.maxeta = 2.5;
  ip.maxmx = 1000.;
  ip.mcut = 2;
  ip.minenergy = 0.;
  // DEBUG
  //ip.itvg = 5;
  //ip.ncvg = 50000;
  //
  ip.minpt = 15.;
  ip.generation = false;
  ip.Dump();

  //tmp.open("tmp/xsec_lpairpp_elastic.dat");
  //tmp.open("tmp/xsec_lpairpp_singleinelastic.dat");
  tmp.open("tmp/xsec_lpairpp_doubleinelastic_v2.dat");
  //tmp.open("tmp/xsec_sqs_lpairpp_elastic_noetacut.dat");
  //tmp.open("tmp.dat");
  //tmp.open("tmp/xsec_sqs_lpairpp_singleinelastic_noetacut.dat");
  //tmp.open("tmp/xsec_sqs_lpairpp_doubleinelastic_noetacut.dat");
  mg = new MCGen(&ip);
  for (int i=0; i<=it; i++) {
    minpt = min+i/(double)it*(max-min);
    mg->parameters->minpt = minpt;
    //sqs = min+(double)i/(double)it*(max-min);
    //ip.in1p = sqs/2.;
    //ip.in2p = sqs/2.;
    mg->ComputeXsection(&xsec, &err);
    std::cout << minpt << "\t" << xsec << "\t" << err << std::endl;
    tmp << minpt << "\t" << xsec << "\t" << err << std::endl;
  }
  delete mg;
  return 0;
}

