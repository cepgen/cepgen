#include "inelastic.h"

InelasticParticle::InelasticParticle()
{
  this->_shower = new std::vector<Particle>();
  _py->readString("Beams:eCM = 8000.");
  _py->init("include/pythia8175/xmldoc/Index.xml");
  _py->init(2212, 2212, 7000.);
  _ev = _py->event;
}

InelasticParticle::~InelasticParticle()
{
  delete this->_shower;
}

void InelasticParticle::PDF2PDG()
{
  /*double dens[8], x, qscale;*/
  
}

void InelasticParticle::Hadronise()
{
  /*//char opt[32];
  std::stringstream s;
  s << "parp(171)=" << 10.;
  //Pythia6Interface::pygive(s.str().c_str());*/
  
  //_py->readString("SoftQCD::singleDiffractive = on");
  _py->readString("Diffraction::Pomflux = 5"); // MBR, as implemented by Robert Ciesielski and Konstantin Goulianos
}
