#include "inelastic.h"

InelasticParticle::InelasticParticle()
{
  //this->_shower = new std::vector<Particle>();
  //_py->readString("Beams:eCM = 8000.");
  //_py->init("include/pythia8175/xmldoc/Index.xml");
  //_py->init(2212, 2212, 7000.);
  //_had = Pythia8::HadronLevel();
  //_ev = _py->event;
}

InelasticParticle::~InelasticParticle()
{
  delete this->_shower;
}

void
InelasticParticle::PDF2PDG()
{
  /*double dens[8], x, qscale;*/
  
}

bool
InelasticParticle::Hadronise(std::string algo_)
{
  if (algo_=="") {
    return false;
  }
  if (algo_=="pythia6") {
    std::cout << "hadronisation using pythia6" << std::endl;
    Pythia6Hadroniser py;
    py.Hadronise(this);
  }
  /*//char opt[32];
  std::stringstream s;
  s << "parp(171)=" << 10.;
  //Pythia6Interface::pygive(s.str().c_str());*/
  
  //_py->readString("SoftQCD::singleDiffractive = on");
  //_py->readString("Diffraction::Pomflux = 5"); // MBR, as implemented by Robert Ciesielski and Konstantin Goulianos
  return true;
}
