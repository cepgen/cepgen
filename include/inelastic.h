#ifndef _INELASTIC_H
#define _INELASTIC_H

#include <vector>
#include <string>

#include "particle.h"
#include "Pythia.h"

/*extern "C" {
  extern void pyinit_(const char*, const char*, const char*, double&, int, int, int);
  extern void pyevnt_();
  extern void pygive_(const char*, int);
  extern void pyfram_(int&);
  extern void pylist_(int&);
  extern void pystat_(int&);
  extern int pycomp_(int&);
  
  extern struct {
	  int n;
	  int npad;
	  int k[5][4000];
	  double p[5][4000];
	  double v[5][4000];
  } pyjets_;
  
  extern struct {
    int mdcy[3][500];
    int mdme[2][8000];
    double brat[8000];
    int kfpd[5][8000];
  } pydat3_;
}*/

class InelasticParticle {
  public:
    InelasticParticle();
    ~InelasticParticle();
    void PDF2PDG();
    void Hadronise();
  private:
    std::vector<Particle> *_shower;
    Pythia8::Pythia *_py;
    Pythia8::Event _ev;    
};

#endif
