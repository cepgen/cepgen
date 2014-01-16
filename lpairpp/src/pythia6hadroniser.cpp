#include "pythia6hadroniser.h"

Pythia6Hadroniser::Pythia6Hadroniser()
{
  _name = "Pythia6";
}

Pythia6Hadroniser::~Pythia6Hadroniser()
{
}

bool
Pythia6Hadroniser::Hadronise(Particle *part_)
{
  pyjets_.p[0][0] = part_->px;
  pyjets_.p[1][0] = part_->py;
  pyjets_.p[2][0] = part_->pz;
  pyjets_.p[3][0] = part_->E();
  pyjets_.p[4][0] = part_->M();

  pyjets_.k[0][0] = 1; // status
  pyjets_.k[1][0] = 2; // particle id
  pyjets_.k[2][0] = 0; // mother
  pyjets_.k[3][0] = 0; // daughter 1
  pyjets_.k[4][0] = 0; // daughter 2

  this->pyexec();
  //pyjets_.v[0][0] = 0;
  std::cout << "[Pythia6Hadroniser::Hadronise] INFO" << std::endl;
  //part_->Dump();
  return true;
}

bool
Pythia6Hadroniser::Hadronise(Event *ev_)
{
  int np;
  std::vector<Particle*> pl = ev_->GetParticles();
  std::vector<Particle*>::iterator p;
  std::ostringstream oss;
  //bool isprimary;
  int njoin, jlpsf[2];

  //std::cout << "[Pythia6Hadroniser::Hadronise] INFO" << std::endl;
  //ev_->Dump();

  // Muting the unneeded information
  //this->pygive("");

  // Filling the common block to propagate to PYTHIA6
  njoin = 0;
  pyjets_.n = pl.size();
  for (p=pl.begin(), np=0; p!=pl.end() && np<4000; p++, np++) {
    pyjets_.p[0][np] = (*p)->px;
    pyjets_.p[1][np] = (*p)->py;
    pyjets_.p[2][np] = (*p)->pz;
    pyjets_.p[3][np] = (*p)->E();
    pyjets_.p[4][np] = (*p)->M();

    pyjets_.k[0][np] = (*p)->status;
    pyjets_.k[1][np] = (*p)->pdgId;
    pyjets_.k[2][np] = 0; // mother
    pyjets_.k[3][np] = 0; // daughter 1
    pyjets_.k[4][np] = 0; // daughter 2
    
    for (int i=0; i<5; i++) {
      pyjets_.v[i][np] = 0.;
    }
    if ((*p)->status==3) {
      jlpsf[njoin] = np+1; //FIXME need to sort this vector<Particle*> !
      njoin++;
    }
  }

  if (njoin==0) return false;
  
#ifdef DEBUG
  std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] Joining " << njoin << " particle(s) in a same string" << std::endl;
  for (int i=0; i<njoin; i++) {
    std::cout << "--> " << jlpsf[i] << " (pdgId=" << pyjets_.k[1][jlpsf[i]-1] << ")" << std::endl;
  }
#endif

  this->pyjoin(njoin, jlpsf);
  //this->pylist(1);
  this->pyexec();
  //this->pylist(2);

  for (int p=0; p<pyjets_.n; p++) {
    /*isprimary = false;
    for (int i=0; i<njoin; i++) {
      if (p==jlpsf[i]-1) {
	isprimary = true;
	break;
      }
    }
    if (isprimary) continue;*/
    //if (p<2*njoin+1) continue;

    Particle pa(p+10, pyjets_.k[1][p]);
    pa.id = p;
    pa.role = p+10;
    pa.status = pyjets_.k[0][p];
    pa.pdgId = pyjets_.k[1][p];
    pa.P(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]);
    pa.M(pyjets_.p[4][p]);

    if (pyjets_.k[2][p]!=0) {
#ifdef DEBUG
      std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] " << "(pdgId=" << pa.pdgId << ") has mother (pdgId=" << pyjets_.k[1][pyjets_.k[2][p]-1] << ")" << std::endl;
#endif
      pa.SetMother(ev_->GetById(pyjets_.k[2][p]-1));
    }

    ev_->AddParticle(&pa);
  }
  return true;
}
