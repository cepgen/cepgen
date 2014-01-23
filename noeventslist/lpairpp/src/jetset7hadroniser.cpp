#include "jetset7hadroniser.h"

Jetset7Hadroniser::Jetset7Hadroniser()
{
  _name = "Jetset7";
}

Jetset7Hadroniser::~Jetset7Hadroniser()
{
}

bool
Jetset7Hadroniser::Hadronise(Particle *part_)
{
  lujets_.p[0][0] = part_->px;
  lujets_.p[1][0] = part_->py;
  lujets_.p[2][0] = part_->pz;
  lujets_.p[3][0] = part_->E();
  lujets_.p[4][0] = part_->M();

  lujets_.k[0][0] = 1; // status
  lujets_.k[1][0] = 2; // particle id
  lujets_.k[2][0] = 0; // mother
  lujets_.k[3][0] = 0; // daughter 1
  lujets_.k[4][0] = 0; // daughter 2

  this->luexec();
  std::cout << "[Jetset7Hadroniser::Hadronise] INFO" << std::endl;
  return true;
}

bool
Jetset7Hadroniser::Hadronise(Event *ev_)
{
  int np;
  std::vector<Particle*> pl = ev_->GetParticles();
  std::vector<Particle*>::iterator p;
  std::ostringstream oss;
  //bool isprimary;
  int njoin, jlpsf[2];

  //std::cout << "[Jetset7Hadroniser::Hadronise] INFO" << std::endl;
  //ev_->Dump();

  // Filling the common block to propagate to JETSET7
  njoin = 0;
  lujets_.n = pl.size();
  for (p=pl.begin(); p!=pl.end(); p++) {
    np = (*p)->id;
    lujets_.p[0][np] = (*p)->px;
    lujets_.p[1][np] = (*p)->py;
    lujets_.p[2][np] = (*p)->pz;
    lujets_.p[3][np] = (*p)->E();
    lujets_.p[4][np] = (*p)->M();

    lujets_.k[0][np] = (*p)->status;
    lujets_.k[1][np] = (*p)->pdgId;
    lujets_.k[2][np] = 0; // mother
    lujets_.k[3][np] = 0; // daughter 1
    lujets_.k[4][np] = 0; // daughter 2
    
    std::cout << "---> " << this->luname((*p)->pdgId) << "\t" << np << std::endl;
    
    /*for (int i=0; i<5; i++) {
      lujets_.v[i][np] = 0.;
      }*/
    if ((*p)->status==3) {
      jlpsf[njoin] = np+1; //FIXME need to sort this vector<Particle*> !
      njoin++;
    }
  }

  if (njoin==0) return false;
  
  //#ifdef DEBUG
  std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] Joining " << njoin << " particle(s) in a same string" << std::endl;
  for (int i=0; i<njoin; i++) {
    std::cout << "--> " << jlpsf[i] << " (pdgId=" << lujets_.k[1][jlpsf[i]-1] << ")" << std::endl;
  }
  //#endif

  this->lujoin(njoin, jlpsf);
  //this->lulist(1);
  this->luexec();
  //this->lulist(2);

  for (int p=0; p<lujets_.n; p++) {
    /*isprimary = false;
    for (int i=0; i<njoin; i++) {
      if (p==jlpsf[i]-1) {
	isprimary = true;
	break;
      }
    }
    if (isprimary) continue;*/
    //if (p<2*njoin+1) continue;

    Particle pa(p+10, lujets_.k[1][p]);
    pa.id = p;
    pa.role = p+10;
    pa.status = lujets_.k[0][p];
    pa.pdgId = lujets_.k[1][p];
    pa.P(lujets_.p[0][p], lujets_.p[1][p], lujets_.p[2][p], lujets_.p[3][p]);
    pa.M(lujets_.p[4][p]);

    if (lujets_.k[2][p]!=0) {
#ifdef DEBUG
      std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] " << "(pdgId=" << pa.pdgId << ") has mother (pdgId=" << lujets_.k[1][lujets_.k[2][p]-1] << ")" << std::endl;
#endif
      pa.SetMother(ev_->GetById(lujets_.k[2][p]-1));
    }
    
    this->_hadrons->push_back(pa);
    //ev_->AddParticle(&pa);
  }
  return true;
}
