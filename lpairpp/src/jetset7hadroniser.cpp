#include "jetset7hadroniser.h"

Jetset7Hadroniser::Jetset7Hadroniser()
{
  _name = "Jetset7";
}

Jetset7Hadroniser::~Jetset7Hadroniser()
{
#ifdef DEBUG
  std::cout << "[Jetset7Hadroniser::~Jetset7Hadroniser] [DEBUG] Destructor called" << std::endl;
#endif
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
  std::vector<int> rl;
  std::vector<int>::iterator r;

  std::vector<Particle*> pr;
  std::vector<Particle*>::iterator p;

  std::vector<Particle*> daug;

  std::ostringstream oss;
  const int max_part_in_str = 3;
  const int max_str_in_evt = 2;
  //bool isprimary;
  int id1, id2;
  int njoin[max_str_in_evt], jlrole[max_str_in_evt], jlpsf[max_str_in_evt][max_part_in_str];

  rl = ev_->GetRoles();

  // First we initialise the string fragmentation variables
  for (int i=0; i<max_str_in_evt; i++) {
    jlrole[i] = -1;
    njoin[i] = 0;
    for (int j=0; j<max_part_in_str; j++) jlpsf[i][j] = -1;
  }

#ifdef DEBUG
  std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] Dump of the event before the hadronisation" << std::endl;
  ev_->Dump();
#endif

  // Filling the common block to propagate to JETSET7
  lujets_.n = 0;

  for (r=rl.begin(), id1=0; r!=rl.end(); r++) {
    pr = ev_->GetByRole(*r);
    for (p=pr.begin(), id2=0; p!=pr.end(); p++) {
      np = (*p)->id;
      
      lujets_.p[0][np] = (float)(*p)->px;
      lujets_.p[1][np] = (float)(*p)->py;
      lujets_.p[2][np] = (float)(*p)->pz;
      lujets_.p[3][np] = (float)(*p)->E();
      lujets_.p[4][np] = (float)(*p)->M();
      
      lujets_.k[0][np] = (*p)->status-1;
      lujets_.k[1][np] = (*p)->pdgId;
      
      if ((*p)->GetMother()!=-1) {
	lujets_.k[2][np] = (*p)->GetMother()+1; // mother
      }
      else {
	lujets_.k[2][np] = 0; // mother
      }

      daug = ev_->GetDaughters(*p);
      if (daug.size()!=0) {
	lujets_.k[3][np] = (*p)->GetDaughters().front()+1; // daughter 1
	lujets_.k[4][np] = (*p)->GetDaughters().back()+1; // daughter 2
      }
      else {
	lujets_.k[3][np] = 0; // daughter 1
	lujets_.k[4][np] = 0; // daughter 2
      }
      
      for (int i=0; i<5; i++) {
	lujets_.v[i][np] = 0.;
      }
            
      if ((*p)->status==3) {
	jlrole[id1] = (*p)->role;
	jlpsf[id1][id2] = (*p)->id+1;
	njoin[id1]++;
	id2++;
      }
      lujets_.n++;
    }
    if (jlrole[id1]!=-1) {
      id1++;
    }
  }

#ifdef DEBUG
  std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] Passed the string construction stage" << std::endl;
#endif

  for (int i=0; i<max_str_in_evt; i++) {
    if (njoin[i]<2) continue;
#ifdef DEBUG
    std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] Joining " << njoin[i] << " particle in a same string (" << i << ") with role " << jlrole[i] << std::endl;
#endif
    for (int j=0; j<max_part_in_str; j++) {
      if (jlpsf[i][j]==-1) continue;
#ifdef DEBUG
      std::cout << " * " << jlpsf[i][j] << " (pdgId=" << lujets_.k[1][jlpsf[i][j]-1] << ")" << std::endl;
#endif
    }
    this->lujoin(njoin[i], jlpsf[i]);
  }
  this->luexec();
  //this->lulist(2);

  for (int p=0; p<lujets_.n; p++) {

    if (lujets_.k[0][p]==0) continue;

    Particle pa;
    pa.id = p;
    pa.pdgId = lujets_.k[1][p];
    if (ev_->GetById(lujets_.k[2][p]-1)!=(Particle*)NULL) {
      pa.role = ev_->GetById(lujets_.k[2][p]-1)->role; // Child particle inherits its mother's role
    }
    pa.status = lujets_.k[0][p];
    pa.P(lujets_.p[0][p], lujets_.p[1][p], lujets_.p[2][p], lujets_.p[3][p]);
    pa.M(lujets_.p[4][p]);
    pa.name = this->luname(pa.pdgId);
    pa.charge = this->luchge(pa.pdgId);

    if (lujets_.k[2][p]!=0) {
#ifdef DEBUG
      std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] "
		<< pa.id << " (pdgId=" << pa.pdgId << ") has mother "
		<< lujets_.k[2][p] << " (pdgId=" << lujets_.k[1][lujets_.k[2][p]-1] << ")"
		<< std::endl;
#endif
      pa.SetMother(ev_->GetById(lujets_.k[2][p]-1));
    }

    ev_->AddParticle(&pa);
  }

  return true;
}
