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
  std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] Dump of the event before the hadronisation" << std::endl;
  ev_->Dump();
#endif

  // Muting the unneeded information
  //this->pygive("");

  // Filling the common block to propagate to PYTHIA6
  pyjets_.n = 0;

  for (r=rl.begin(), id1=0; r!=rl.end(); r++) {
    pr = ev_->GetByRole(*r);
    for (p=pr.begin(), id2=0; p!=pr.end(); p++) {
      np = (*p)->id;
      
      pyjets_.p[0][np] = (*p)->px;
      pyjets_.p[1][np] = (*p)->py;
      pyjets_.p[2][np] = (*p)->pz;
      pyjets_.p[3][np] = (*p)->E();
      pyjets_.p[4][np] = (*p)->M();
      
      pyjets_.k[0][np] = (*p)->status;
      pyjets_.k[1][np] = (*p)->pdgId;
      
      if ((*p)->GetMother()!=(Particle*)NULL) {
	pyjets_.k[2][np] = (*p)->GetMother()->id+1; // mother
      }
      else {
	pyjets_.k[2][np] = 0; // mother
      }
      
      daug = (*p)->GetDaughters();
      if (daug.size()!=0) {
	pyjets_.k[3][np] = (*p)->GetDaughters().front()->id+1; // daughter 1
	pyjets_.k[4][np] = (*p)->GetDaughters().back()->id+1; // daughter 2
      }
      else {
	pyjets_.k[3][np] = 0; // daughter 1
	pyjets_.k[4][np] = 0; // daughter 2
      }
      
      for (int i=0; i<5; i++) {
	pyjets_.v[i][np] = 0.;
      }
      
      if ((*p)->status==3) {
	jlrole[id1] = (*p)->role;
	jlpsf[id1][id2] = (*p)->id+1;
	njoin[id1]++;
	id2++;
      }
      pyjets_.n++;
    }
    if (jlrole[id1]!=-1) {
      id1++;
    }
  }

#ifdef DEBUG
  std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] Passed the string construction stage" << std::endl;
#endif

  for (int i=0; i<max_str_in_evt; i++) {
    if (njoin[i]<2) continue;
#ifdef DEBUG
    std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] Joining " << njoin[i] << " particle in a same string (" << i << ") with role " << jlrole[i] << std::endl;
#endif
    for (int j=0; j<max_part_in_str; j++) {
      if (jlpsf[i][j]==-1) continue;
#ifdef DEBUG
      std::cout << " * " << jlpsf[i][j] << " (pdgId=" << pyjets_.k[1][jlpsf[i][j]-1] << ")" << std::endl;
#endif
    }
    this->pyjoin(njoin[i], jlpsf[i]);
  }

  //this->pylist(2);

  this->pyexec();
  //this->pylist(2);

  for (int p=0; p<pyjets_.n; p++) {

    if (pyjets_.k[0][p]==0) continue;

    Particle pa;
    pa.id = p;
    if (ev_->GetById(pyjets_.k[2][p]-1)!=(Particle*)NULL) {
      pa.role = ev_->GetById(pyjets_.k[2][p]-1)->role; // Child particle inherits its mother's role
    }
    pa.status = pyjets_.k[0][p];
    pa.pdgId = pyjets_.k[1][p];
    pa.P(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]);
    pa.M(pyjets_.p[4][p]);
    pa.name = this->pyname(pa.pdgId);
    pa.charge = (float)(this->pyp(p+1,6));

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
