#include "pythia6hadroniser.h"

Pythia6Hadroniser::Pythia6Hadroniser()
{
  _name = "Pythia6";
  this->pygive("MSTU(21)=1");
}

Pythia6Hadroniser::~Pythia6Hadroniser()
{
#ifdef DEBUG
  std::cout << "[Pythia6Hadroniser::~Pythia6Hadroniser] [DEBUG] Destructor called" << std::endl;
#endif
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
  std::cout << "[Pythia6Hadroniser::Hadronise] INFO" << std::endl;
  return true;
}

bool
Pythia6Hadroniser::Hadronise(Event *ev_)
{
  int np;
  bool quarks_built;
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
  
  quarks_built = this->PrepareHadronisation(ev_);
  if (!quarks_built) return quarks_built;

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
      
      if ((*p)->GetMother()!=-1) pyjets_.k[2][np] = (*p)->GetMother()+1; // mother
      else pyjets_.k[2][np] = 0; // mother
      
      daug = ev_->GetDaughters(*p);
      if (daug.size()!=0) {
	pyjets_.k[3][np] = (*p)->GetDaughters().front()+1; // daughter 1
	pyjets_.k[4][np] = (*p)->GetDaughters().back()+1; // daughter 2
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
  this->pyexec();
  //this->pylist(2);

  for (int p=0; p<pyjets_.n; p++) {

    // First we filter the particles with status <= 0 :
    //  Status code = -1 : CLPAIR "internal" particles (not to be interacted with)
    //                 0 : Pythia6 empty lines
    //if (pyjets_.k[0][p]<=0 or pyjets_.k[0][p]==13 or pyjets_.k[0][p]==12) continue;
    if (pyjets_.k[0][p]<=0 or pyjets_.k[0][p]==13) continue;

    Particle pa;
    pa.id = p;
    pa.pdgId = pyjets_.k[1][p];
    if (ev_->GetById(pyjets_.k[2][p]-1)!=(Particle*)NULL) {
      pa.role = ev_->GetById(pyjets_.k[2][p]-1)->role; // Child particle inherits its mother's role
    }
    pa.status = pyjets_.k[0][p];
    pa.P(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]);
    pa.M(pyjets_.p[4][p]);
    pa.name = this->pyname(pa.pdgId);
    pa.charge = (float)(this->pyp(p+1,6));

    if (pyjets_.k[2][p]!=0) {
#ifdef DEBUG
      std::cout << "[Pythia6Hadroniser::Hadronise] [DEBUG] "
		<< pa.id << " (pdgId=" << pa.pdgId << ") has mother "
		<< pyjets_.k[2][p] << " (pdgId=" << pyjets_.k[1][pyjets_.k[2][p]-1] << ")"
		<< std::endl;
#endif
      pa.SetMother(ev_->GetById(pyjets_.k[2][p]-1));
    }

    ev_->AddParticle(&pa);
  }

  return true;
}

bool
Pythia6Hadroniser::PrepareHadronisation(Event *ev_)
{
  int singlet_id, doublet_id;
  double ranudq, ulmdq, ulmq;
  double ranmxp, ranmxt;
  double pmxp;
  double pmxda[4];
  double partpb[4];

#ifdef DEBUG
  std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Hadronisation preparation called !" << std::endl;
#endif

  Particles pp;
  Particles::iterator p;
  
  pp = ev_->GetParticles();
  for (p=pp.begin(); p!=pp.end(); p++) {
    if ((*p)->status==-2) { // One proton to be fragmented
      ranudq = (double)rand()/RAND_MAX;
      if (ranudq<1./9.) {
        singlet_id = 1;
        doublet_id = 2203;
      }
      else if (ranudq<5./9.) {
        singlet_id = 2;
        doublet_id = 2101;
      }
      else {
        singlet_id = 2;
        doublet_id = 2103;
      }
      ulmdq = pymass(doublet_id);
      ulmq = pymass(singlet_id);

      // Choose random direction in MX frame
      ranmxp = 2.*pi*(double)rand()/RAND_MAX;
      ranmxt = acos(2.*(double)rand()/RAND_MAX-1.);

      // Compute momentum of decay particles from MX
      pmxp = std::sqrt(std::pow( (*p)->M2() - std::pow(ulmdq, 2) + std::pow(ulmq, 2), 2) / (4.*(*p)->M2()) - std::pow(ulmq, 2));

      // Build 4-vectors and boost decay particles
      pmxda[0] = sin(ranmxt)*cos(ranmxp)*pmxp;
      pmxda[1] = sin(ranmxt)*sin(ranmxp)*pmxp;
      pmxda[2] = cos(ranmxt)*pmxp;
      pmxda[3] = sqrt(pow(pmxp, 2)+pow(ulmdq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);

      if (!(partpb[0]<0) and !(partpb[0]>0)) {
        return false;
      }

      Particle singlet((*p)->role, singlet_id);
      singlet.status = 3;
      singlet.SetMother(ev_->GetOneByRole((*p)->role));
      if (!singlet.P(partpb)) {
    #ifdef ERROR
        std::cerr << "[GamGam::PrepareHadronisation] ERROR while setting the 4-momentum of singlet" << std::endl;
    #endif
      }
      ev_->AddParticle(&singlet);

      pmxda[0] = -pmxda[0];
      pmxda[1] = -pmxda[1];
      pmxda[2] = -pmxda[2];
      pmxda[3] = sqrt(pow(pmxp, 2)+pow(ulmq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);
      
      Particle doublet((*p)->role, doublet_id);
      doublet.status = 3;
      doublet.SetMother(ev_->GetOneByRole((*p)->role));
      if (!doublet.P(partpb)) {
    #ifdef ERROR
        std::cout << "[GamGam::PrepareHadronisation] ERROR while setting the 4-momentum of doublet" << std::endl;
    #endif
      }
      ev_->AddParticle(&doublet);
    }
  }
  return true;
}

