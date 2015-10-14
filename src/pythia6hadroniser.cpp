#include "pythia6hadroniser.h"

Pythia6Hadroniser::Pythia6Hadroniser()
{
  _name = "Pythia6";
  //this->pygive("MSTU(21)=1");
}

Pythia6Hadroniser::~Pythia6Hadroniser()
{
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Destructor called" << std::endl;
#endif
}

bool
Pythia6Hadroniser::Hadronise(Particle *part_)
{
  pyjets_.p[0][0] = part_->Px();
  pyjets_.p[1][0] = part_->Py();
  pyjets_.p[2][0] = part_->Pz();
  pyjets_.p[3][0] = part_->E();
  pyjets_.p[4][0] = part_->M();

  pyjets_.k[0][0] = 1; // status
  pyjets_.k[1][0] = 2; // particle id
  pyjets_.k[2][0] = 0; // mother
  pyjets_.k[3][0] = 0; // daughter 1
  pyjets_.k[4][0] = 0; // daughter 2

  this->pyexec();
  std::cout << __PRETTY_FUNCTION__ << " INFO" << std::endl;
  return true;
}

bool
Pythia6Hadroniser::Hadronise(Event *ev_)
{
  int np, status;
  std::vector<int> rl;
  std::vector<int>::iterator r;

  ParticlesRef pr, daug;
  ParticlesRef::iterator p;

  const int max_part_in_str = 3;
  const int max_str_in_evt = 2;

  int str_in_evt, part_in_str, num_part_in_str[max_str_in_evt];
  int jlrole[max_str_in_evt], jlpsf[max_str_in_evt][max_part_in_str];
  int oldnpart, criteria; //FIXME find an other name...
  
  this->PrepareHadronisation(ev_);

  rl = ev_->GetRoles();

  // First we initialise the string fragmentation variables
  for (int i=0; i<max_str_in_evt; i++) {
    jlrole[i] = -1;
    num_part_in_str[i] = 0;
    for (int j=0; j<max_part_in_str; j++) jlpsf[i][j] = -1;
  }

#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Dump of the event before the hadronisation" << std::endl;
  ev_->Dump();
#endif

  // Filling the common block to propagate to PYTHIA6
  pyjets_.n = 0;
  str_in_evt = 0;

  for (r=rl.begin(); r!=rl.end(); r++) {
    pr = ev_->GetByRole(*r);
    part_in_str = 0;
    for (p=pr.begin(); p!=pr.end(); p++) {
      np = (*p)->id;
      
      pyjets_.p[0][np] = (float)(*p)->Px();
      pyjets_.p[1][np] = (float)(*p)->Py();
      pyjets_.p[2][np] = (float)(*p)->Pz();
      pyjets_.p[3][np] = (float)(*p)->E();
      pyjets_.p[4][np] = (float)(*p)->M();

      if ((*p)->status<=0) status = 21;
      else status = (*p)->status;
      
      pyjets_.k[0][np] = status;
      pyjets_.k[1][np] = (int)((*p)->pdgId);
      
      //if ((*p)->GetMother()!=-1) pyjets_.k[2][np] = (*p)->GetMother()+1; // mother
      if ((*p)->GetMothersIds().size()>0) pyjets_.k[2][np] = *((*p)->GetMothersIds().begin())+1; // mother
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
        pyjets_.k[0][np] = 1; //FIXME PYTHIA/JETSET workaround
	jlrole[str_in_evt] = (*p)->role;
	jlpsf[str_in_evt][part_in_str] = (*p)->id+1;
	num_part_in_str[str_in_evt]++;
	part_in_str++;
      }
      pyjets_.n++;
    }
    if (jlrole[str_in_evt]>0) {
      str_in_evt++;
    }
  }

  oldnpart = pyjets_.n;

#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Passed the string construction stage." << std::endl
	    << "  " << str_in_evt << " string objects were identified and constructed" << std::endl;
#endif

  for (int i=0; i<str_in_evt; i++) {
    if (num_part_in_str[i]<2) continue;
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Joining " << num_part_in_str[i] << " particle in a same string (" << i << ") with role " << jlrole[i] << std::endl;
#endif
    for (int j=0; j<num_part_in_str[i]; j++) {
      if (jlpsf[i][j]==-1) continue;
#ifdef DEBUG
      std::cout << " * " << jlpsf[i][j] << " (pdgId=" << pyjets_.k[1][jlpsf[i][j]-1] << ")" << std::endl;
#endif
    }
    this->pyjoin(num_part_in_str[i], jlpsf[i]);
    //this->pyexec();//FIXME FIXME FIXME
  }
  this->pyexec();//FIXME FIXME FIXME
  //this->pylist(2);

  criteria = oldnpart+1;
  for (int i=0; i<max_str_in_evt; i++) criteria += num_part_in_str[i];

  if (pyjets_.k[1][criteria]==2212 and pyjets_.k[0][criteria]==1) {
    //this->pylist(2);
#ifdef ERROR
    std::cerr << __PRETTY_FUNCTION__ << " [ERROR] System non-inelastic" << std::endl;
#endif
    return false;
  }

  for (int p=0; p<pyjets_.n; p++) {

    //FIXME FIXME FIXME FIXME need to reimplement this first filter under this philosophy
    // First we filter the particles with status <= 0 :
    //  Status code = -1 : CLPAIR "internal" particles (not to be interacted with)
    //                 0 : Pythia6 empty lines
    //if (pyjets_.k[0][p]<=0) continue;
    //FIXME FIXME FIXME FIXME

    // We filter the first particles already present in the event
    if (p<oldnpart) continue;

    Particle pa;
    pa.id = p;
    pa.pdgId = static_cast<Particle::ParticleCode>(pyjets_.k[1][p]);
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
      std::cout << __PRETTY_FUNCTION__ << " [DEBUG] "
		<< pa.id << " (pdgId=" << pa.pdgId << ") has mother "
		<< pyjets_.k[2][p] << " (pdgId=" << pyjets_.k[1][pyjets_.k[2][p]-1] << ")"
		<< std::endl;
#endif
      pa.SetMother(ev_->GetById(pyjets_.k[2][p]-1));
    }

    ev_->AddParticle(pa);
  }

  return true;
}

bool
Pythia6Hadroniser::PrepareHadronisation(Event *ev_)
{
  Particle::ParticleCode singlet_id, doublet_id;
  double ranudq, ulmdq, ulmq;
  double ranmxp, ranmxt;
  double pmxp;
  double pmxda[4];
  double partpb[4];

#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Hadronisation preparation called !" << std::endl;
#endif

  ParticlesRef pp;
  ParticlesRef::iterator p;
  
  pp = ev_->GetParticles();
  for (p=pp.begin(); p!=pp.end(); p++) {
    if ((*p)->status==-2) { // One proton to be fragmented
      ranudq = (double)rand()/RAND_MAX;
      if (ranudq<1./9.) {
        singlet_id = Particle::QUARK_D;
        doublet_id = Particle::DIQUARK_UU1;
      }
      else if (ranudq<5./9.) {
        singlet_id = Particle::QUARK_U;
        doublet_id = Particle::DIQUARK_UD0;
      }
      else {
        singlet_id = Particle::QUARK_U;
        doublet_id = Particle::DIQUARK_UD1;
      }
      ulmdq = pymass(doublet_id);
      ulmq = pymass(singlet_id);

      // Choose random direction in MX frame
      ranmxp = 2.*pi*drand();       // phi angle
      ranmxt = acos(2.*drand()-1.); // theta angle

      // Compute momentum of decay particles from MX
      pmxp = std::sqrt(std::pow( (*p)->M2() - std::pow(ulmdq, 2) + std::pow(ulmq, 2), 2) / (4.*(*p)->M2()) - std::pow(ulmq, 2));

      // Build 4-vectors and boost decay particles
      pmxda[0] = pmxp*sin(ranmxt)*cos(ranmxp);
      pmxda[1] = pmxp*sin(ranmxt)*sin(ranmxp);
      pmxda[2] = pmxp*cos(ranmxt);
      pmxda[3] = sqrt(pow(pmxp, 2)+pow(ulmq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);

      if (!(partpb[0]<0) and !(partpb[0]>0)) {
        return false;
      }

      Particle singlet((*p)->role, singlet_id);
      singlet.status = 3;
      if (!singlet.P(partpb)) {
#ifdef ERROR
        std::cerr << __PRETTY_FUNCTION__ << " ERROR while setting the 4-momentum of singlet" << std::endl;
#endif
      }
      singlet.M(-1); //FIXME

      pmxda[0] = -pmxda[0];
      pmxda[1] = -pmxda[1];
      pmxda[2] = -pmxda[2];
      pmxda[3] = sqrt(pow(pmxp, 2)+pow(ulmq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);
      
      Particle doublet((*p)->role, doublet_id);
      doublet.status = 3;
      if (!doublet.P(partpb)) {
#ifdef ERROR
        std::cout << __PRETTY_FUNCTION__ << " ERROR while setting the 4-momentum of doublet" << std::endl;
#endif
      }
      //std::cout << "doublet, mass = " << doublet.M() << std::endl;
      doublet.M(-1); //FIXME

      if ((*p)->NumDaughters()==0) {
        singlet.SetMother(ev_->GetById((*p)->id));
        doublet.SetMother(ev_->GetById((*p)->id));

        ev_->AddParticle(singlet);
        ev_->AddParticle(doublet);
#ifdef DEBUG
        std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Quark/diquark content succesfully added to the event!" << std::endl;
#endif
      }
      else { // Quark/diquark content already present in the event
	std::vector<int> daugh;
	std::vector<int>::iterator did;

#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Quark/diquark content already present in the event!" << std::endl
		  << "  Role of these particles: " << (*p)->role << std::endl;
#endif
	daugh = (*p)->GetDaughters();
	for (did=daugh.begin(); did!=daugh.end(); did++) {
	  if (ev_->GetById(*did)->pdgId==1 or ev_->GetById(*did)->pdgId==2) { // Quark
	    singlet.SetMother(ev_->GetById((*p)->id));
	    *(ev_->GetById(*did)) = singlet;
#ifdef DEBUG
      std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Singlet replaced" << std::endl;
#endif
	  }
	  else { // Diquark
	    doublet.SetMother(ev_->GetById((*p)->id));
	    *(ev_->GetById(*did)) = doublet;
#ifdef DEBUG
      std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Doublet replaced" << std::endl;
#endif
	  }
	}
      }
    }
  }
  return true;
}

