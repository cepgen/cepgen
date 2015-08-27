#include "jetset7hadroniser.h"

Jetset7Hadroniser::Jetset7Hadroniser()
{
  _name = "Jetset7";
  //this->lugive("MSTU(21)=1");
  //this->lugive("MSTJ(1)=1");
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
  lujets_.p[0][0] = part_->Px();
  lujets_.p[1][0] = part_->Py();
  lujets_.p[2][0] = part_->Pz();
  lujets_.p[3][0] = part_->E();
  lujets_.p[4][0] = part_->M();

  lujets_.k[0][0] = 1; // status
  lujets_.k[1][0] = 2; // particle id
  lujets_.k[2][0] = 0; // mother 1
  lujets_.k[3][0] = 0; // mother 2
  lujets_.k[4][0] = 0; // daughter

  this->luexec();
  std::cout << "[Jetset7Hadroniser::Hadronise] INFO" << std::endl;
  return true;
}

bool
Jetset7Hadroniser::Hadronise(Event *ev_)
{
  int np, oldnpart;
  bool quarks_built;
  std::vector<int> rl;
  std::vector<int>::iterator r;

  ParticlesRef pr, daug;
  ParticlesRef::iterator p;

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
  std::cout << "[Jetset7Hadroniser::Hadronise] [DEBUG] Dump of the event before the hadronisation" << std::endl;
  ev_->Dump();
#endif

  // Filling the common block to propagate to JETSET7
  lujets_.n = 0;

  int status;

  for (r=rl.begin(), id1=0; r!=rl.end(); r++) {
    pr = ev_->GetByRole(*r);
    //std::cout << "--> role " << *r << " contains " << pr.size() << " particles" << std::endl;
    for (p=pr.begin(), id2=0; p!=pr.end(); p++) {
      np = (*p)->id;
      
      //(*p)->Dump();

      lujets_.p[0][np] = (float)(*p)->Px();
      lujets_.p[1][np] = (float)(*p)->Py();
      lujets_.p[2][np] = (float)(*p)->Pz();
      lujets_.p[3][np] = (float)(*p)->E();
      lujets_.p[4][np] = (float)(*p)->M();

      if ((*p)->status==-1 or (*p)->status==0) status = 21;
      else status = (*p)->status;
      
      lujets_.k[0][np] = status;
      lujets_.k[1][np] = (ParticleId)((*p)->pdgId);
      
      if ((*p)->GetMothersIds().size()>0) lujets_.k[2][np] = *((*p)->GetMothersIds().begin())+1; // mother
      else lujets_.k[2][np] = 0; // no mother registered
      
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
        //FIXME workaround
        lujets_.k[0][np] = 1;
        //FIXME
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

  oldnpart = lujets_.n;
  //this->lulist(2);

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

    //FIXME FIXME FIXME FIXME need to reimplement this first filter under this philosophy
    // First we filter the particles with status <= 0 :
    //  Status code = -1 : CLPAIR "internal" particles (not to be interacted with)
    //                 0 : Jetset7 empty lines
    //if (lujets_.k[0][p]<=0) continue;
    //FIXME FIXME FIXME FIXME

    // We filter the first particles already present in the event
    if (p<oldnpart) continue;

    Particle pa;
    pa.id = p;
    pa.pdgId = (ParticleId)(lujets_.k[1][p]);
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

    ev_->AddParticle(pa);
  }

  ev_->Dump();

  return true;
}

bool
Jetset7Hadroniser::PrepareHadronisation(Event *ev_)
{
  ParticleId singlet_id, doublet_id;
  double ranudq, ulmdq, ulmq;
  double ranmxp, ranmxt;
  double pmxp;
  double pmxda[4];
  double partpb[4];

#ifdef DEBUG
  std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Hadronisation preparation called !" << std::endl;
#endif

  ParticlesRef pp;
  ParticlesRef::iterator p;
  
  //ev_->Dump();

  pp = ev_->GetParticles();
  for (p=pp.begin(); p!=pp.end(); p++) {
    if ((*p)->status==-2) { // One proton to be fragmented

      //(*p)->Dump();
      //std::cout << "m2 = " << (*p)->M2() << std::endl;
      
      ranudq = drand();
      if (ranudq<1./9.) {
        singlet_id = QUARK_D;
        doublet_id = DIQUARK_UU1;
      }
      else if (ranudq<5./9.) {
        singlet_id = QUARK_U;
        doublet_id = DIQUARK_UD0;
      }
      else {
        singlet_id = QUARK_U;
        doublet_id = DIQUARK_UD1;
      }
      ulmdq = ulmass(doublet_id);
      ulmq = ulmass(singlet_id);

      // Choose random direction in MX frame
      ranmxp = 2.*pi*drand();       // phi angle
      ranmxt = acos(2.*drand()-1.); // theta angle

      // Compute momentum of decay particles from MX
      pmxp = std::sqrt(std::pow( (*p)->M2() - std::pow(ulmdq, 2) + std::pow(ulmq, 2), 2) / (4.*(*p)->M2()) - std::pow(ulmq, 2));

      /*if (!(pmxda[0]<0) and !(pmxda[0]>0)) { //NaN
	std::cout << "-----> " << pmxp << ", " << ranmxt << ", " << ranmxp << std::endl;
	}*/

      // Build 4-vectors and boost decay particles
      pmxda[0] = pmxp*sin(ranmxt)*cos(ranmxp);
      pmxda[1] = pmxp*sin(ranmxt)*sin(ranmxp);
      pmxda[2] = pmxp*cos(ranmxt);
      pmxda[3] = std::sqrt(std::pow(pmxp, 2)+std::pow(ulmq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);

      if (!(partpb[0]<0) and !(partpb[0]>0)) {
	/*std::cout << "=== " << pmxp << "\t" << ulmdq << "\t" << ulmq << "\t" << (*p)->M() << std::endl;
	for (int i=0; i<4; i++) {
	  std::cout << "(" << i << ")-> " << pmxda[i] << " --> " << partpb[i] << std::endl;
	  }*/
	return false;
      }

      Particle singlet((*p)->role, singlet_id);
      singlet.status = 3;
      //singlet.SetMother(ev_->GetOneByRole((*p)->role));
      if (!singlet.P(partpb)) {
	//#ifdef ERROR
        std::cerr << "[GamGam::PrepareHadronisation] ERROR while setting the 4-momentum of singlet" << std::endl;
	//#endif
      }
      //std::cout << "singlet, mass = " << singlet.M() << std::endl;
      //singlet.Dump();
      singlet.M(-1); //FIXME
      //ev_->AddParticle(singlet);

      pmxda[0] = -pmxda[0];
      pmxda[1] = -pmxda[1];
      pmxda[2] = -pmxda[2];
      pmxda[3] = std::sqrt(std::pow(pmxp, 2)+std::pow(ulmdq, 2));

      Lorenb((*p)->M(), (*p)->P4(), pmxda, partpb);
      
      Particle doublet((*p)->role, doublet_id);
      doublet.status = 3;
      doublet.SetMother(ev_->GetOneByRole((*p)->role));
      if (!doublet.P(partpb)) {
	//#ifdef ERROR
        std::cout << "[GamGam::PrepareHadronisation] ERROR while setting the 4-momentum of doublet" << std::endl;
	//#endif
      }
      //std::cout << "doublet, mass = " << doublet.M() << std::endl;
      doublet.M(-1); //FIXME

      if ((*p)->NumDaughters()==0) {
        singlet.SetMother(ev_->GetById((*p)->id));
        doublet.SetMother(ev_->GetById((*p)->id));

        ev_->AddParticle(singlet);
        ev_->AddParticle(doublet);
#ifdef DEBUG
        std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Quark/diquark content succesfully added to the event!" << std::endl;
#endif
      }
      else { // Quark/diquark content already present in the event
	      std::vector<int> daugh;
	      std::vector<int>::iterator did;

#ifdef DEBUG
	      std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Quark/diquark content already present in the event!" << std::endl
		        << "  Role of these particles: " << (*p)->role << std::endl;
#endif
	      daugh = (*p)->GetDaughters();
	      for (did=daugh.begin(); did!=daugh.end(); did++) {
	        if (ev_->GetById(*did)->pdgId==1 or ev_->GetById(*did)->pdgId==2) { // Quark
	          singlet.SetMother(ev_->GetById((*p)->id));
	          *(ev_->GetById(*did)) = singlet;
#ifdef DEBUG
            std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Singlet replaced" << std::endl;
#endif
	        }
	        else { // Diquark
	          doublet.SetMother(ev_->GetById((*p)->id));
	          *(ev_->GetById(*did)) = doublet;
#ifdef DEBUG
            std::cout << "[GamGam::PrepareHadronisation] [DEBUG] Doublet replaced" << std::endl;
#endif
	        }
	      }
      }
    }
  }
  return true;
}

