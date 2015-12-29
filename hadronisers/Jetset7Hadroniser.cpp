#include "Jetset7Hadroniser.h"

Jetset7Hadroniser::Jetset7Hadroniser() : GenericHadroniser("Jetset7")
{
  //this->lugive("MSTU(21)=1");
  //this->lugive("MSTJ(1)=1");
}

Jetset7Hadroniser::~Jetset7Hadroniser()
{
  //Debug("Destructor called");
}

bool
Jetset7Hadroniser::Hadronise(Particle *part_)
{
  lujets_.p[0][0] = part_->GetMomentum().Px();
  lujets_.p[1][0] = part_->GetMomentum().Py();
  lujets_.p[2][0] = part_->GetMomentum().Pz();
  lujets_.p[3][0] = part_->E();
  lujets_.p[4][0] = part_->M();

  lujets_.k[0][0] = 1; // status
  lujets_.k[1][0] = 2; // particle id
  lujets_.k[2][0] = 0; // mother 1
  lujets_.k[3][0] = 0; // mother 2
  lujets_.k[4][0] = 0; // daughter

  this->luexec();
  return true;
}

bool
Jetset7Hadroniser::Hadronise(Event *ev_)
{
  int np, oldnpart;
  bool quarks_built;
  ParticleRoles rl;
  ParticleRoles::iterator r;

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
  
  if (Logger::GetInstance()->Level>=Logger::Debug) {
    Debug("Dump of the event before the hadronisation");
    ev_->Dump();
  }
  
  // Filling the common block to propagate to JETSET7
  lujets_.n = 0;

  int status;

  for (r=rl.begin(), id1=0; r!=rl.end(); r++) {
    pr = ev_->GetByRole(*r);
    //std::cout << "--> role " << *r << " contains " << pr.size() << " particles" << std::endl;
    for (p=pr.begin(), id2=0; p!=pr.end(); p++) {
      np = (*p)->id;
      
      //(*p)->Dump();

      lujets_.p[0][np] = (float)(*p)->GetMomentum().Px();
      lujets_.p[1][np] = (float)(*p)->GetMomentum().Py();
      lujets_.p[2][np] = (float)(*p)->GetMomentum().Pz();
      lujets_.p[3][np] = (float)(*p)->E();
      lujets_.p[4][np] = (float)(*p)->M();

      if ((*p)->status==-1 or (*p)->status==0) status = 21;
      else status = (*p)->status;
      
      lujets_.k[0][np] = status;
      lujets_.k[1][np] = static_cast<Particle::ParticleCode>((*p)->GetPDGId());
      
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

  std::ostringstream dbg;

  oldnpart = lujets_.n;
  //this->lulist(2);

#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Passed the string construction stage" << std::endl;
#endif

  for (int i=0; i<max_str_in_evt; i++) {
    if (njoin[i]<2) continue;
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Joining " << njoin[i] << " particle in a same string (" << i << ") with role " << jlrole[i] << std::endl;
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
    pa.SetPDGId(static_cast<Particle::ParticleCode>(lujets_.k[1][p]));
    if (ev_->GetById(lujets_.k[2][p]-1)!=(Particle*)NULL) {
      pa.role = ev_->GetById(lujets_.k[2][p]-1)->role; // Child particle inherits its mother's role
    }
    pa.status = static_cast<Particle::Status>(lujets_.k[0][p]);
    pa.SetMomentum(Particle::Momentum(lujets_.p[0][p], lujets_.p[1][p], lujets_.p[2][p], lujets_.p[3][p]));
    pa.SetM(lujets_.p[4][p]);
    pa.name = this->luname(pa.GetPDGId());
    pa.charge = this->luchge(pa.GetPDGId());

    if (lujets_.k[2][p]!=0) {
      dbg << Form("\n\t%2d (pdgId=%4d) has mother %2d (pdgId=%4d)", pa.id, pa.GetPDGId(), lujets_.k[2][p], lujets_.k[1][lujets_.k[2][p]-1]);
      pa.SetMother(ev_->GetById(lujets_.k[2][p]-1));
    }

    ev_->AddParticle(pa);
  }
  Debug(Form("Passed the string construction stage.\n\t %d string objects were identified and constructed",
             "%s", max_str_in_evt, dbg.str().c_str()));

  return true;
}

bool
Jetset7Hadroniser::PrepareHadronisation(Event *ev_)
{
  Particle::ParticleCode singlet_id, doublet_id;
  double ranudq, ulmdq, ulmq;
  double ranmxp, ranmxt;
  double pmxp;
  double pmxda[4];
  double partpb[4];

  Debug("Hadronisation preparation called!");

  ParticlesRef pp;
  ParticlesRef::iterator p;
  
  pp = ev_->GetParticles();
  for (p=pp.begin(); p!=pp.end(); p++) {
    if ((*p)->status==Particle::Undecayed) continue;
    // One proton to be fragmented
    ranudq = drand();
    if (ranudq<1./9.) {
      singlet_id = Particle::dQuark;
      doublet_id = Particle::uu1Diquark;
    }
    else if (ranudq<5./9.) {
      singlet_id = Particle::uQuark;
      doublet_id = Particle::ud0Diquark;
    }
    else {
      singlet_id = Particle::uQuark;
      doublet_id = Particle::ud1Diquark;
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
    
    Lorenb((*p)->M(), (*p)->GetMomentum(), pmxda, partpb);
    
    if (!(partpb[0]<0) and !(partpb[0]>0)) return false;
    
    Particle singlet((*p)->role, singlet_id);
    singlet.status = Particle::DebugResonance;
    //singlet.SetMother(ev_->GetOneByRole((*p)->role));
    if (!singlet.SetMomentum(partpb)) {
      throw Exception(__PRETTY_FUNCTION__, "ERROR while setting the 4-momentum of singlet", JustWarning);
    }
    //std::cout << "singlet, mass = " << singlet.M() << std::endl;
    //singlet.Dump();
    singlet.SetM(); //FIXME
    //ev_->AddParticle(singlet);
    
    pmxda[0] = -pmxda[0];
    pmxda[1] = -pmxda[1];
    pmxda[2] = -pmxda[2];
    pmxda[3] = std::sqrt(std::pow(pmxp, 2)+std::pow(ulmdq, 2));
    
    Lorenb((*p)->M(), (*p)->GetMomentum(), pmxda, partpb);
    
    Particle doublet((*p)->role, doublet_id);
    doublet.status = Particle::DebugResonance;
    doublet.SetMother(ev_->GetOneByRole((*p)->role));
    if (!doublet.SetMomentum(partpb)) {
      throw Exception(__PRETTY_FUNCTION__, "ERROR while setting the 4-momentum of doublet", JustWarning);
    }
    //std::cout << "doublet, mass = " << doublet.M() << std::endl;
    doublet.SetM(); //FIXME
    
    if ((*p)->NumDaughters()==0) {
      singlet.SetMother(ev_->GetById((*p)->id));
      doublet.SetMother(ev_->GetById((*p)->id));

      ev_->AddParticle(singlet);
      ev_->AddParticle(doublet);
        
      Debug("Quark/diquark content succesfully added to the event!");
    }
    else { // Quark/diquark content already present in the event
      std::vector<int> daugh;
      std::vector<int>::iterator did;
      
      Debug(Form("Quark/diquark content already present in the event!\n\tRole of these particles: %d", (*p)->role));
      
      daugh = (*p)->GetDaughters();
      for (did=daugh.begin(); did!=daugh.end(); did++) {
        if (ev_->GetById(*did)->GetPDGId()==Particle::uQuark
         or ev_->GetById(*did)->GetPDGId()==Particle::dQuark) { // Quark
          singlet.SetMother(ev_->GetById((*p)->id));
          *(ev_->GetById(*did)) = singlet;
          Debug("Singlet replaced");
        }
        else { // Diquark
          doublet.SetMother(ev_->GetById((*p)->id));
          *(ev_->GetById(*did)) = doublet;
          Debug("Doublet replaced");
        }
      }
    }
  }
  return true;
}

