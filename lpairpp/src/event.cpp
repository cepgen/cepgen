#include "event.h"

Event::Event() :
  num_hadronisation_trials(0),
  time_generation(-1.), time_total(-1.)
{
  //this->_part = new ParticlesMap();
}

Event::~Event()
{
  //delete this->_part;
}

Event&
Event::operator=(const Event &ev_)
{
  this->_part = ev_._part;
  this->time_generation = ev_.time_generation;
  this->time_total = ev_.time_total;
  this->num_hadronisation_trials = ev_.num_hadronisation_trials;
  return *this;
}

ParticlesRef
Event::GetByRole(int role_)
{
  int i;
  ParticlesRef out;
  std::pair<ParticlesMap::iterator,ParticlesMap::iterator> ret = this->_part.equal_range(role_);
  ParticlesMap::iterator it;

  for (it=ret.first, i=0; it!=ret.second && i<100; it++, i++) {
    out.push_back(&(it->second));
  }
  return out;
}

Particle*
Event::GetById(int id_)
{
  ParticlesMap::iterator out;
  for (out=this->_part.begin(); out!=this->_part.end(); out++) {
    if (out->second.id==id_) {
      return &out->second;
    }
  }
  return 0;
}

const Particle
Event::GetConstById(int id_) const
{
  ParticlesMap::const_iterator out;
  for (out=this->_part.begin(); out!=this->_part.end(); out++) {
    if (out->second.id==id_) {
      return static_cast<Particle>(out->second);
    }
  }
  return 0;
}

std::vector<int>
Event::GetRoles() const
{
  ParticlesMap::const_iterator it, end;
  std::vector<int> out;
  for (it=this->_part.begin(), end=this->_part.end(); it!=end; it=this->_part.upper_bound(it->first)) {
    out.push_back(it->first);
  }
  return out;
}

int
Event::AddParticle(Particle part_, bool replace_)
{
#ifdef DEBUG
  std::cout << "[Event::AddParticle] [DEBUG] Particle with PDGid = " << part_.pdgId << " has role " << part_.role << std::endl;
#endif
  if (part_.role<=0) {
    return -1;
  }
  ParticlesRef part_with_same_role = this->GetByRole(part_.role);
  part_.id = this->_part.size(); //FIXME is there any better way of introducing this id ?
  if (replace_ and part_with_same_role.size()!=0) {
    part_with_same_role.at(0) = &part_;
    return 0;
  }
  this->_part.insert(std::pair<int,Particle>(part_.role, part_));
  return 1;
}

int
Event::AddParticle(int role_, bool replace_)
{
  int out;
  if (role_<=0) {
    return -1;
  }
  np = new Particle();
  np->role = role_;
  out = this->AddParticle(*np, replace_);

  delete np;
  return out;
}

/*HEPEUP
Event::GetHEPEUP() const
{
  const Particles particles = this->GetConstParticles();
  Particles::const_iterator p;
  ParticlesRef daughters;
  ParticlesRef::iterator dg;
  int np, status;

  std::cout << "particles = " << this->NumParticles() << std::endl;

  HEPEUP out;
  out.nup = this->NumParticles();
  out.idprup = 1; //FIXME what about multiple subprocesses ?

  std::cout << "aaa -> " << out.GetLHERecord() << std::endl;

  np = 0;
  for (p=particles.begin(); p!=particles.end(); p++) {
    std::cout << "--> " << p->pdgId << std::endl;
    out.idup[np] = p->pdgId;
    status = (p->status==0) ? 1 : p->status;
    out.istup[np] = status;
    if (!p->Primary()) {
      out.mothup[np][0] = *(p->GetMothersIds().begin())+1;
      out.mothup[np][1] = *(p->GetMothersIds().end())+1;
    }
    for (int i=0; i<5; i++) {
      out.pup[np][i] = p->P(i);
    }
    out.spinup[np] = p->helicity;
    np++;
    std::cout << out.istup[np] << ", " << out.pup[np][0] << std::endl;
  }
  return out;
  }*/

/*std::string
Event::GetLHERecord(const double weight_)
{
  std::stringstream ss;
  ParticlesRef particles, daughters;
  ParticlesRef::iterator p, dg;
  int min_id, max_id;

  //FIXME need to fetch the vector (not the multimap), so that we can sort on the particle unique identifier (also TODO!!!)
  ss << "<event>" << std::endl;
  ss << this->NumParticles() << "\t"
     << this->event_info.idprup << "\t"
     << this->event_info.xwgtup << "\t"
     << this->event_info.scalup << "\t" // scale of the event, in GeV
     << this->event_info.aqedup << "\t" // alphaQED
     << this->event_info.aqcdup // alphaQCD
     << std::endl;
  particles = this->GetParticles();
  for (p=particles.begin(); p!=particles.end(); p++) {
    if ((*p)->status==0) (*p)->status = 1;
    
    ss << std::setw(4) << (*p)->id+1 << "  "
       << std::setw(3) << (*p)->status << "  "
       << std::setw(5) << (*p)->pdgId << "  ";
    if (!(*p)->Primary()) {
      ss << std::setw(2) << *((*p)->GetMothersIds().begin())+1 << "  ";
    }
    else {
      ss << std::setw(2) << "0" << "  ";
    }
    daughters = this->GetDaughters(*p);
    max_id = 0;
    min_id = 999;
    if (daughters.size()>0) {
      for (dg=daughters.begin(); dg!=daughters.end(); dg++) {
	      if ((*dg)->id>this->NumParticles() or (*dg)->id<0) continue; //FIXMEEEEEEEEEEEEEEEEEEEEEE
	      if ((*dg)->id>max_id) max_id = (*dg)->id;
	      if ((*dg)->id<min_id) min_id = (*dg)->id;
      }
      if (min_id==max_id) 
      	ss << std::setw(4) << min_id+1 << "  " << std::setw(4) << "0";
      else 
      	ss << std::setw(4) << min_id+1 << "  " << std::setw(4) << max_id+1;
    }
    else {
      ss << std::setw(4) << "0" << "  " << std::setw(4) << "0";
    }
    ss << std::setw(4) << "  " << "0" << "  "; 
    ss << std::setw(12) << (*p)->Px() << "  " 
       << std::setw(12) << (*p)->Py() << "  " 
       << std::setw(12) << (*p)->Pz() << "  " 
       << std::setw(12) << (*p)->E() << "  " 
       << std::setw(12) << (*p)->M() << "  " 
       << std::setw(4) << weight_ 
       << std::endl; 
  }
  ss << "</event>" << std::endl;

  return ss.str();
}
*/
void
Event::Store(std::ofstream *of_, double weight_)
{
  Particle *l1 = this->GetByRole(6).at(0);
  Particle *l2 = this->GetByRole(7).at(0);
  
  *of_ << std::setw(8) << l1->E() << "\t"
       << std::setw(8) << l1->Px() << "\t"
       << std::setw(8) << l1->Py() << "\t"
       << std::setw(8) << l1->Pz() << "\t"
       << std::setw(8) << l1->Pt() << "\t"
       << std::setw(8) << l1->M() << "\t"
       << std::setw(8) << l1->Eta() << "\t"
       << std::setw(8) << l1->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
  *of_ << std::setw(8) << l2->E() << "\t"
       << std::setw(8) << l2->Px() << "\t"
       << std::setw(8) << l2->Py() << "\t"
       << std::setw(8) << l2->Pz() << "\t"
       << std::setw(8) << l2->Pt() << "\t"
       << std::setw(8) << l2->M() << "\t"
       << std::setw(8) << l2->Eta() << "\t"
       << std::setw(8) << l2->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
}

ParticlesRef
Event::GetParticles()
{
  ParticlesRef out;
  ParticlesMap::iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    out.push_back(&it->second);
  }
  std::sort(out.begin(), out.end(), compareParticlePtrs);
  return out;
}

Particles
Event::GetConstParticles() const
{
  Particles out;
  ParticlesMap::const_iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    out.push_back(static_cast<Particle>(it->second));
  }
  std::sort(out.begin(), out.end(), compareParticle);
  return out;
}

ParticlesRef
Event::GetStableParticles()
{
  ParticlesRef out;
  ParticlesMap::iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    if (it->second.status==0 or it->second.status==1) {
      out.push_back(&it->second);
    }
  }
  std::sort(out.begin(), out.end());
  return out;  
}

void
Event::Dump(bool stable_)
{
  ParticlesRef particles;
  ParticlesRef::iterator p;
  double pxtot, pytot, pztot, etot;
  int sign;

  pxtot = pytot = pztot = etot = 0.;
  particles = this->GetParticles();
  std::cout << "[Event::Dump]" << std::endl;
  std::cout << std::left;
  std::cout << "Particle" << "\t" << "PDG id" << "\t\t" << "Charge" << "\t" << "Role" << "\t" << "Status" << "\t" << "Mother" << "\t\t\t" << "4-Momentum [GeV]" << std::endl;
  std::cout << "--------" << "\t" << "------" << "\t\t" << "------" << "\t" << "----" << "\t" << "------" << "\t" << "------" << "\t" << "---------------------------------------" << std::endl;
  for (p=particles.begin(); p!=particles.end(); p++) {
    if (stable_ and (*p)->status!=1) continue;
    std::cout << std::setfill(' ') << std::setw(8) << (*p)->id
	      << "\t" << std::setw(6) << (*p)->pdgId;
    if ((*p)->name!="") {
      std::cout << std::setw(6) << (*p)->name;
    }
    else std::cout << "\t";
    std::cout << "\t";
    if ((*p)->charge!=999.)
    	std::cout << std::setprecision(2) << std::setw(6) << (*p)->charge;
    std::cout << "\t" << std::setw(4) << (*p)->role
	      << "\t" << std::setw(6) << (*p)->status << "\t";
    if ((*p)->GetMothersIds().size()>0) {
      std::cout << std::setw(2) << *((*p)->GetMothersIds().begin()) 
		<< " (" << std::right << std::setw(2) << this->GetById(*((*p)->GetMothersIds().begin()))->role << std::left << ") ";
    }
    else std::cout << std::setw(8) << "";
    std::cout << std::right;
    std::cout << std::setprecision(3) << std::setw(9) << (*p)->Px() << " ";
    std::cout << std::setprecision(3) << std::setw(9) << (*p)->Py() << " ";
    std::cout << std::setprecision(3) << std::setw(9) << (*p)->Pz() << " ";
    std::cout << std::setprecision(3) << std::setw(9) << (*p)->E() << " ";
    std::cout << std::left;
    std::cout << std::endl;
    if ((*p)->status>=0 and (*p)->status<=1) {
      sign = ((*p)->role==1 or (*p)->role==2) ? -1 : 1;
      pxtot += sign*(*p)->Px();
      pytot += sign*(*p)->Py();
      pztot += sign*(*p)->Pz();
      etot += sign*(*p)->E();
    }
  }
  // We set a threshold to the computation precision
  if (fabs(pxtot)<1.e-12) pxtot = 0.;
  if (fabs(pytot)<1.e-12) pytot = 0.;
  if (fabs(pztot)<1.e-12) pztot = 0.;
  if (fabs(etot)<1.e-12) etot = 0.;
  //
  std::cout << std::setfill('-') << std::setw(103) << "" << std::endl
	    << std::setfill(' ') << "Total:\t\t\t\t\t\t\t\t"
	    << std::right
	    << std::setprecision(2) << std::setw(9) << pxtot << " "
	    << std::setprecision(2) << std::setw(9) << pytot << " "
	    << std::setprecision(2) << std::setw(9) << pztot << " "
	    << std::setprecision(2) << std::setw(9) << etot << " "
	    << std::left << std::endl;
}
