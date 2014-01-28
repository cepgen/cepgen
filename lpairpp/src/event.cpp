#include "event.h"

Event::Event()
{
  //this->_part = new std::multimap<int,Particle>();
}

Event::~Event()
{
  //delete this->_part;
}

/*Event&
Event::operator=(const Event &ev_)
{
  this->_part = ev_._part;
  return *this;
  }*/

Particles
Event::GetByRole(int role_)
{
  int i;
  Particles out;
  std::pair<std::multimap<int,Particle>::iterator,std::multimap<int,Particle>::iterator> ret = this->_part.equal_range(role_);
  std::multimap<int,Particle>::iterator it;

  for (it=ret.first, i=0; it!=ret.second && i<100; it++, i++) {
    out.push_back(&(it->second));
  }
  return out;
}

Particle*
Event::GetById(int id_)
{
  std::multimap<int,Particle>::iterator out;
  for (out=this->_part.begin(); out!=this->_part.end(); out++) {
    if (out->second.id==id_) {
      return &out->second;
    }
  }
  return (Particle*)NULL;
}

std::vector<int>
Event::GetRoles()
{
  std::multimap<int,Particle>::iterator it, end;
  std::vector<int> out;
  for (it=this->_part.begin(), end=this->_part.end(); it!=end; it=this->_part.upper_bound(it->first)) {
    out.push_back(it->first);
  }
  return out;
}

int
Event::AddParticle(Particle *part_, bool replace_)
{
#ifdef DEBUG
  std::cout << "[Event::AddParticle] [DEBUG] Particle with PDGid = " << part_->pdgId << " has role " << part_->role << std::endl;
#endif
  if (part_->role<=0) {
    return -1;
  }
  Particles part_with_same_role = this->GetByRole(part_->role);
  part_->id = this->_part.size(); //FIXME is there any better way of introducing this id ?
  if (replace_ and part_with_same_role.size()!=0) {
    part_with_same_role.at(0) = part_;
    return 0;
  }
  this->_part.insert(std::pair<int,Particle>(part_->role, *part_));
  return 1;
}

std::string
Event::GetLHERecord(const double weight_)
{
  //std::multimap<int,Particle>::iterator p;
  std::stringstream ss;
  Particles particles, daughters;
  Particles::iterator p, dg;
  int min_id, max_id;

  //FIXME need to fetch the vector (not the multimap), so that we can sort on the particle unique identifier (also TODO!!!)
  ss << "<event>" << std::endl;
  /*ss << this->NumParticles() << "\t"
     << this->event_info.idprup << "\t"
     << this->event_info.xwgtup << "\t"
     << this->event_info.scalup << "\t" // scale of the event, in GeV
     << this->event_info.aqedup << "\t" // alphaQED
     << this->event_info.aqcdup // alphaQCD
     << std::endl;*/
  particles = this->GetParticles();
  for (p=particles.begin(); p!=particles.end(); p++) {
    if ((*p)->status==0) (*p)->status = 1;
    
    ss << std::setw(4) << (*p)->id+1 << "  "
       << std::setw(3) << (*p)->status << "  "
       << std::setw(5) << (*p)->pdgId << "  ";
    if (!(*p)->Primary()) {
      ss << std::setw(2) << (*p)->GetMother()+1 << "  ";
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
	ss << std::setw(4) << min_id << "  " << std::setw(4) << "0";
      else 
	ss << std::setw(4) << min_id << "  " << std::setw(4) << max_id;
    }
    else {
      ss << std::setw(4) << "0" << "  " << std::setw(4) << "0";
    }
    ss << std::setw(4) << "  " << "0" << "  "; 
    ss << std::setw(12) << (*p)->px << "  " 
       << std::setw(12) << (*p)->py << "  " 
       << std::setw(12) << (*p)->pz << "  " 
       << std::setw(12) << (*p)->E() << "  " 
       << std::setw(12) << (*p)->M() << "  " 
       << std::setw(4) << weight_ 
       << std::endl; 
  }
  ss << "</event>" << std::endl;

  return ss.str();
}

void
Event::Store(std::ofstream *of_, double weight_)
{
  Particle *l1 = this->GetByRole(6).at(0);
  Particle *l2 = this->GetByRole(7).at(0);
  
  *of_ << std::setw(8) << l1->E() << "\t"
       << std::setw(8) << l1->px << "\t"
       << std::setw(8) << l1->py << "\t"
       << std::setw(8) << l1->pz << "\t"
       << std::setw(8) << l1->Pt() << "\t"
       << std::setw(8) << l1->M() << "\t"
       << std::setw(8) << l1->Eta() << "\t"
       << std::setw(8) << l1->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
  *of_ << std::setw(8) << l2->E() << "\t"
       << std::setw(8) << l2->px << "\t"
       << std::setw(8) << l2->py << "\t"
       << std::setw(8) << l2->pz << "\t"
       << std::setw(8) << l2->Pt() << "\t"
       << std::setw(8) << l2->M() << "\t"
       << std::setw(8) << l2->Eta() << "\t"
       << std::setw(8) << l2->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
}

Particles
Event::GetParticles()
{
  Particles out;
  std::multimap<int,Particle>::iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    out.push_back(&it->second);
  }
  std::sort(out.begin(), out.end());
  return out;
}

Particles
Event::GetStableParticles()
{
  Particles out;
  std::multimap<int,Particle>::iterator it;
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
  //std::multimap<int,Particle>::iterator it;
  Particles particles;
  Particles::iterator p;

  particles = this->GetParticles();
  std::cout << "[Event::Dump]" << std::endl;
  std::cout << "Particle" << "\t" << "PDG id" << "\t\t" << "Charge" << "\t" << "Role" << "\t" << "Status" << "\t" << "Mother" << std::endl;
  std::cout << "--------" << "\t" << "------" << "\t\t" << "------" << "\t" << "----" << "\t" << "------" << "\t" << "------" << std::endl;
  for (p=particles.begin(); p!=particles.end(); p++) {
    if (stable_ and (*p)->status!=1) continue;
    std::cout << std::setw(8) << (*p)->id
	      << "\t" << std::setw(6) << (*p)->pdgId;
    if ((*p)->name!="")
      std::cout << " " << std::setw(6) << (*p)->name;
    else
      std::cout << "\t";
    std::cout << "\t";
    if ((*p)->charge!=999.)
	std::cout << std::setprecision(2) << std::setw(6) << (*p)->charge;
    std::cout << "\t" << std::setw(4) << (*p)->role
	      << "\t" << std::setw(6) << (*p)->status;
      //<< std::endl;
    if ((*p)->GetMother()!=-1)
      std::cout << "\t" << std::setw(6) << (*p)->GetMother();
    std::cout << std::endl;
  }
}
