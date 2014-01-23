#include "event.h"

Event::Event()
{
}

Event::~Event()
{
}

Event&
Event::operator=(const Event &ev_)
{
  this->_part = ev_._part;
  return *this;
}

std::vector<Particle*>
Event::GetByRole(int role_)
{
  int i;
  std::vector<Particle*> out;
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
  std::vector<Particle*> part_with_same_role = this->GetByRole(part_->role);
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
  std::multimap<int,Particle>::iterator p;

  std::stringstream ss;

  ss << "<event>" << std::endl;
  ss << this->NumParticles() << "\t0\t0.2983460E-04\t0.9118800E+02\t0.7546772E-02\t0.1300000E+00" << std::endl;
  for (p=this->_part.begin(); p!=this->_part.end(); p++) {
    if (p->second.status==0) p->second.status = 1;
    ss << std::setw(4) << p->second.id+1 << "\t"
       << std::setw(4) << p->second.status << "\t"
       << std::setw(4) << p->second.pdgId << "\t";
    if (p->second.GetMother()!=(Particle*)NULL) {
      ss << std::setw(4) << p->second.GetMother()->id+1 << "\t";
    }
    else {
      ss << std::setw(4) << "0" << "\t";
    }
    ss << std::setw(4) << "0" << "\t";
    ss << std::setw(4) << "0" << "\t";
    ss << std::setw(4) << "0" << "\t";
    ss << std::setw(8) << p->second.px << "\t"
       << std::setw(8) << p->second.py << "\t"
       << std::setw(8) << p->second.pz << "\t"
       << std::setw(8) << p->second.E() << "\t"
       << std::setw(8) << p->second.M() << "\t"
       << std::setw(8) << weight_
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

std::vector<Particle*>
Event::GetParticles()
{
  std::vector<Particle*> out;
  std::multimap<int,Particle>::iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    out.push_back(&it->second);
  }
  return out;
}

std::vector<Particle*>
Event::GetStableParticles()
{
  std::vector<Particle*> out;
  std::multimap<int,Particle>::iterator it;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    if (it->second.status==1) {
      out.push_back(&it->second);
    }
  }
  return out;  
}

void
Event::Dump(bool stable_)
{
  std::multimap<int,Particle>::iterator it;
  std::cout << "[Event::Dump]" << std::endl;
  std::cout << "Particle" << "\t" << "PDG id" << "\t\t" << "Charge" << "\t" << "Role" << "\t" << "Status" << "\t" << "Mother" << std::endl;
  std::cout << "--------" << "\t" << "------" << "\t\t" << "------" << "\t" << "----" << "\t" << "------" << "\t" << "------" << std::endl;
  for (it=this->_part.begin(); it!=this->_part.end(); it++) {
    if (stable_ and it->second.status!=1) continue;
    std::cout << std::setw(8) << it->second.id
	      << "\t" << std::setw(6) << it->second.pdgId;
    if (it->second.name!="")
      std::cout << " (" << it->second.name << ")";
    else std::cout << "\t";
    if (it->second.charge!=999.)
      std::cout << "\t" << std::setprecision(2) << std::setw(6) << it->second.charge;
    else std::cout << "\t";
    std::cout << "\t" << std::setw(4) << it->second.role
	      << "\t" << std::setw(6) << it->second.status;
    if (it->second.GetMother()!=(Particle*)NULL)
      std::cout << "\t" << std::setw(6) << it->second.GetMother()->id;
    std::cout << std::endl;
  }
}
