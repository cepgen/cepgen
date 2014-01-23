#include "event.h"

Event::Event()
{
  this->_null = new Particle();
  this->_part = new std::multimap<int,Particle>();
}

Event::~Event()
{
  delete this->_null;
  delete this->_part;
}

std::vector<Particle*>
Event::GetByRole(int role_)
{
  int i;
  std::vector<Particle*> out;
  std::pair<std::multimap<int,Particle>::iterator,std::multimap<int,Particle>::iterator> ret = this->_part->equal_range(role_);
  std::multimap<int,Particle>::iterator it;

  for (it=ret.first, i=0; it!=ret.second && i<100; it++, i++) {
    out.push_back(&(it->second));
  }
  return out;
  /*  if (ret!=this->_part->end()) {
#ifdef DEBUG
    std::cout << "[Event::GetByRole] [DEBUG] Particle with role " << role_ << " successfully returned (pdgId=" << ret->second.pdgId << ")" << std::endl;
#endif
    return &ret->second;
  }
  else {
#ifdef DEBUG
    std::cout << "[Event::GetByRole] [DEBUG] Failed to locate particle with role " << role_ << " !" << std::endl;
#endif
    return this->_null;
    }*/
}

Particle*
Event::GetById(int id_)
{
  std::multimap<int,Particle>::iterator out;
  for (out=this->_part->begin(); out!=this->_part->end(); out++) {
    if (out->second.id==id_) {
      return &out->second;
    }
  }
  return this->_null;
}

std::vector<int>
Event::GetRoles()
{
  std::multimap<int,Particle>::iterator it, end;
  std::vector<int> out;
  for (it=this->_part->begin(), end=this->_part->end(); it!=end; it=this->_part->upper_bound(it->first)) {
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
  part_->id = this->_part->size(); //FIXME is there any better way of introducing this id ?
  if (replace_ and part_with_same_role.size()!=0) {
    part_with_same_role.at(0) = part_;
    return 0;
  }
  this->_part->insert(std::pair<int,Particle>(part_->role, *part_));
  return 1;
}

void
Event::StoreLHERecord(std::ofstream *of_, const double weight_)
{
  std::multimap<int,Particle>::iterator p;
  for (p=this->_part->begin(); p!=this->_part->end(); p++) {
    *of_ << std::setw(8) << p->second.E() << "\t"
         << std::setw(8) << p->second.px << "\t"
         << std::setw(8) << p->second.py << "\t"
         << std::setw(8) << p->second.pz << "\t"
         << std::setw(8) << p->second.Pt() << "\t"
         << std::setw(8) << p->second.M() << "\t"
         << std::setw(8) << p->second.Eta() << "\t"
         << std::setw(8) << p->second.pdgId << "\t"
         << std::setw(8) << weight_
         << std::endl;
  }
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
  for (it=this->_part->begin(); it!=this->_part->end(); it++) {
    out.push_back(&it->second);
  }
  return out;
}

std::vector<Particle*>
Event::GetStableParticles()
{
  std::vector<Particle*> out;
  std::multimap<int,Particle>::iterator it;
  for (it=this->_part->begin(); it!=this->_part->end(); it++) {
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
  for (it=this->_part->begin(); it!=this->_part->end(); it++) {
    if (stable_ and it->second.status!=1) continue;
    std::cout << std::setw(8) << it->second.id
	      << "\t" << std::setw(6) << it->second.pdgId;
    if (it->second.name!="")
      std::cout << " (" << it->second.name << ")";
    else std::cout << "\t";
    if (it->second.charge!=999.)
      std::cout << "\t" << std::setprecision(2) << std::setw(6) << it->second.charge;
    else std::cout << "\t";
    //std::cout << std::endl;
    std::cout << "\t" << std::setw(4) << it->second.role// << std::endl
	      << "\t" << std::setw(6) << it->second.status;
      //<< std::endl;
    if (it->second.GetMother()!=(Particle*)NULL)
      std::cout << "\t" << std::setw(6) << it->second.GetMother()->id;
    std::cout << std::endl;
    //it->second.Dump();
    //std::cout << "=========================" << std::endl;
  }
}
