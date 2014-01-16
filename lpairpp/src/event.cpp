#include "event.h"

Event::Event()
{
  this->_null = new Particle();
  this->_part = new std::map<int,Particle>();
}

Event::~Event()
{
  delete this->_null;
  delete this->_part;
}

Particle*
Event::GetByRole(int role_)
{
  std::map<int,Particle>::iterator out = this->_part->find(role_);
  if (out!=this->_part->end()) {
    //std::cout << "[Event::GetByRole] [DEBUG] Particle with role " << role_ << " successfully returned (pdgId=" << out->second.pdgId << ")" << std::endl;
    return &out->second;
  }
  else {
    //std::cout << "[Event::GetByRole] [DEBUG] Failed to locate particle with role " << role_ << " !" << std::endl;
    return this->_null;
  }
}

Particle*
Event::GetById(int id_)
{
  std::map<int,Particle>::iterator out;
  for (out=this->_part->begin(); out!=this->_part->end(); out++) {
    if (out->second.id==id_) {
      return &out->second;
    }
  }
  return this->_null;
}

int
Event::AddParticle(Particle *part_)
{
#ifdef DEBUG
  std::cout << "[Event::AddParticle] [DEBUG] Particle with PDGid = " << part_->pdgId << " has role " << part_->role << std::endl;
#endif
  if (part_->role<=0) {
    return -1;
  }
  if (!this->GetByRole(part_->role)->Valid()) {
    part_->id = this->_part->size(); //FIXME is there any better way of introducing this id ?
    this->_part->insert(std::pair<int,Particle>(part_->role, *part_));
    //this->GetByRole(part_->role)->id = this->_part->size();
    return 0;
  }
  else {
#ifdef DEBUG
    std::cout << "[Event::AddParticle] [DEBUG] Replacing an existing particle : " 
	      << part_->role << " (pdgId=" << part_->pdgId << ", p=" << part_->P() << ") --> " 
	      << tmp->role << " (pdgId=" << tmp->pdgId << ", p=" << part_->P() << ")" << std::endl;
#endif
    this->_part->at(part_->role) = *part_;
    return 1;
  }
}

void
Event::StoreLHERecord(std::ofstream *of_, const double weight_)
{
  std::map<int,Particle>::iterator p;
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
  Particle *l1 = this->GetByRole(6);
  Particle *l2 = this->GetByRole(7);
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
  std::map<int,Particle>::iterator it;
  for (it=this->_part->begin(); it!=this->_part->end(); it++) {
    out.push_back(&it->second);
    //std::cout << it->second.id << "\t" << it->second.role << "\t" << it->second.pdgId << "\t" << it->second.P() << std::endl;
    //std::cout << it->first << std::endl;
  }
  return out;
}

void
Event::Dump()
{
  std::map<int,Particle>::iterator out;
  std::cout << "[Event::Dump]" << std::endl;
  for (out=this->_part->begin(); out!=this->_part->end(); out++) {
    out->second.Dump();
    std::cout << "=========================" << std::endl;
  }
}
