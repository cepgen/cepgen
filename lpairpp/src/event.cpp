#include "event.h"

Event::Event()
{
  this->_part = new std::map<int,Particle>;
  this->_null = new Particle();
}

Event::~Event()
{
  delete this->_part;
  delete this->_null;
}

Particle*
Event::GetByRole(int role_)
{
  std::map<int,Particle>::iterator out = this->_part->find(role_);
  if (out!=this->_part->end()) {
    return &(out->second);
  }
  else {
    return this->_null;
  }
}

int
Event::SetParticle(Particle *part_)
{
#ifdef DEBUG
  std::cout << "[Event::SetParticle] [DEBUG] Particle with PDGid = " << part_->pdgId << " has role " << part_->role << std::endl;
#endif
  if (part_->role<=0) {
    return -1;
  }
  Particle *tmp = this->GetByRole(part_->role);
  if (!tmp->isValid) {
    this->_part->insert(std::pair<int,Particle>(part_->role, *part_));
    return 0;
  }
  else {
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
         << std::setw(8) << p->second.pt << "\t"
         << std::setw(8) << p->second.M() << "\t"
         << std::setw(8) << p->second.eta << "\t"
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
       << std::setw(8) << l1->pt << "\t"
       << std::setw(8) << l1->M() << "\t"
       << std::setw(8) << l1->eta << "\t"
       << std::setw(8) << l1->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
  *of_ << std::setw(8) << l2->E() << "\t"
       << std::setw(8) << l2->px << "\t"
       << std::setw(8) << l2->py << "\t"
       << std::setw(8) << l2->pz << "\t"
       << std::setw(8) << l2->pt << "\t"
       << std::setw(8) << l2->M() << "\t"
       << std::setw(8) << l2->eta << "\t"
       << std::setw(8) << l2->pdgId << "\t"
       << std::setw(8) << weight_
       << std::endl;
}

void
Event::Dump()
{
  std::map<int,Particle>::iterator out;
  std::cout << "[Event::Dump]" << std::endl;
  for (out=this->_part->begin(); out!=this->_part->end(); out++) {
    /*std::cout << "-> Particle with role " << out->second.role
              << "\n\tPDG = " << out->second.pdgId
              << "\n\tP = (" << out->second.px << ", " << out->second.py << ", " << out->second.pz << ")"
              << "\n\tE = " << out->second.e
              << "\n\tM = " << out->second.m
              << std::endl;*/
    out->second.Dump();
    std::cout << "=========================" << std::endl;
  }
}
