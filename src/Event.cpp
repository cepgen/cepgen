#include "Event.h"

Event::Event() :
  num_hadronisation_trials(0),
  time_generation(-1.), time_total(-1.)
{
  //fParticles = new ParticlesMap();
}

Event::~Event()
{
  //delete fParticles;
}

Event&
Event::operator=(const Event &ev_)
{
  fParticles = ev_.fParticles;
  time_generation = ev_.time_generation;
  time_total = ev_.time_total;
  num_hadronisation_trials = ev_.num_hadronisation_trials;
  return *this;
}

void
Event::clear()
{
  fParticles.clear();
  time_generation = -1.;
  time_total = -1.;
}

void
Event::Init()
{
  fLastParticle = fParticles.end();
}

void
Event::Restore()
{
  fParticles.erase(fLastParticle, fParticles.end());
}

ParticlesRef
Event::GetByRole(int role_)
{
  int i;
  ParticlesRef out;
  std::pair<ParticlesMap::iterator,ParticlesMap::iterator> ret = fParticles.equal_range(role_);
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
  for (out=fParticles.begin(); out!=fParticles.end(); out++) {
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
  for (out=fParticles.begin(); out!=fParticles.end(); out++) {
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
  for (it=fParticles.begin(), end=fParticles.end(); it!=end; it=fParticles.upper_bound(it->first)) {
    out.push_back(it->first);
  }
  return out;
}

int
Event::AddParticle(Particle part_, bool replace_)
{
  DebugInsideLoop(Form("Particle with PDGid = %d has role %d", part_.GetPDGId(), part_.role));
  
  if (part_.role<=0) {
    return -1;
  }
  ParticlesRef part_with_same_role = GetByRole(part_.role);
  part_.id = fParticles.size(); //FIXME is there any better way of introducing this id ?
  if (replace_ and part_with_same_role.size()!=0) {
    part_with_same_role.at(0) = &part_;
    return 0;
  }
  fParticles.insert(std::pair<int,Particle>(part_.role, part_));
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
  out = AddParticle(*np, replace_);

  delete np;
  return out;
}

void
Event::Store(std::ofstream *of_, double weight_)
{
  Particle *l1 = GetByRole(6).at(0);
  Particle *l2 = GetByRole(7).at(0);
  
  *of_ << std::setw(8) << l1->E() << "\t"
       << std::setw(8) << l1->Px() << "\t"
       << std::setw(8) << l1->Py() << "\t"
       << std::setw(8) << l1->Pz() << "\t"
       << std::setw(8) << l1->Pt() << "\t"
       << std::setw(8) << l1->M() << "\t"
       << std::setw(8) << l1->Eta() << "\t"
       << std::setw(8) << l1->GetPDGId() << "\t"
       << std::setw(8) << weight_
       << std::endl;
  *of_ << std::setw(8) << l2->E() << "\t"
       << std::setw(8) << l2->Px() << "\t"
       << std::setw(8) << l2->Py() << "\t"
       << std::setw(8) << l2->Pz() << "\t"
       << std::setw(8) << l2->Pt() << "\t"
       << std::setw(8) << l2->M() << "\t"
       << std::setw(8) << l2->Eta() << "\t"
       << std::setw(8) << l2->GetPDGId() << "\t"
       << std::setw(8) << weight_
       << std::endl;
}

ParticlesRef
Event::GetParticles()
{
  ParticlesRef out;
  ParticlesMap::iterator it;
  for (it=fParticles.begin(); it!=fParticles.end(); it++) {
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
  for (it=fParticles.begin(); it!=fParticles.end(); it++) {
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
  for (it=fParticles.begin(); it!=fParticles.end(); it++) {
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
  std::ostringstream os;

  pxtot = pytot = pztot = etot = 0.;
  particles = GetParticles();
  for (p=particles.begin(); p!=particles.end(); p++) {
    if (stable_ and (*p)->status!=1) continue;
    os << Form("\n %2d\t%6d", (*p)->id, (*p)->GetPDGId());
    if ((*p)->name!="") os << Form("%6s", (*p)->name.c_str());
    //else                os << std::setw(6) << Particle::ParticleCode(abs((*p)->GetPDGId()));
    else                os << "\t";
    os << "\t";
    if ((*p)->charge!=999.) os << Form("%6.2f\t", (*p)->charge);
    else                    os << "\t";
    os << Form("%4d\t%6d\t", (*p)->role, (*p)->status);
    if ((*p)->GetMothersIds().size()>0) 
      os << Form("%2d(%2d)", *((*p)->GetMothersIds().begin()), GetById(*((*p)->GetMothersIds().begin()))->role);
    else
      os << "      ";
    os << Form("% 9.3f % 9.3f % 9.3f % 9.3f", (*p)->Px(), (*p)->Py(), (*p)->Pz(), (*p)->E());
    if ((*p)->status==0 or (*p)->status==1) {
      sign = ((*p)->status==0) ? -1 : 1;
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
  Info(Form("Events content:\n"
  "Part.\tPDG id\t\tCharge\tRole\tStatus\tMother\t\t4-Momentum [GeV]\n"
  "----\t------\t\t------\t----\t------\t------\t-------------------------------------"
  "%s\n"
  "---------------------------------------------------------------------------------------------\n"
  "Total:\t\t\t\t\t\t      % 9.4f % 9.4f % 9.4f % 9.4f", os.str().c_str(), pxtot, pytot, pztot, etot));
}
