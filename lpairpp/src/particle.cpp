#include "particle.h"

Particle::Particle() :
  pdgId(-1), role(-1),
  px(0.), py(0.), pz(0.), pt(-1.), eta(0.),
  isValid(false),
  e(-1.), m(-1.),
  isPrimary(true)
{
  //this->SetMother(new Particle());
}

Particle::Particle(int role_, int pdgId_) :
  role(-1),
  px(0.), py(0.), pz(0.), pt(-1.), eta(0.),
  isValid(false),
  e(-1.), m(-1.),
  isPrimary(true)
{
  this->role = role_;
  this->pdgId = pdgId_;
  if (this->pdgId!=0) {
    this->M(-1.);
  }
  //this->_mother = new Particle();
}

Particle::~Particle() {
  //delete this->_mother;
}

std::string Particle::GetLHEline(bool revert_)
{
  std::stringstream line;
  /*line << pdgId << "\t"
       << px << "\t"
       << py << "\t"
       << pz << "\t"
       << e << "\t"
       << m << std::endl;*/
  line << pdgId << "\t";
  line << "1 1 2 0 0" << "\t";
  line << px << "\t";
  line << py << "\t";
  if (revert_) {
    pz = -pz;
  }
  line << pz << "\t";
  line << e << "\t";
  line << m << "\t";
  line << "0." << "\t";
  line << "0."; //FIXME iz!!!
  return line.str();
}

bool
Particle::P(double p_[3], double E_=-1.)
{
  if (E_<0.) return this->P(p_[0], p_[1], p_[2]);
  else return this->P(p_[0], p_[1], p_[2], E_);
}

bool
Particle::M(double m_)
{
  double mass;
  if (m_>=0.) this->m = m_;
  else if (this->pdgId!=0) {
    mass = GetMassFromPDGId(this->pdgId);
    if (mass<0.) return false;
    this->m = mass;
  }
  else return false;
  return true;
}

Particle*
Particle::GetMother()
{
  if (!this->isPrimary) return this->_mother;
  return (Particle*)NULL;
}

void
Particle::Dump()
{
  if (this->isValid) {
    std::cout << "[Particle] DUMP"
	      << "\n\tRole = " << this->role
	      << "\n\tPDG id = " << this->pdgId
	      << "\n\tP = (" << this->px << ", " << this->py << ", " << this->pz << ") GeV"
	      << "\n\tPt = " << this->pt << " GeV"
	      << "\n\tE = " << this->e << " GeV"
	      << "\n\tM = " << this->m << " GeV"
	      << "\n\teta = " << this->eta
	      << "\n\tIs valid ? " << this->isValid << std::endl;
    if (!this->isPrimary) {
      std::cout << "\tMother = " << this->GetMother()->role 
		<< " (pdgId=" << this->GetMother()->pdgId << ")"
		<< std::endl;
    }
  }
}
