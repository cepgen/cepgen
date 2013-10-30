#include "particle.h"

Particle::Particle() :
  pdgId(-1), role(-1),
  e(-1.), m(-1.),
  px(0.), py(0.), pz(0.), pt(-1.), eta(0.),
  isValid(false), isPrimary(true)
{
  //this->SetMother(new Particle());
}

Particle::Particle(int role_, int pdgId_) :
  role(-1),
  e(-1.), m(-1.),
  px(0.), py(0.), pz(0.), pt(-1.), eta(0.),
  isValid(false), isPrimary(true)
{
  this->role = role_;
  this->pdgId = pdgId_;
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
Particle::SetP(double p_[3], double E_=-1.)
{
  if (E_<0.) {
    this->SetP(p_[0], p_[1], p_[2]);
  }
  else {
    this->SetP(p_[0], p_[1], p_[2], E_);
  }
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
