#include "particle.h"

Particle::Particle() :
  pdgId(-1), role(-1), e(-1.), m(-1.), px(0.), py(0.), pz(0.), pt(-1.)
{}

Particle::~Particle() {}

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
