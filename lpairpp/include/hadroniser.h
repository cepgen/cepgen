#include <string>
#include <vector>
#include "particle.h"

/**
 * Class template to define any hadroniser as a general object with defined methods
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date January 2014
 */
class Hadroniser
{
 public:
  Hadroniser();
  ~Hadroniser();
  /**
   * @brief Main caller to hadronise a particle
   */
  bool Hadronise(Particle *part_);
  inline std::vector<Particle>* GetHadrons() { return this->_hadrons; };
 protected:
  std::string _name;
  std::vector<Particle> *_hadrons;
};
