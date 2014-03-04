#ifndef _HADRONISER_H
#define _HADRONISER_H

#include <string>
#include <vector>

#include "event.h"
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
  inline virtual bool Hadronise(Particle *part_) { return (part_!=(Particle*)NULL && part_->status!=2); };
  /**
   * Launches the hadroniser on the full event information
   * @brief Hadronises a full event
   * @param[inout] ev_ The event to hadronise
   * @return A boolean stating whether or not the hadronisation occured successfully
   */
  inline virtual bool Hadronise(Event *ev_) { ev_->Dump(); return false; };
  /**
   * Gets the full list of hadrons (as Particle objects) produced by the hadronisation
   * @return A vector of Particle containing all the hadrons produced
   */
  inline std::vector<Particle> GetHadrons() { return *(this->_hadrons); };
  inline std::string GetName() { return this->_name; };
 protected:
  /** @brief Name of the hadroniser */
  std::string _name;
  /** @brief List of hadrons produced by this hadronisation process */
  std::vector<Particle> *_hadrons;
};

#endif
