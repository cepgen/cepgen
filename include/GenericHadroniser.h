#ifndef GenericHadroniser_h
#define GenericHadroniser_h

#include <string>
#include <vector>

#include "Event.h"
#include "Physics.h"

/**
 * Class template to define any hadroniser as a general object with defined methods
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date January 2014
 */
class GenericHadroniser
{
 public:
  GenericHadroniser(std::string name_="unnamed_hadroniser");
  virtual ~GenericHadroniser();
  /**
   * @brief Main caller to hadronise a single particle
   * @param[in] part_ The Particle object which will be hadronised
   * @return A boolean stating whether or not the hadronisation occured successfully
   */
  inline virtual bool Hadronise(Particle *part_) { return (part_!=(Particle*)NULL and part_->status!=2); };
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
  inline Particles GetHadrons() { return *fHadrons; };
  /**
   * @brief Returns the human-readable name of the hadroniser used
   */
  inline std::string GetName() { return fName; };
 protected:
  /**
   * @brief Name of the hadroniser
   */
  std::string fName;
  /**
   * @brief List of hadrons produced by this hadronisation process
   */
  Particles *fHadrons;
};

#endif
