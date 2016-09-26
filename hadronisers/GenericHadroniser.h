#ifndef GenericHadroniser_h
#define GenericHadroniser_h

#include "physics/Event.h"
#include "physics/Physics.h"

/**
 * Class template to define any hadroniser as a general object with defined methods
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date January 2014
 */
class GenericHadroniser
{
 public:
  /// Default constructor for an undefined hadroniser
  GenericHadroniser(std::string name_="unnamed_hadroniser");
  virtual ~GenericHadroniser();
  /// Main caller to hadronise a single particle
  /// \param[in] part_ The Particle object which will be hadronised
  /// \return A boolean stating whether or not the hadronisation occured successfully
  inline virtual bool Hadronise(Particle *part_) { return (part_!=(Particle*)NULL and part_->status!=2); };
  /// Hadronise a full event
  /// \param[inout] ev_ Event to hadronise
  /// \return Boolean stating whether or not the hadronisation occured successfully
  inline virtual bool Hadronise(Event *ev_) { ev_->Dump(); return false; };
  /// Get the full list of hadrons produced in the hadronisation
  /// \return Vector of Particle containing all the hadrons produced
  inline Particles GetHadrons() { return *fHadrons; };
  /// Return a human-readable name for this hadroniser
  inline std::string GetName() { return fName; };
 protected:
  /// Name of the hadroniser
  std::string fName;
  /// List of hadrons produced by this hadronisation process
  Particles *fHadrons;
};

#endif
