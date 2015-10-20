#ifndef _EVENT_H
#define _EVENT_H

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "particle.h"

/**
 * @brief Convention to simplify the user interface while fetching a list of particles in the event
 */
typedef std::vector<Particle> Particles;
typedef std::vector<Particle*> ParticlesRef;
typedef std::multimap<int,Particle> ParticlesMap;

/**
 * Class containing all the information on the in- and outgoing particles' kinematics
 * @brief Kinematic information on the particles in the event
 */
class Event {
  public:
    Event();
    ~Event();
    /**
     * @brief Copies all the relevant quantities from one Event object to another
     */
    Event& operator=(const Event&);
    /**
     * @brief Empties the whole event content
     */
    inline void clear() { fParticles.clear(); time_generation=-1.; time_total=-1.; };
    /**
     * Returns the list of pointers to the Particle objects corresponding to a certain role in the process kinematics
     * @brief Gets a list of particles by their role in the event
     * @param[in] role_ The role the particles have to play in the process
     * @return A vector of pointers to the requested Particle objects
     */
    ParticlesRef GetByRole(int role_);
    /**
     * Returns the first Particle object in the particles list whose role corresponds to the given argument
     * @param[in] role_ The role the particle has to play in the event
     * @return A Particle object corresponding to the first particle found in this event
     */
    inline Particle* GetOneByRole(int role_) { 
      ParticlesRef out = GetByRole(role_);
      if (out.size()==0) return 0;
      else return out.at(0);
    };
    /**
     * Returns the pointer to the Particle object corresponding to a unique identifier in the event
     * @brief Gets one particle by its unique identifier in the event
     * @param[in] id_ The unique identifier to this particle in the event
     * @return A pointer to the requested Particle object
     */
    Particle* GetById(int id_);
    const Particle GetConstById(int id_) const;
    /**
     * Returns the pointers to the Particle objects corresponding to the unique identifiers in the event
     * @brief Gets a vector of particles by their unique identifier in the event
     * @param[in] ids_ The unique identifiers to the particles to be selected in the event
     * @return A vector of pointers to the requested Particle objects
     */
    inline ParticlesRef GetByIds(std::vector<int> ids_) {
      std::vector<int>::iterator id;
      ParticlesRef out;
      for (id=ids_.begin(); id!=ids_.end(); id++) out.push_back(GetById(*id));
      return out;
    }
    /**
     * Returns the pointer to the mother particle of any given Particle object in this event
     * @param[in] part_ The pointer to the Particle object from which we want to extract the mother particle
     * @return A pointer to the mother Particle object
     */
    inline ParticlesRef GetMothers(Particle* part_) {
      ParticlesRef out;
      const ParticlesIds moth = part_->GetMothersIds();
      ParticlesIds::iterator m;
      for (m=moth.begin(); m!=moth.end(); m++) {
      	out.push_back(GetById(*m));
      }
      return out;
    }; // FIXME
    inline Particles GetConstMothers(Particle* part_) const {
      Particles out;
      const ParticlesIds moth = part_->GetMothersIds();
      ParticlesIds::iterator m;
      for (m=moth.begin(); m!=moth.end(); m++) {
      	out.push_back(GetConstById(*m));
      }
      return out;
    }; // FIXME
    /**
     * @brief Gets a vector containing all the daughters from a particle
     * @param[in] part_ The particle for which the daughter particles have to be retrieved
     * @return A Particle objects vector containing all the daughters' kinematic information
     */
    inline ParticlesRef GetDaughters(Particle* part_) { return GetByIds(part_->GetDaughters()); };
    /**
     * Gets a list of roles for the given event (really process-dependant for the central system)
     * @return A vector of integers corresponding to all the roles the particles can play in the event
     */
    std::vector<int> GetRoles() const;
    /**
     * Sets the information on one particle in the process
     * @brief Add a particle to the event
     * @param[in] part_ The Particle object to insert or modify in the event
     * @param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
     * @return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int AddParticle(Particle part_, bool replace_=false);
    /**
     * @brief Creates a new particle in the event, with no kinematic information but the role it has to play in the process
     * @param[in] role_ The role the particle will play in the process
     * @param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
     * @return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int AddParticle(int role_, bool replace_=false);
    //HEPEUP GetHEPEUP() const;
    /**
     * Returns an event block in a LHE format (a XML-style) with all the information on the particles composing this event
     * @brief Gets the LHE block for this event
     * @param[in] weight_ The weight of the event
     * @return A string containing the kinematic quantities for each of the particles in the event, formatted as the LHE standard requires.
     */
    //std::string GetLHERecord(const double weight_=1.);
    /**
     * Stores in a file (raw format) all the kinematics on the outgoing leptons
     * @param[in] weight_ The weight of the event
     */
    void Store(std::ofstream*,double weight_=1.);
    //void Hadronise(std::string algo_="");
    /**
     * Dumps all the known information on every Particle object contained in this Event container in the output stream
     * @param[in] stable_ Do we only show the stable particles in this event ?
     */
    void Dump(bool stable_=false);
    /**
     * @brief Gets a vector of particles in the event
     * @return A vector containing all the pointers to the Particle objects contained in the event
     */
    ParticlesRef GetParticles();
    Particles GetConstParticles() const;
    /**
     * @brief Gets a vector of stable particles in the event
     * @return A vector containing all the pointers to the stable Particle objects contained in the event
     */
    ParticlesRef GetStableParticles();
    /**
     * @brief Number of particles in the event
     * @return The number of particles in the event, as an integer
     */
    inline int NumParticles() const { return fParticles.size(); };
    /**
     * @brief Number of trials before the event was "correctly" hadronised
     */
    int num_hadronisation_trials;
    /**
     * The time took by the generator to build the event without hadronising it, in seconds
     * @brief Time needed to generate the event at parton level
     */
    float time_generation;
    /**
     * The time took by the generator to build and hadronise the event, in seconds
     * @brief Time needed to generate the hadronised (if needed) event
     */
    float time_total;
    //HEPEUP event_info;
  private:
    /**
     * List of particles in the event, mapped to their role in this event
     */
    ParticlesMap fParticles;
    /**
     * Empty particle returned to the get-ers if no particle matches the requirements
     */
    Particle* np;
};

#endif
