#ifndef _EVENT_H
#define _EVENT_H

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

#include "particle.h"

//#include "pythia6hadroniser.h"

/**
 * Class containing all the information on the in- and outgoing particles' kinematics
 * @brief Kinematic information on the particles in the event
 */
class Event {
  public:
    Event();
    ~Event();
    /**
     * Returns the list of pointers to the Particle objects corresponding to a certain role in the process kinematics
     * @brief Gets a list of particles by their role in the event
     * @param role_ The role the particles have to play in the process
     * @return A vector of pointers to the requested Particle objects
     */
    std::vector<Particle*> GetByRole(int role_);
    Particle* GetOneByRole(int role_) { 
      std::vector<Particle*> out = this->GetByRole(role_);
      if (out.size()==0) return _null;
      else return out.at(0);
    };
    /**
     * Returns the pointer to the Particle object corresponding to a unique identifier in the event
     * @brief Gets one particle by its unique identifier in the event
     * @param id_ The unique identifier to this particle in the event
     * @return A pointer to the requested Particle object
     */
    Particle* GetById(int id_);
    std::vector<int> GetRoles();
    /**
     * Sets the information on one particle in the process
     * @brief Add a particle to the event
     * @param part_ The Particle object to insert or modify in the event
     * @param replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
     * @return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int AddParticle(Particle* part_, bool replace_=false);
    /**
     * Stores in a LHE format (a XML-style) all the information on the particles composing this event
     * @brief Stores the LHE block for this event
     * @param of_ The file stream on which the event record has to be saved
     * @param weight_ The weight of the event
     */
    void StoreLHERecord(std::ofstream* of_, const double weight_=1.);
    /**
     * Stores in a file (raw format) all the kinematics on the outgoing leptons
     * @param weight_ The weight of the event
     */
    void Store(std::ofstream*,double weight_=1.);
    //void Hadronise(std::string algo_="");
    /**
     * Dumps all the known information on every Particle object contained in this Event container in the output stream
     */
    void Dump();
    /**
     * @brief Gets a vector of particles in the event
     * @return A vector containing all the pointers to the Particle objects contained in the event
     */
    std::vector<Particle*> GetParticles();
    std::vector<Particle*> GetStableParticles();
    /**
     * @brief Number of particles in the event
     * @return The number of particles in the event, as an integer
     */
    inline int NumParticles() { return this->_part->size(); };
  private:
    std::multimap<int,Particle> *_part;
    Particle *_null;
};

#endif
