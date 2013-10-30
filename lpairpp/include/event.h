#ifndef _EVENT_H
#define _EVENT_H

#include <map>
#include <string>
#include <fstream>
#include <iomanip>

#include "particle.h"

/**
 * Class containing all the information on the in- and outgoing particles' kinematics
 * @brief Kinematic information on the particles in the event
 */
class Event {
  public:
    Event();
    ~Event();
    /**
     * Returns the pointer to the Particle object corresponding to a certain role in the process kinematics
     * @param role_ The role the particle has to play in the process
     * @return A pointer to the requested Particle object
     */
    Particle *GetByRole(int role_);
    /**
     * Sets the information on one particle in the process
     * @param part_ The Particle object to insert or modify in the event
     * @return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int SetParticle(Particle* part_);
    /**
     * Stores in a LHE format (a XML-style) all the information on the particles composing this event
     * @brief Stores the LHE block for this event
     * @param of_ The file stream on which the event record has to be saved
     * @param weight_ The weight of the event
     */
    void StoreLHERecord(std::ofstream*, const double weight_=1.);
    /**
     * Stores in a file (raw format) all the kinematics on the outgoing leptons
     * @param weight_ The weight of the event
     */
    void Store(std::ofstream*,double weight_=1.);
    /**
     * Dumps all the known information on every Particle object contained in this Event container in the output stream
     */
    void Dump();
  private:
    std::map<int,Particle> *_part;
};

#endif
