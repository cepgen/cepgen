#ifndef _INELASTIC_H
#define _INELASTIC_H

#include <vector>
#include <string>

#include "particle.h"
#include "pythia6hadroniser.h"

/**
 * Class containing the information on a particle supposed to decay or
 * fragment in the process
 */
class InelasticParticle : public Particle {
  public:
    InelasticParticle();
    ~InelasticParticle();
    void PDF2PDG();
    /**
     * Hadronises the particle with Pythia, and builds the shower
     * (list of Particle objects) embedded in this object
     * @param algo_ Algorithm in use to hadronise the particle
     * @brief Hadronises the particle using Pythia
     */
    bool Hadronise(std::string algo_);
  private:
    /**
     * @brief List of particles produced with this decay
     * Shower of particles arising from the decay of this
     * inelastic particle
     */
    std::vector<Particle> *_shower;
};

#endif
