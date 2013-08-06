#ifndef _LHE_H
#define _LHE_H

#include "gamgam.h"

/**
 * @brief Stores the event according to the Les Houches Event format
 * @param symm_ Is the symmetrisation of the event necessary ? If set to true,
 * the third component of the outgoing leptons' momentum is reverted to simulate
 * a symmetric proton-proton or electron-electron collision.
 * @param g_ The GamGam object related to the events generation requested.
 */
std::string GetLHEvent(GamGam*,bool symm_=false);


#endif
