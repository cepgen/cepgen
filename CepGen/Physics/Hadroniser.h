#ifndef CepGen_Physics_Hadroniser_h
#define CepGen_Physics_Hadroniser_h

#include "CepGen/Core/EventModifier.h"

namespace cepgen {
  /// Location for all hadronisers to be run downstream to the events generation
  namespace hadr {
    /**
     * \brief Class template to define any hadroniser as a general object with defined methods
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date January 2014
     */
    class Hadroniser : public EventModifier {
    public:
      /// Default constructor for an undefined hadroniser
      explicit Hadroniser(const ParametersList&);

      /// Specify whether the beam remnants are to be fragmented
      bool fragmentRemnants() const { return remn_fragm_; }

    protected:
      /// Switch on/off the remnants fragmentation where applicable
      const bool remn_fragm_;
    };
  }  // namespace hadr
}  // namespace cepgen

#endif
