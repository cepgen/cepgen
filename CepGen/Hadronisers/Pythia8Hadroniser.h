#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include "GenericHadroniser.h"

#ifdef PYTHIA8
#include <Pythia8/Pythia.h>
#include <memory>
#endif

namespace CepGen
{
  class Parameters;
  namespace Hadroniser
  {
    /**
     * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia8 hadronisation algorithm
     */
    class Pythia8Hadroniser : public GenericHadroniser
    {
      public:
        explicit Pythia8Hadroniser( const Parameters& );
        ~Pythia8Hadroniser();

        bool decay( Event& ev, double& weight ) override;
        bool hadronise( Event& ev, double& weight ) override;
        void setSeed( long long seed ) override;

#ifdef PYTHIA8
        bool init();
        void readString( const char* param );
        void readString( const std::string& param ) { readString( param.c_str() ); }
#endif

      private:
        static constexpr unsigned short invalid_idx_ = 999;
        unsigned short max_attempts_;
        std::vector<unsigned short> min_ids_;
        std::map<short,short> py_cg_corresp_, cg_py_corresp_;
#ifdef PYTHIA8
        bool launchPythia( Event& ev );
        void fragmentState( unsigned short idx, double xbj = 0. );
        void updateEvent( Event& ev, double& weight );
        /// A Pythia8 core to be wrapped
        std::unique_ptr<Pythia8::Pythia> pythia_;
#endif
    };
  }
}

#endif
