#ifndef CepGen_Hadronisers_Pythia8Hadroniser_h
#define CepGen_Hadronisers_Pythia8Hadroniser_h

#include "CepGen/Hadronisers/GenericHadroniser.h"

#ifdef PYTHIA8
#include <Pythia8/Pythia.h>
#include <memory>
#endif

#include <unordered_map>
#include <vector>

namespace Pythia8 { class CepGenEvent; }
namespace cepgen
{
  class Particle;
  class Parameters;
  class ParametersList;
  enum class KinematicsMode;

  namespace hadr
  {
    /**
     * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia8 hadronisation algorithm
     */
    class Pythia8Hadroniser : public GenericHadroniser
    {
      public:
        explicit Pythia8Hadroniser( const ParametersList& );
        ~Pythia8Hadroniser();

        void setParameters( const Parameters& ) override;
        void readString( const char* param ) override;
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override;

      private:
        std::vector<unsigned short> min_ids_;
        std::unordered_map<short,short> py_cg_corresp_;
#ifdef PYTHIA8
        unsigned short findRole( const Event& ev, const Pythia8::Particle& p ) const;
        void updateEvent( Event& ev, double& weight ) const;
        Particle& addParticle( Event& ev, const Pythia8::Particle&, const Pythia8::Vec4& mom, unsigned short ) const;
        /// A Pythia8 core to be wrapped
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::unique_ptr<Pythia8::CepGenEvent> cg_evt_;
#endif
        static constexpr unsigned short PYTHIA_STATUS_IN_BEAM = 12;
        static constexpr unsigned short PYTHIA_STATUS_IN_PARTON_KT = 61;
        bool correct_central_;
        bool enable_hadr_;
        unsigned short offset_;
        bool first_evt_;
    };
  }
}

#endif
