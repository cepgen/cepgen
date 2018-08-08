#ifndef CepGen_Hadronisers_Pythia8Hadroniser_h
#define CepGen_Hadronisers_Pythia8Hadroniser_h

#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Physics/Kinematics.h"

#ifdef PYTHIA8
#include <Pythia8/Pythia.h>
#include <memory>
#endif

#include <unordered_map>
#include <vector>

namespace CepGen
{
  class Particle;
#ifdef PYTHIA8
  class LHAEvent : public Pythia8::LHAup
  {
    public:
      explicit LHAEvent( const Parameters* );
      void feedEvent( const Event& ev, bool full, const Kinematics::Mode& );
      bool setInit() override;
      bool setEvent( int ) override;
      void setCrossSection( int id, double xsec, double xsec_err );
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

      unsigned short cgPart( unsigned short py_id ) const;
      unsigned short pyPart( unsigned short cg_id ) const;
      void addCorresp( unsigned short py_id, unsigned short cg_id );
      void dumpCorresp() const;

      static constexpr unsigned short invalid_id = 999;
    private:
      static const double mp_, mp2_;
      std::vector<std::pair<unsigned short, unsigned short> > py_cg_corresp_;
      const Parameters* params_;
  };
#endif

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

        bool run( Event& ev, double& weight, bool full ) override;
        void setSeed( long long seed ) override;
        void setCrossSection( double xsec, double xsec_err ) override;

        bool fullEvent() const { return full_evt_; }
        void setFullEvent( bool full = true ) { full_evt_ = full; }

        bool init();
        void readString( const char* param ) override;

      private:
        static constexpr unsigned short invalid_idx_ = 999;
        unsigned short max_attempts_;
        std::vector<unsigned short> min_ids_;
        std::unordered_map<short,short> py_cg_corresp_;
#ifdef PYTHIA8
        unsigned short findRole( const Event& ev, const Pythia8::Particle& p ) const;
        void updateEvent( Event& ev, double& weight, bool full ) const;
        Particle& addParticle( Event& ev, const Pythia8::Particle&, const Pythia8::Vec4& mom, unsigned short ) const;
        /// A Pythia8 core to be wrapped
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::shared_ptr<LHAEvent> lhaevt_;
#endif
        bool full_evt_;
        unsigned short offset_;
        bool first_evt_;
        const Parameters* params_; // not owning
    };
  }
}

#endif

