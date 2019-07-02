#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/ParametersList.h" //FIXME
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <algorithm>

extern "C"
{
  /// Get the particle's mass in GeV from the Pythia6 module
  extern double pymass_( int& );
  /// Get the resonant particle's width in GeV from the Pythia6 module
  //extern double pywidt_( int& );
  /// Launch the Pythia6 fragmentation
  extern void pyexec_();
  /// Set a parameter value to the Pythia6 module
  extern void pygive_( const char*, int );
  extern void pyckbd_();
  /// List all the particles in the event in a human-readable format
  extern void pylist_( int& );
  /// Join two coloured particles in a colour singlet
  extern void pyjoin_( int&, int& );
  /// Get a particle's human-readable name from the Pythia6 module
  extern void pyname_( int&, char*, int );
  /// Get integer-valued event information from the Pythia6 module
  extern int pyk_( int&, int& );
  /// Get real-valued event information from the Pythia6 module
  extern double pyp_( int&, int& );
  /// Purely virtual method to call at the end of the run
  void pystop_() { CG_INFO( "Pythia6Hadroniser" ) << "End of run"; }

  /// Particles content of the event
  extern struct
  {
    /// Number of particles in the event
    int n;
    int npad;
    /// Particles' general information (status, PDG id, mother, daughter 1, daughter 2)
    int k[5][4000];
    /// Particles' kinematics, in GeV (px, py, pz, E, M)
    double p[5][4000];
    /// Primary vertex for the particles
    double v[5][4000];
  } pyjets_;
}

namespace cepgen
{
  namespace hadr
  {
    /**
     * Full interface to the Pythia6 @cite Sjostrand:2006za algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * @brief Pythia6 hadronisation algorithm
     */
    class Pythia6Hadroniser : public GenericHadroniser
    {
      public:
        Pythia6Hadroniser( const ParametersList& );

        void setParameters( const Parameters& ) override {}
        inline void readString( const char* param ) override { pygive( param ); }
        void init() override {}
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}
        //bool hadronise( const Particle* );

      private:
        static constexpr unsigned short MAX_PART_STRING = 3;
        static constexpr unsigned short MAX_STRING_EVENT = 2;
        /// Maximal number of characters to fetch for the particle's name
        static constexpr unsigned short NAME_CHR = 16;

        inline static double pymass(int pdgid_) { return pymass_(pdgid_); }
        //inline static double pywidt(int pdgid_) { return pywidt_(pdgid_); }
        inline static void pyckbd() { pyckbd_(); }
        inline static void pygive( const std::string& line ) { pygive_( line.c_str(), line.length() ); }
        inline static void pylist( int mlist ) { pylist_( mlist ); }
        inline static int pyk( int id, int qty ) { return pyk_( id, qty ); }
        inline static double pyp( int id, int qty ) { return pyp_( id, qty ); }
        inline static std::string pyname( int pdgid ) {
          char out[NAME_CHR];
          std::string s;
          pyname_( pdgid, out, NAME_CHR );
          s = std::string( out, NAME_CHR );
          s.erase( remove( s.begin(), s.end(), ' ' ), s.end() );
          return s;
        }
        /**
         * \brief Connect entries with colour flow information
         * \param[in] njoin Number of particles to join in the colour flow
         * \param[in] ijoin List of particles unique identifier to join in the colour flow
         */
        inline static void pyjoin( int njoin, int ijoin[2] ) { return pyjoin_( njoin, *ijoin ); }
        bool prepareHadronisation( Event& );
        std::pair<short,short> pickPartonsContent() const;
        struct EventProperties
        {
          unsigned int str_in_evt = 0;
          unsigned int num_part_in_str[MAX_STRING_EVENT] = { 0 };
        };
        EventProperties fillParticles( const Event& ) const;
    };

    Pythia6Hadroniser::Pythia6Hadroniser( const ParametersList& plist ) :
      GenericHadroniser( plist, "pythia6" )
    {}

    bool
    Pythia6Hadroniser::run( Event& ev, double& weight, bool full )
    {
      full = true; //FIXME

      weight = 1.;
      if ( full )
        prepareHadronisation( ev );

      if ( utils::Logger::get().level >= utils::Logger::Level::debug ) {
        CG_DEBUG( "Pythia6Hadroniser" ) << "Dump of the event before the hadronisation:";
        ev.dump();
      }

      //--- fill Pythia 6 common blocks
      auto prop = fillParticles( ev );

      CG_DEBUG( "Pythia6Hadroniser" )
        << "Passed the string construction stage.\n\t"
        << " " << prop.str_in_evt << " string objects were identified and constructed.";

      unsigned int oldnpart = pyjets_.n;

      pyexec_();

      int criteria = oldnpart+1;
      for ( unsigned int i = 0; i < MAX_STRING_EVENT; ++i )
        criteria += prop.num_part_in_str[i];

      if ( pyjets_.k[1][criteria] == 2212
        && pyjets_.k[0][criteria] == 1 ) {
        CG_WARNING( "Pythia6Hadroniser" ) << "System is non-inelastic.";
        return false;
      }

      // We filter the first particles already present in the event
      for ( unsigned int p = oldnpart; p < (unsigned int)pyjets_.n; ++p ) {
        const pdgid_t pdg_id = abs( pyjets_.k[1][p] );
        const short charge = pyjets_.k[1][p]/(short)pdg_id;
        ParticleProperties prop;
        if ( full )
          try { prop = PDG::get()( pdg_id ); } catch ( const Exception& ) {
            prop = ParticleProperties{ pdg_id,
              pyname( pdg_id ), pyname( pdg_id ),
              (short)pyk( p+1, 12 ), // colour factor
              pymass( pdg_id ), -1., //pmas( pdg_id, 2 ),
              (short)pyk( p+1, 6 ), // charge
              false
            };
            PDG::get().define( prop );
          }

        const Particle::Role role = pyjets_.k[2][p] != 0
          ? ev[pyjets_.k[2][p]-1].role() // child particle inherits its mother's role
          : Particle::Role::UnknownRole;

        auto& pa = ev.addParticle( role );
        pa.setId( p );
        pa.setPdgId( pdg_id, charge );
        pa.setStatus( (Particle::Status)pyjets_.k[0][p] );
        pa.setMomentum( Particle::Momentum( pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p] ) );
        pa.setMass( pyjets_.p[4][p] );
        if ( role != Particle::Role::UnknownRole ) {
          auto& moth = ev[pyjets_.k[2][p]-1];
          moth.setStatus( role == Particle::Role::CentralSystem
            ? Particle::Status::Resonance
            : Particle::Status::Fragmented );
          pa.addMother( moth );
        }
      }
      return true;
    }

    bool
    Pythia6Hadroniser::prepareHadronisation( Event& ev )
    {
      CG_DEBUG( "Pythia6Hadroniser" ) << "Hadronisation preparation called.";

      for ( const auto& part : ev.particles() ) {
        if ( part.status() != Particle::Status::Unfragmented )
          continue;
        //--- only loop over all protons to be fragmented

part.dump();

        const auto partons = pickPartonsContent();
        const double mx2 = part.mass2();
        const double mq = pymass( partons.first ), mq2 = mq*mq;
        const double mdq = pymass( partons.second ), mdq2 = mdq*mdq;

        //--- choose random direction in MX frame
        const double phi = 2.*M_PI*drand(), theta = acos( 2.*drand()-1. ); // theta angle

        //--- compute momentum of decay particles from MX
        const double px = std::sqrt( std::pow( mx2-mdq2+mq2, 2 ) / ( 4.*mx2 ) - mq2 );

        //--- build 4-vectors and boost decay particles
//        const auto pq = Particle::Momentum::fromPThetaPhi( px, theta, phi, std::hypot( px, mq ) );
        const auto pq = Particle::Momentum::fromPThetaPhi( px, theta, phi, part.mass() );
        CG_INFO("")<<pq << "|" << -pq;

        //--- singlet
        auto singl_mom = part.momentum();
        singl_mom.lorentzBoost( pq );

        auto& quark = ev.addParticle( part.role() );
        quark.addMother( ev[part.id()] );
        quark.setPdgId( partons.first );
        quark.setStatus( Particle::Status::FinalState );
        quark.setMomentum( singl_mom );

        //--- doublet
        auto doubl_mom = part.momentum();
        doubl_mom.lorentzBoost( -pq );

        auto& diquark = ev.addParticle( part.role() );
        diquark.addMother( ev[part.id()] );
        diquark.setPdgId( partons.second );
        diquark.setStatus( Particle::Status::FinalState );
        diquark.setMomentum( doubl_mom );

        std::cout << part.momentum()-(singl_mom+doubl_mom) << "|" << part.numDaughters() << std::endl;
        ev[part.id()].setStatus( Particle::Status::Fragmented );
      }
ev.dump();

      return true;
    }

    Pythia6Hadroniser::EventProperties
    Pythia6Hadroniser::fillParticles( const Event& ev ) const
    {
      pyjets_.n = 0;

      //--- initialising the string fragmentation variables
      EventProperties out;
      out.str_in_evt = 0;
      int jlpsf[MAX_STRING_EVENT][MAX_PART_STRING] = { 0 };

      for ( const auto& role : ev.roles() ) { // loop on roles
        unsigned int part_in_str = 0;
        bool role_has_string = false;
        for ( const auto& part : ev[role] ) {
          unsigned int np = part.id();

          pyjets_.p[0][np] = part.momentum().px();
          pyjets_.p[1][np] = part.momentum().py();
          pyjets_.p[2][np] = part.momentum().pz();
          pyjets_.p[3][np] = part.energy();
          pyjets_.p[4][np] = part.mass();
          pyjets_.k[0][np] = part.status() <= Particle::Status::Undefined
            ? 21 // incoming beam
            : (int)part.status();
          pyjets_.k[1][np] = part.integerPdgId();
          pyjets_.k[2][np] = part.mothers().empty()
            ? 0 // no mother
            : *( part.mothers().begin() )+1; // mother
          const auto& daug = part.daughters();
          if ( daug.empty() )
            pyjets_.k[3][np] = pyjets_.k[4][np] = 0; // no daughters
          else {
            pyjets_.k[3][np] = *daug.begin()+1; // daughter 1
            pyjets_.k[4][np] = *daug.rbegin()+1; // daughter 2
          }
          for ( int i = 0; i < 5; ++i )
            pyjets_.v[i][np] = 0.;

          if ( part.status() == Particle::Status::DebugResonance ) {
            pyjets_.k[0][np] = 1; //FIXME PYTHIA/JETSET workaround
            jlpsf[out.str_in_evt][part_in_str++] = part.id()+1;
            out.num_part_in_str[out.str_in_evt]++;
            role_has_string = true;
          }
          else if ( part.status() == Particle::Status::Undecayed )
            pyjets_.k[0][np] = 2; // intermediate resonance
          pyjets_.n++;
        }
        //--- at most one string per role
        if ( role_has_string )
          out.str_in_evt++;
      }

      //--- loop over the strings to bind everything together
      for ( unsigned short i = 0; i < out.str_in_evt; ++i ) {
        if ( out.num_part_in_str[i] < 2 )
          continue;

        std::ostringstream dbg;
        for ( unsigned short j = 0; j < out.num_part_in_str[i]; ++j )
          if ( jlpsf[i][j] != -1 )
            dbg << Form( "\n\t * %2d (pdgId=%4d)", jlpsf[i][j], pyjets_.k[1][jlpsf[i][j]-1] );

        CG_INFO( "Pythia6Hadroniser" )
          << "Joining " << out.num_part_in_str[i] << " particle" << utils::s( out.num_part_in_str[i] )
          << " with " << ev[jlpsf[i][0]].role() << " role"
          << " in a same string (id=" << i << ")"
          << dbg.str();

        pyjoin( out.num_part_in_str[i], jlpsf[i] );
      }
      return out;
    }

    std::pair<short,short>
    Pythia6Hadroniser::pickPartonsContent() const
    {
      const double ranudq = drand();
      if ( ranudq < 1./9. )
        return { PDG::down, 2203 }; // (d,uu1)
      if ( ranudq < 5./9. )
        return { PDG::up, 2101 };   // (u,ud0)
      return { PDG::up, 2103 };     // (u,ud1)
    }

  }
}

// register hadroniser and define alias
REGISTER_HADRONISER( pythia6, Pythia6Hadroniser )

