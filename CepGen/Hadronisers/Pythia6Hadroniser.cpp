#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <algorithm>

/// Maximal number of characters to fetch for the particle's name
#define NAME_CHR 16

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
  /// Get kinematic information on a particle from the Pythia6 module
  extern double pyp_( int&, int& );
  /// Store one parton/particle in the PYJETS common block
  extern void py1ent_( int&, int&, double&, double&, double& );

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
        void readString( const char* param ) override;
        void init() override {}
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}
        //bool hadronise( const Particle* );

      private:
        inline static double pymass(int pdgid_) { return pymass_(pdgid_); }
        //inline static double pywidt(int pdgid_) { return pywidt_(pdgid_); }
        inline static void pyexec() { pyexec_(); }
        inline static void pyckbd() { pyckbd_(); }
        inline static void pygive( const std::string& line ) { pygive_( line.c_str(), line.length() ); }
        inline static void pylist( int mlist ) { pylist_( mlist ); }
        inline static double pyp( int role, int qty ) { return pyp_( role, qty ); }
        //inline static void py1ent( int* kf, double* pe, double theta, double phi ) { int one=1; py1ent_( &one, kf, pe, theta, phi ); }
        //inline static void py1ent( int* kf, double* pe, double theta, double phi ) { py1ent_( 1, kf, pe, theta, phi ); }
        inline static std::string pyname( int pdgid ) {
          char out[NAME_CHR];
          std::string s;
          pyname_( pdgid, out, NAME_CHR );
          s = std::string( out, NAME_CHR );
          //s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
          s.erase( remove( s.begin(), s.end(), ' ' ), s.end() );
          return s;
        }
        /**
         * \brief Connect entries with colour flow information
         * \param[in] njoin_ Number of particles to join in the colour flow
         * \param[in] ijoin_ List of particles unique identifier to join in the colour flow
         */
        inline static void pyjoin( int njoin, int ijoin[2] ) { return pyjoin_( njoin, *ijoin ); }
        bool prepareHadronisation( Event& );
    };

    Pythia6Hadroniser::Pythia6Hadroniser( const ParametersList& plist ) :
      GenericHadroniser( plist, "pythia6" )
    {
      //pygive( "MSTU(21)=1" );
    }

    void
    Pythia6Hadroniser::readString( const char* param )
    {
      pygive( param );
    }

    /*void
    Pythia6Hadroniser::run( Event& ev, double& weight, bool full )
    {
      pyjets_.p[0][0] = part->momentum().px();
      pyjets_.p[1][0] = part->momentum().py();
      pyjets_.p[2][0] = part->momentum().pz();
      pyjets_.p[3][0] = part->energy();
      pyjets_.p[4][0] = part->mass();

      pyjets_.k[0][0] = 1; // status
      pyjets_.k[1][0] = 2; // particle id
      pyjets_.k[2][0] = 0; // mother
      pyjets_.k[3][0] = 0; // daughter 1
      pyjets_.k[4][0] = 0; // daughter 2

      this->pyexec();
      return true;
    }*/

    bool
    Pythia6Hadroniser::run( Event& ev, double& weight, bool full )
    {
      Particles::iterator p;

      const unsigned int max_part_in_str = 3, max_str_in_evt = 2;

      unsigned int num_part_in_str[max_str_in_evt];
      int jlrole[max_str_in_evt], jlpsf[max_str_in_evt][max_part_in_str];
      int criteria; //FIXME find an other name...

      prepareHadronisation( ev );

      ParticleRoles rl = ev.roles();

      // First we initialise the string fragmentation variables
      for ( unsigned int i=0; i<max_str_in_evt; i++ ) {
        jlrole[i] = -1;
        num_part_in_str[i] = 0;
        for ( unsigned int j=0; j<max_part_in_str; j++ ) jlpsf[i][j] = -1;
      }

      if ( utils::Logger::get().level >= utils::Logger::Level::debug ) {
        CG_DEBUG( "Pythia6Hadroniser" ) << "Dump of the event before the hadronisation:";
        ev.dump();
      }

      // Filling the common block to propagate to PYTHIA6
      pyjets_.n = 0;
      unsigned int str_in_evt = 0;

      for ( ParticleRoles::iterator r=rl.begin(); r!=rl.end(); r++ ) {
        Particles& pr = ev[*r];
        unsigned int part_in_str = 0;
        for ( Particles::iterator part_it=pr.begin(); part_it!=pr.end(); ++part_it ) {
          Particle& part = *part_it;

          unsigned int np = part.id();

          pyjets_.p[0][np] = (double)part.momentum().px();
          pyjets_.p[1][np] = (double)part.momentum().py();
          pyjets_.p[2][np] = (double)part.momentum().pz();
          pyjets_.p[3][np] = (double)part.energy();
          pyjets_.p[4][np] = (double)part.mass();
          part.dump();

          pylist(2);
          pyjets_.k[0][np] = ( part.status() <= Particle::Status::Undefined )
            ? 21 // incoming beam
            : (int)part.status();
          pyjets_.k[1][np] = (int)part.pdgId();

          //if ( part.mother()!=-1 ) pyjets_.k[2][np] = part.mother()+1; // mother
          if ( part.mothers().size()>0 ) pyjets_.k[2][np] = *( part.mothers().begin() )+1; // mother
          else pyjets_.k[2][np] = 0; // mother

          const Particles daug = ev.daughters( part );
          if ( !daug.empty() ) {
            pyjets_.k[3][np] = *part.daughters().begin()+1; // daughter 1
            pyjets_.k[4][np] = *part.daughters().end()+1; // daughter 2
          }
          else {
            pyjets_.k[3][np] = 0; // daughter 1
            pyjets_.k[4][np] = 0; // daughter 2
          }

          for ( int i=0; i<5; i++ ) {
            pyjets_.v[i][np] = 0.;
          }

          if ( part.status() == Particle::Status::DebugResonance ) {
            pyjets_.k[0][np] = 1; //FIXME PYTHIA/JETSET workaround
            jlrole[str_in_evt] = part.role();
            jlpsf[str_in_evt][part_in_str] = part.id()+1;
            num_part_in_str[str_in_evt]++;
            part_in_str++;
          }
          pyjets_.n++;
        }
        if ( jlrole[str_in_evt]>0 ) {
          str_in_evt++;
        }
      }
      unsigned int oldnpart = pyjets_.n;

      std::ostringstream dbg;

      for ( unsigned int i=0; i<str_in_evt; i++ ) {
        if ( num_part_in_str[i]<2 ) continue;

        std::ostringstream os;
        for ( unsigned int j=0; j<num_part_in_str[i]; j++ ) {
          if ( jlpsf[i][j]!=-1 ) os << Form( "\n\t * %2d (pdgId=%4d)", jlpsf[i][j], pyjets_.k[1][jlpsf[i][j]-1] );
        }
        dbg << Form( "Joining %d particle(s) in a same string (%d) with role %d"
                     "%s", num_part_in_str[i], i, jlrole[i], os.str().c_str() );

        this->pyjoin( num_part_in_str[i], jlpsf[i] );
        this->pyexec();//FIXME FIXME FIXME
      }
      //this->pyexec();//FIXME FIXME FIXME
      this->pylist( 2 );

      criteria = oldnpart+1;
      for ( unsigned int i=0; i<max_str_in_evt; i++ ) criteria += num_part_in_str[i];
      this->pylist( 2 );
      if ( pyjets_.k[1][criteria] == 2212 && pyjets_.k[0][criteria] == 1 ) {
        //this->pylist( 2 );
        throw CG_WARNING( "Pythia6Hadroniser" ) << "System is non-inelastic.";
      }

      for ( unsigned int p=0; p<(unsigned int)pyjets_.n; p++ ) {

        //FIXME FIXME FIXME FIXME need to reimplement this first filter under this philosophy
        // First we filter the particles with status <= 0 :
        //  Status code = -1: CLPAIR "internal" particles (not to be interacted with)
        //                 0: Pythia6 empty lines
        //if (pyjets_.k[0][p]<=0) continue;
        //FIXME FIXME FIXME FIXME

        // We filter the first particles already present in the event
        if ( p < oldnpart ) continue;

        Particle pa;
        pa.setId( p );
        pa.setPdgId( (short)pyjets_.k[1][p] );
        if ( ev[pyjets_.k[2][p]-1].valid() ) {
          pa.setRole( ev[pyjets_.k[2][p]-1].role() ); // Child particle inherits its mother's role
        }
        pa.setStatus( static_cast<Particle::Status>( pyjets_.k[0][p] ) );
        pa.setMomentum( Particle::Momentum( pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p] ) );
        pa.setMass( pyjets_.p[4][p] );
        //pa.setCharge( (float)this->pyp( p+1, 6 ) );

        if ( pyjets_.k[2][p] != 0 ) {
          dbg << Form( "\n\t%2d (pdgId=%4d) has mother %2d (pdgId=%4d)", pa.id(), pa.pdgId(), pyjets_.k[2][p], pyjets_.k[1][pyjets_.k[2][p]-1] );
          pa.addMother( ev[pyjets_.k[2][p]-1] );
        }

        ev.addParticle( pa );
      }
      CG_DEBUG( "Pythia6Hadroniser" )
        << "Passed the string construction stage.\n\t"
        << " " << str_in_evt << " string objects were identified and constructed:"
        << dbg.str();

      return true;
    }

    bool
    Pythia6Hadroniser::prepareHadronisation( Event& ev )
    {
      short singlet_id, doublet_id;
      double ranudq, ulmdq, ulmq;
      double ranmxp, ranmxt;
      double pmxp;

      CG_DEBUG( "Pythia6Hadroniser" ) << "Hadronisation preparation called.";

      Particles pp = ev.particles();
      for ( Particles::iterator part_it=pp.begin(); part_it!=pp.end(); part_it++ ) {

        Particle& part = *part_it;

        if ( part.status() != Particle::Status::Undecayed ) continue;
        // One proton to be fragmented
        ranudq = drand();
        if ( ranudq < 1./9. ) {
          singlet_id = PDG::down;
          doublet_id = 2203; // uu1 diquark
        }
        else if ( ranudq < 5./9. ) {
          singlet_id = PDG::up;
          doublet_id = 2101; // ud0 diquark
        }
        else {
          singlet_id = PDG::up;
          doublet_id = 2103; // ud1 diquark
        }
        ulmdq = pymass( doublet_id );
        ulmq = pymass( singlet_id );

        // Choose random direction in MX frame
        ranmxp = 2.*M_PI*drand();       // phi angle
        ranmxt = acos( 2.*drand()-1. ); // theta angle

        // Compute momentum of decay particles from MX
        pmxp = std::sqrt( std::pow( part.mass2() - ulmdq*ulmdq + ulmq*ulmq, 2 ) / ( 4.*part.mass2() ) - ulmq*ulmq );

        // Build 4-vectors and boost decay particles

        // Start with the singlet
        Particle::Momentum pmxda( pmxp*sin( ranmxt )*cos( ranmxp ), pmxp*sin( ranmxt )*sin( ranmxp ), pmxp*cos( ranmxt ), sqrt( pmxp*pmxp + ulmq*ulmq ) );
        Particle::Momentum part_mom = part.momentum();
        part_mom.lorentzBoost( pmxda );

        if ( !( part_mom.px() < 0 ) && !( part_mom.px() > 0 ) ) return false;

        Particle singlet( part.role(), singlet_id );
        singlet.setStatus( Particle::Status::DebugResonance );
        singlet.setMomentum( part_mom );
        //singlet.setMass(); //FIXME

        // Continue with the doublet
        pmxda = Particle::Momentum( -pmxp*sin( ranmxt )*cos( ranmxp ), -pmxp*sin( ranmxt )*sin( ranmxp ), -pmxp*cos( ranmxt ), sqrt( pmxp*pmxp + ulmq*ulmq ) );
        part_mom = part.momentum();
        part_mom.lorentzBoost( pmxda );

        Particle doublet( part.role(), doublet_id );
        doublet.setStatus( Particle::Status::DebugResonance );
        doublet.setMomentum( part_mom );
        //std::cout << "doublet, mass = " << doublet.mass() << std::endl;
        //doublet.setMass(); //FIXME

        if ( part.numDaughters() == 0 ) {
          singlet.addMother( ev[part.id()] );
          doublet.addMother( ev[part.id()] );

          ev.addParticle( singlet );
          ev.addParticle( doublet );

          CG_DEBUG( "Pythia6Hadroniser" ) << "Quark/diquark content succesfully added to the event.";
        }
        else { // Quark/diquark content already present in the event
          CG_DEBUG( "Pythia6Hadroniser" )
            << "Quark/diquark content already present in the event!\n\t"
            << "Role of these particles: " << part.role() << ".";

          ParticlesIds daugh = part.daughters();
          for ( ParticlesIds::const_iterator did=daugh.begin(); did!=daugh.end(); did++ ) {
            if ( ev[*did].pdgId() == PDG::up
              || ev[*did].pdgId() == PDG::down ) { // Quark
              singlet.addMother( ev[part.id()] );
              ev[*did] = singlet;
              CG_DEBUG( "Pythia6Hadroniser" ) << "Singlet replaced.";
            }
            else { // Diquark
              doublet.addMother( ev[part.id()] );
              ev[*did] = doublet;
              CG_DEBUG( "Pythia6Hadroniser" ) << "Doublet replaced";
            }
          }
        }
      }
      return true;
    }
  }
}

// register hadroniser and define alias
REGISTER_HADRONISER( pythia6, Pythia6Hadroniser )

