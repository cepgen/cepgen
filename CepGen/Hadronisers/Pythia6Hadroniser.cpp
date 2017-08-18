#include "Pythia6Hadroniser.h"

namespace CepGen
{
  namespace Hadroniser
  {

#ifdef PYTHIA6

    Pythia6Hadroniser::Pythia6Hadroniser() : GenericHadroniser( "Pythia6" )
    {
      //this->pygive("MSTU(21)=1");
    }

    Pythia6Hadroniser::~Pythia6Hadroniser()
    {
      //Debugging("Destructor called");
    }

    bool
    Pythia6Hadroniser::hadronise( const Particle* part )
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
    }

    bool
    Pythia6Hadroniser::hadronise( Event* ev )
    {
      Particle::Status status;

      Particles::iterator p;

      const unsigned int max_part_in_str = 3, max_str_in_evt = 2;

      unsigned int num_part_in_str[max_str_in_evt];
      int jlrole[max_str_in_evt], jlpsf[max_str_in_evt][max_part_in_str];
      int criteria; //FIXME find an other name...

      try { prepareHadronisation( ev ); } catch ( Exception& e ) { e.dump(); throw e; }

      ParticleRoles rl = ev->roles();

      // First we initialise the string fragmentation variables
      for ( unsigned int i=0; i<max_str_in_evt; i++ ) {
        jlrole[i] = -1;
        num_part_in_str[i] = 0;
        for ( unsigned int j=0; j<max_part_in_str; j++ ) jlpsf[i][j] = -1;
      }

      if ( Logger::get().level >= Logger::Debug ) {
        Debugging( "Dump of the event before the hadronisation" );
        ev->dump();
      }

      // Filling the common block to propagate to PYTHIA6
      pyjets_.n = 0;
      unsigned int str_in_evt = 0;

      for ( ParticleRoles::iterator r=rl.begin(); r!=rl.end(); r++ ) {
        Particles& pr = ev->getByRole( *r );
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

          if ( part.status() <= 0 ) status = Particle::PythiaHIncoming;
          else status = part.status();
          pylist(2);
          pyjets_.k[0][np] = status;
          pyjets_.k[1][np] = (int)part.pdgId();

          //if ( part.mother()!=-1 ) pyjets_.k[2][np] = part.mother()+1; // mother
          if ( part.mothersIds().size()>0 ) pyjets_.k[2][np] = *( part.mothersIds().begin() )+1; // mother
          else pyjets_.k[2][np] = 0; // mother

          const Particles daug = ev->daughters( part );
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

          if ( part.status() == Particle::DebugResonance ) {
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
        throw Exception( __PRETTY_FUNCTION__, "System is non-inelastic", JustWarning );
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
        pa.setPdgId( static_cast<Particle::ParticleCode>( pyjets_.k[1][p] ) );
        if ( ev->getById( pyjets_.k[2][p]-1 ).valid() ) {
          pa.setRole( ev->getById( pyjets_.k[2][p]-1 ).role() ); // Child particle inherits its mother's role
        }
        pa.setStatus( static_cast<Particle::Status>( pyjets_.k[0][p] ) );
        pa.setMomentum( Particle::Momentum( pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p] ) );
        pa.setMass( pyjets_.p[4][p] );
        pa.setCharge( (float)this->pyp( p+1, 6 ) );

        if ( pyjets_.k[2][p] != 0 ) {
          dbg << Form( "\n\t%2d (pdgId=%4d) has mother %2d (pdgId=%4d)", pa.id(), pa.pdgId(), pyjets_.k[2][p], pyjets_.k[1][pyjets_.k[2][p]-1] );
          pa.setMother( ev->getById( pyjets_.k[2][p]-1 ) );
        }

        ev->addParticle( pa );
      }
      Debugging( Form( "Passed the string construction stage.\n\t %d string objects were identified and constructed",
                       "%s", str_in_evt, dbg.str().c_str() ) );

      return true;
    }

    bool
    Pythia6Hadroniser::prepareHadronisation( Event* ev )
    {
      Particle::ParticleCode singlet_id, doublet_id;
      double ranudq, ulmdq, ulmq;
      double ranmxp, ranmxt;
      double pmxp;
      double pmxda[4];
      double partpb[4];

      Debugging("Hadronisation preparation called!");

      Particles pp = ev->particles();
      for ( Particles::iterator part_it=pp.begin(); part_it!=pp.end(); part_it++ ) {

        Particle& part = *part_it;

        if ( part.status() != Particle::Undecayed ) continue;
        // One proton to be fragmented
        ranudq = drand();
        if ( ranudq < 1./9. ) {
          singlet_id = Particle::dQuark;
          doublet_id = Particle::uu1Diquark;
        }
        else if ( ranudq < 5./9. ) {
          singlet_id = Particle::uQuark;
          doublet_id = Particle::ud0Diquark;
        }
        else {
          singlet_id = Particle::uQuark;
          doublet_id = Particle::ud1Diquark;
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
        pmxda[0] = pmxp*sin( ranmxt )*cos( ranmxp );
        pmxda[1] = pmxp*sin( ranmxt )*sin( ranmxp );
        pmxda[2] = pmxp*cos( ranmxt );
        pmxda[3] = std::sqrt( pmxp*pmxp + ulmq*ulmq );

        Lorenb( part.mass(), part.momentum(), pmxda, partpb );

        if ( !( partpb[0] < 0 ) && !( partpb[0] > 0 ) ) return false;

        Particle singlet( part.role(), singlet_id );
        singlet.setStatus( Particle::DebugResonance );
        singlet.setMomentum( partpb );
        //singlet.setMass(); //FIXME

        // Continue with the doublet
        pmxda[0] = -pmxda[0];
        pmxda[1] = -pmxda[1];
        pmxda[2] = -pmxda[2];
        pmxda[3] = std::sqrt( pmxp*pmxp + ulmdq*ulmdq );

        Lorenb( part.mass(), part.momentum(), pmxda, partpb );

        Particle doublet( part.role(), doublet_id );
        doublet.setStatus( Particle::DebugResonance );
        doublet.setMomentum( partpb );
        //std::cout << "doublet, mass = " << doublet.mass() << std::endl;
        //doublet.setMass(); //FIXME

        if ( part.numDaughters() == 0 ) {
          singlet.setMother( ev->getById( part.id() ) );
          doublet.setMother( ev->getById( part.id() ) );

          ev->addParticle( singlet );
          ev->addParticle( doublet );

          Debugging( "Quark/diquark content succesfully added to the event!" );
        }
        else { // Quark/diquark content already present in the event

          Debugging( Form( "Quark/diquark content already present in the event!\n\tRole of these particles: %d", part.role() ) );

          ParticlesIds daugh = part.daughters();
          for ( ParticlesIds::const_iterator did=daugh.begin(); did!=daugh.end(); did++ ) {
            if ( ev->getById( *did ).pdgId() == Particle::uQuark
              || ev->getById( *did ).pdgId() == Particle::dQuark ) { // Quark
              singlet.setMother( ev->getById( part.id() ) );
              ev->getById( *did ) = singlet;
              Debugging( "Singlet replaced" );
            }
            else { // Diquark
              doublet.setMother( ev->getById( part.id() ) );
              ev->getById( *did ) = doublet;
              Debugging( "Doublet replaced" );
            }
          }
        }
      }
      return true;
    }

#endif

  }
}
