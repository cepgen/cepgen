#include "Pythia6Hadroniser.h"

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
Pythia6Hadroniser::Hadronise( Particle *part_ )
{
  pyjets_.p[0][0] = part_->GetMomentum().Px();
  pyjets_.p[1][0] = part_->GetMomentum().Py();
  pyjets_.p[2][0] = part_->GetMomentum().Pz();
  pyjets_.p[3][0] = part_->E();
  pyjets_.p[4][0] = part_->M();

  pyjets_.k[0][0] = 1; // status
  pyjets_.k[1][0] = 2; // particle id
  pyjets_.k[2][0] = 0; // mother
  pyjets_.k[3][0] = 0; // daughter 1
  pyjets_.k[4][0] = 0; // daughter 2

  this->pyexec();
  return true;
}

bool
Pythia6Hadroniser::Hadronise( Event *ev_ )
{
  Particle::Status status;

  ParticlesRef daug;
  ParticlesRef::iterator p;

  const unsigned int max_part_in_str = 3,
                     max_str_in_evt = 2;

  unsigned int num_part_in_str[max_str_in_evt];
  int jlrole[max_str_in_evt], jlpsf[max_str_in_evt][max_part_in_str];
  int criteria; //FIXME find an other name...
  
  try { PrepareHadronisation( ev_ ); } catch ( Exception& e ) { e.Dump(); throw e; }

  ParticleRoles rl = ev_->GetRoles();

  // First we initialise the string fragmentation variables
  for ( unsigned int i=0; i<max_str_in_evt; i++ ) {
    jlrole[i] = -1;
    num_part_in_str[i] = 0;
    for ( unsigned int j=0; j<max_part_in_str; j++ ) jlpsf[i][j] = -1;
  }
  
  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    Debugging( "Dump of the event before the hadronisation" );
    ev_->Dump();
  }
  
  // Filling the common block to propagate to PYTHIA6
  pyjets_.n = 0;
  unsigned int str_in_evt = 0;

  for ( ParticleRoles::iterator r=rl.begin(); r!=rl.end(); r++ ) {
    ParticlesRef pr = ev_->GetByRole( *r );
    unsigned int part_in_str = 0;
    for ( ParticlesRef::iterator part=pr.begin(); part!=pr.end(); part++ ) {
      Particle* p = *part;
      
      unsigned int np = p->id;
      
      pyjets_.p[0][np] = (double)p->GetMomentum().Px();
      pyjets_.p[1][np] = (double)p->GetMomentum().Py();
      pyjets_.p[2][np] = (double)p->GetMomentum().Pz();
      pyjets_.p[3][np] = (double)p->E();
      pyjets_.p[4][np] = (double)p->M();
      p->Dump();

      if ( p->status<=0 ) status = Particle::PythiaHIncoming;
      else status = p->status;
      pylist(2);
      pyjets_.k[0][np] = status;
      pyjets_.k[1][np] = (int)p->GetPDGId();
      
      //if ( p->GetMother()!=-1 ) pyjets_.k[2][np] = p->GetMother()+1; // mother
      if ( p->GetMothersIds().size()>0 ) pyjets_.k[2][np] = *( p->GetMothersIds().begin() )+1; // mother
      else pyjets_.k[2][np] = 0; // mother
      
      daug = ev_->GetDaughters( p );
      if ( daug.size()!=0 ) {
        pyjets_.k[3][np] = *p->GetDaughters().begin()+1; // daughter 1
        pyjets_.k[4][np] = *p->GetDaughters().end()+1; // daughter 2
      }
      else {
        pyjets_.k[3][np] = 0; // daughter 1
        pyjets_.k[4][np] = 0; // daughter 2
      }
      
      for (int i=0; i<5; i++) {
        pyjets_.v[i][np] = 0.;
      }
      
      if ( p->status==Particle::DebugResonance ) {
        pyjets_.k[0][np] = 1; //FIXME PYTHIA/JETSET workaround
        jlrole[str_in_evt] = p->role;
        jlpsf[str_in_evt][part_in_str] = p->id+1;
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
  if ( pyjets_.k[1][criteria]==2212 and pyjets_.k[0][criteria]==1 ) {
    //this->pylist( 2 );
    throw Exception( __PRETTY_FUNCTION__, "System is non-inelastic", JustWarning );
  }

  for ( unsigned int p=0; p<(unsigned int)pyjets_.n; p++ ) {

    //FIXME FIXME FIXME FIXME need to reimplement this first filter under this philosophy
    // First we filter the particles with status <= 0 :
    //  Status code = -1 : CLPAIR "internal" particles (not to be interacted with)
    //                 0 : Pythia6 empty lines
    //if (pyjets_.k[0][p]<=0) continue;
    //FIXME FIXME FIXME FIXME

    // We filter the first particles already present in the event
    if (p<oldnpart) continue;

    Particle pa;
    pa.id = p;
    pa.SetPDGId( static_cast<Particle::ParticleCode>( pyjets_.k[1][p] ) );
    if ( ev_->GetById( pyjets_.k[2][p]-1 )!=(Particle*)NULL ) {
      pa.role = ev_->GetById( pyjets_.k[2][p]-1 )->role; // Child particle inherits its mother's role
    }
    pa.status = static_cast<Particle::Status>( pyjets_.k[0][p] );
    pa.SetMomentum( Particle::Momentum( pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p] ) );
    pa.SetM( pyjets_.p[4][p] );
    pa.name = this->pyname( pa.GetPDGId() );
    pa.charge = (float)this->pyp( p+1, 6 );

    if ( pyjets_.k[2][p]!=0 ) {
      dbg << Form( "\n\t%2d (pdgId=%4d) has mother %2d (pdgId=%4d)", pa.id, pa.GetPDGId(), pyjets_.k[2][p], pyjets_.k[1][pyjets_.k[2][p]-1] );
      pa.SetMother( ev_->GetById( pyjets_.k[2][p]-1 ) );
    }

    ev_->AddParticle( pa );
  }
  Debugging( Form( "Passed the string construction stage.\n\t %d string objects were identified and constructed",
                   "%s", str_in_evt, dbg.str().c_str() ) );

  return true;
}

bool
Pythia6Hadroniser::PrepareHadronisation( Event *ev_ )
{
  Particle::ParticleCode singlet_id, doublet_id;
  double ranudq, ulmdq, ulmq;
  double ranmxp, ranmxt;
  double pmxp;
  double pmxda[4];
  double partpb[4];

  Debugging("Hadronisation preparation called!");

  ParticlesRef pp = ev_->GetParticles();
  for ( ParticlesRef::iterator part=pp.begin(); part!=pp.end(); part++ ) {

    Particle* p = *part;

    if ( p->status!=Particle::Undecayed ) continue;
    // One proton to be fragmented
    ranudq = drand();
    if ( ranudq<1./9. ) {
      singlet_id = Particle::dQuark;
      doublet_id = Particle::uu1Diquark;
    }
    else if ( ranudq<5./9. ) {
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
    ranmxp = 2.*Constants::Pi*drand();       // phi angle
    ranmxt = acos( 2.*drand()-1. ); // theta angle
    
    // Compute momentum of decay particles from MX
    pmxp = std::sqrt( std::pow( p->M2() - ulmdq*ulmdq + ulmq*ulmq, 2 ) / ( 4.*p->M2() ) - ulmq*ulmq );
    
    // Build 4-vectors and boost decay particles

    // Start with the singlet
    pmxda[0] = pmxp*sin( ranmxt )*cos( ranmxp );
    pmxda[1] = pmxp*sin( ranmxt )*sin( ranmxp );
    pmxda[2] = pmxp*cos( ranmxt );
    pmxda[3] = std::sqrt( pmxp*pmxp + ulmq*ulmq );
    
    Lorenb( p->M(), p->GetMomentum(), pmxda, partpb );
    
    if ( !( partpb[0]<0 ) and !( partpb[0]>0 ) ) return false;
    
    Particle singlet( p->role, singlet_id );
    singlet.status = Particle::DebugResonance;
    if ( !singlet.SetMomentum( partpb ) ) {
      throw Exception( __PRETTY_FUNCTION__, "ERROR while setting the 4-momentum of singlet", JustWarning );
    }
    //singlet.SetM(); //FIXME
    
    // Continue with the doublet
    pmxda[0] = -pmxda[0];
    pmxda[1] = -pmxda[1];
    pmxda[2] = -pmxda[2];
    pmxda[3] = std::sqrt( pmxp*pmxp + ulmdq*ulmdq );
    
    Lorenb( p->M(), p->GetMomentum(), pmxda, partpb );
    
    Particle doublet( p->role, doublet_id );
    doublet.status = Particle::DebugResonance;
    if ( !doublet.SetMomentum( partpb ) ) {
      throw Exception( __PRETTY_FUNCTION__, "ERROR while setting the 4-momentum of doublet", JustWarning );
    }
    //std::cout << "doublet, mass = " << doublet.M() << std::endl;
    //doublet.SetM(); //FIXME
    
    if ( p->NumDaughters()==0 ) {
      singlet.SetMother( ev_->GetById( p->id ) );
      doublet.SetMother( ev_->GetById( p->id ) );
      
      ev_->AddParticle( singlet );
      ev_->AddParticle( doublet );
      
      Debugging( "Quark/diquark content succesfully added to the event!" );
    }
    else { // Quark/diquark content already present in the event
      
      Debugging(Form("Quark/diquark content already present in the event!\n\tRole of these particles: %d", p->role));
      
      ParticlesIds daugh = p->GetDaughters();
      for ( ParticlesIds::const_iterator did=daugh.begin(); did!=daugh.end(); did++ ) {
        if ( ev_->GetById( *did )->GetPDGId()==Particle::uQuark
          or ev_->GetById( *did )->GetPDGId()==Particle::dQuark ) { // Quark
          singlet.SetMother( ev_->GetById( p->id ) );
          *( ev_->GetById( *did ) ) = singlet;
          Debugging( "Singlet replaced" );
        }
        else { // Diquark
          doublet.SetMother( ev_->GetById( p->id ) );
          *( ev_->GetById( *did ) ) = doublet;
          Debugging( "Doublet replaced" );
        }
      }
    }
  }
  return true;
}

#endif
