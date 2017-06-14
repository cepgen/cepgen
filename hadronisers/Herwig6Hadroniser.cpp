#include "Herwig6Hadroniser.h"

#ifdef HERWIG6

Herwig6Hadroniser::Herwig6Hadroniser() : GenericHadroniser( "Herwig6" )
{
  //Debug( Form( "Constructor called" ) );
}

Herwig6Hadroniser::~Herwig6Hadroniser()
{
  //Debug( Form( "Destructor called" ) );
}

bool
Herwig6Hadroniser::Hadronise( Event* ev )
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  ev->Dump();
  
  const ParticlesRef pp = ev->GetParticles();
  unsigned int i=0;
  for ( ParticlesRef::const_iterator p=pp.begin(); p!=pp.end() && i<NMXHEP; p++, i++ ) {
    Particle* part = *p;

    if ( part->status==Particle::Undecayed )
      part->status = Particle::HerwigFragment; //FIXME workaround for cluster fragmentation
  
    hepevt_.idhep[i] = part->GetPDGId();
    hepevt_.isthep[i] = static_cast<int>( part->status );
    hepevt_.phep[i][0] = part->GetMomentum().Px();
    hepevt_.phep[i][1] = part->GetMomentum().Py();
    hepevt_.phep[i][2] = part->GetMomentum().Pz();
    hepevt_.phep[i][3] = part->E();
    hepevt_.phep[i][4] = part->M();
    for ( unsigned int j=0; j<4; j++ ) {
      hepevt_.vhep[i][j] = 0.; // Particle produced at the detector center
    }
    std::cout << "(" << i << ")--> " << part->GetPDGId() << std::endl;
  }
  hepevt_.nhep = i;
  
  this->hwdhad();
  
  std::cout << "after hadronisation" << std::endl;
  for ( unsigned int i=0; i<hepevt_.nhep; i++ ) {
    std::cout << "--> " << hepevt_.idhep[i] << std::endl;
  }
  
  return true;
}

#endif
