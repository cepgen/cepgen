#include "LHEFHandler.h"

#ifdef HEPMC_VERSION3

namespace OutputHandler
{
  LHEFHandler::LHEFHandler( const char* filename ) :
    ExportHandler( ExportHandler::LHE ),
    lhe_output_( std::make_unique<LHEF::Writer>( filename ) )
  {
    //lhe_output_->headerBlock() << "";
    lhe_output_->initComments() << "Sample created using CepGen.\n";
  }

  void
  LHEFHandler::initialise( const Parameters& params )
  {
    lhe_output_->initComments() << " Input parameters:\n";
    params.dump( lhe_output_->initComments(), false );
    LHEF::HEPRUP run = lhe_output_->heprup;
    run.IDBMUP = std::pair<int,int>( params.in1pdg, params.in2pdg );
    run.EBMUP = std::pair<double,double>( params.in1p, params.in2p );
    run.NPRUP = 1;
    run.resize();
    run.XSECUP[0] = cross_sect_;
    run.XERRUP[0] = cross_sect_err_;
    run.XMAXUP[0] = 1.;
    run.LPRUP[0] = 1;
    lhe_output_->heprup = run;
    lhe_output_->init();
  }

  void
  LHEFHandler::operator<<( const Event* ev )
  {
    LHEF::HEPEUP out;
    out.heprup = &lhe_output_->heprup;
    out.XWGTUP = 1.;
    out.XPDWUP = std::pair<double,double>( 0., 0. );
    out.SCALUP = 0.;
    out.AQEDUP = Constants::alphaEM;
    out.AQCDUP = Constants::alphaQCD;
    out.NUP = ev->numParticles();
    out.resize();
    for ( unsigned short ip=0; ip<ev->numParticles(); ip++ ) {
      const Particle part = ev->getConstById( ip );
      out.IDUP[ip] = part.integerPdgId(); // PDG id
      out.ISTUP[ip] = part.status; // status code
      out.MOTHUP[ip] = std::pair<int,int>( *part.mothersIds().begin(), ( part.mothersIds().size()>0 ) ? *( part.mothersIds().end()-- ) : 0 ); // mothers
      out.ICOLUP[ip] = std::pair<int,int>( 0, 0 );
      out.PUP[ip] = std::vector<double>( { part.momentum().px(), part.momentum().py(), part.momentum().pz(), part.energy(), part.mass() } ); // momentum
      out.VTIMUP[ip] = 0.; // invariant lifetime
      out.SPINUP[ip] = 0.;
    }
    lhe_output_->eventComments() << "haha";
    lhe_output_->hepeup = out;
    lhe_output_->writeEvent();
  }
}

#endif
