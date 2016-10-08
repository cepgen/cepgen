#ifndef OutputHandler_ExportHandler_h
#define OutputHandler_ExportHandler_h

#include "physics/Event.h"

namespace OutputHandler
{
  class ExportHandler
  {
   public:

    enum OutputType {
      HepMC, LHE
    };

   public:
    ExportHandler( const OutputType& type ) : fType( type ) {;}
    virtual ~ExportHandler() {;}
    
    void SetCrossSection( const float& xsec, const float& err_xsec ) {
      fCrossSect = xsec;
      fCrossSectErr = err_xsec;
    }
    void SetEventNumber( const unsigned int& ev_id ) { fEventNum = ev_id; }
    virtual void operator<<( const Event* ) {;}

   protected:

    OutputType fType;
    float fCrossSect, fCrossSectErr;
    unsigned int fEventNum;
    
  };
}

#endif
